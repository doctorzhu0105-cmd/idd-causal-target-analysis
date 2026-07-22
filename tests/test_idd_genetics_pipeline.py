import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_genetics_pipeline import (  # noqa: E402
    _leave_one_out_results,
    _sensitivity_results,
    _enrich_eqtl,
    _extract_af,
    _extract_eqtl,
    _extract_outcome,
    load_rows,
    run_clump,
    split_rsids,
)


EQTL_HEADER = (
    "Pvalue\tSNP\tSNPChr\tSNPPos\tAssessedAllele\tOtherAllele\tZscore\tGene\t"
    "GeneSymbol\tGeneChr\tGenePos\tNrCohorts\tNrSamples\tFDR\tBonferroniP\n"
)


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_candidate_eqtl_af_and_effect_extraction(self):
        source = self.root / "eqtl.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(EQTL_HEADER)
            handle.write("1e-10\trs1\t1\t110\tA\tG\t6\tENSG1\tGENE1\t1\t100\t30\t30000\t0.01\t1e-8\n")
            handle.write("1e-12\trs2\t1\t5000000\tA\tC\t7\tENSG1\tGENE1\t1\t100\t30\t30000\t0.01\t1e-9\n")
            handle.write("1e-20\trs3\t1\t110\tA\tC\t8\tENSG2\tGENE2\t1\t100\t30\t30000\t0.01\t1e-10\n")
        raw = self.root / "candidate.raw.tsv.gz"
        _extract_eqtl(
            source,
            raw,
            {"ENSG1": "GENE1"},
            pvalue_max=5e-8,
            f_min_exclusive=10,
            window_bp=1_000_000,
        )
        rows = load_rows(raw)
        self.assertEqual([row["snp"] for row in rows], ["rs1"])

        af_source = self.root / "af.gz"
        with gzip.open(af_source, "wt", encoding="utf-8") as handle:
            handle.write("SNP\thg19_chr\thg19_pos\tAlleleA\tAlleleB\tallA_total\tallAB_total\tallB_total\tAlleleB_all\n")
            handle.write("rs1\t1\t110\tG\tA\t1\t1\t1\t0.25\n")
        af = _extract_af(af_source, {"rs1"}, self.root / "af_extract.tsv.gz")
        enriched = self.root / "candidate.tsv.gz"
        _enrich_eqtl(raw, enriched, af)
        result = load_rows(enriched)[0]
        self.assertEqual(result["af_status"], "matched_allele_b")
        self.assertAlmostEqual(float(result["eaf"]), 0.25)
        self.assertAlmostEqual(float(result["beta"]) / float(result["se"]), 6.0)

    def test_finngen_rsid_intersection(self):
        source = self.root / "outcome.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write("#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls\n")
            handle.write("1\t100\tG\tA\trs1,rs9\tGENE1\t0.01\t2\t0.1\t0.02\t0.2\t0.2\t0.2\n")
        destination = self.root / "outcome.tsv.gz"
        _extract_outcome(source, destination, "TEST", {"rs1", "rs2"})
        rows = load_rows(destination)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["snp"], "rs1")
        self.assertEqual(rows[0]["effect_allele"], "A")

    def test_split_rsids(self):
        self.assertEqual(split_rsids("rs1,rs2;rs3|rs4"), {"rs1", "rs2", "rs3", "rs4"})


class ClumpTests(unittest.TestCase):
    def test_missing_reference_is_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = SimpleNamespace(
                plink=root / "plink.exe",
                ld_prefix=root / "reference",
                instrument_p=5e-8,
                clump_r2=0.001,
                clump_kb=10000,
            )
            selected, status = run_clump(
                [{"snp": "rs1", "pvalue": "1e-10"}],
                p_column="pvalue",
                label="test",
                args=args,
                clump_dir=root,
            )
            self.assertEqual(selected, set())
            self.assertEqual(status, "plink_missing")


class SensitivityTests(unittest.TestCase):
    def setUp(self):
        self.candidate = {"gene_symbol": "GENE1", "ensembl_gene_id": "ENSG1"}
        self.rows = [
            {
                "snp": f"rs{index}",
                "beta_exposure": bx,
                "se_exposure": 0.005,
                "beta_outcome": bx * 0.5,
                "se_outcome": 0.01,
            }
            for index, bx in enumerate((0.1, 0.2, 0.3), start=1)
        ]

    def test_sensitivity_methods_are_estimable_with_three_instruments(self):
        result = _sensitivity_results(self.candidate, "TEST", "expression_to_disease", self.rows)
        self.assertEqual({row["method"] for row in result}, {"weighted_median", "mr_egger"})
        self.assertTrue(all(row["status"] == "estimable" for row in result))

    def test_leave_one_out_has_one_row_per_exclusion(self):
        result = _leave_one_out_results(self.candidate, "TEST", "expression_to_disease", self.rows)
        self.assertEqual(len(result), 3)
        self.assertEqual({row["excluded_snp"] for row in result}, {"rs1", "rs2", "rs3"})


if __name__ == "__main__":
    unittest.main()
