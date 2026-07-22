import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_run_smr import apply_smr_fdr, map_candidate_probes, parse_smr_output, write_ma  # noqa: E402


class SmrInputTests(unittest.TestCase):
    def test_probe_mapping_strips_ensembl_version(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "test.epi"
            path.write_text("1\tENSG000001.4\t0\t100\tGENE1\t+\n", encoding="ascii")
            result = map_candidate_probes(
                path,
                [{"ensembl_gene_id": "ENSG000001", "gene_symbol": "GENE1"}],
            )
            self.assertEqual(result, {"ENSG000001": "ENSG000001.4"})

    def test_ma_uses_total_manifest_sample_size(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "outcome.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["snp", "effect_allele", "other_allele", "eaf", "beta", "se", "pvalue"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "snp": "rs1",
                        "effect_allele": "A",
                        "other_allele": "G",
                        "eaf": "0.2",
                        "beta": "0.1",
                        "se": "0.02",
                        "pvalue": "1e-6",
                    }
                )
            destination = root / "outcome.ma"
            write_ma(source, destination, 400)
            lines = destination.read_text(encoding="ascii").splitlines()
            self.assertEqual(lines[0], "SNP\tA1\tA2\tfreq\tb\tse\tp\tN")
            self.assertTrue(lines[1].endswith("\t400"))

    def test_parse_official_output_columns(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "test.smr"
            path.write_text(
                "ProbeID\ttopSNP\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HEIDI\tnsnp_HEIDI\n"
                "ENSG1\trs1\tA\tG\t0.2\t0.1\t0.02\t1e-6\t0.2\t0.03\t1e-10\t0.5\t0.1\t1e-5\t0.2\t12\n",
                encoding="ascii",
            )
            result = parse_smr_output(path)
            self.assertEqual(result["ENSG1"]["top_snp"], "rs1")
            self.assertEqual(result["ENSG1"]["p_heidi"], "0.2")

    def test_smr_fdr_is_adjusted_within_outcome(self):
        rows = [
            {"outcome_id": "A", "p_smr": "0.01"},
            {"outcome_id": "A", "p_smr": "0.02"},
            {"outcome_id": "B", "p_smr": "0.02"},
            {"outcome_id": "B", "p_smr": ""},
        ]
        apply_smr_fdr(rows)
        self.assertAlmostEqual(rows[0]["p_smr_fdr"], 0.02)
        self.assertAlmostEqual(rows[1]["p_smr_fdr"], 0.02)
        self.assertAlmostEqual(rows[2]["p_smr_fdr"], 0.02)
        self.assertEqual(rows[3]["p_smr_fdr"], "")


if __name__ == "__main__":
    unittest.main()
