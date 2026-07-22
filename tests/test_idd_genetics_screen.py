import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_genetics_screen import (  # noqa: E402
    apply_bh,
    choose_lead_instruments,
    read_finan_universe,
    read_reference_snps,
)


class DruggableScreenTests(unittest.TestCase):
    def test_finan_reader_excludes_missing_hgnc(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "finan.xlsx"
            frame = pd.DataFrame(
                [
                    {"ensembl_gene_id": "ENSG1", "hgnc_names": "GENE1", "druggability_tier": "Tier 1"},
                    {"ensembl_gene_id": "ENSG2", "hgnc_names": None, "druggability_tier": "Tier 2"},
                ]
            )
            with pd.ExcelWriter(path) as writer:
                frame.to_excel(writer, sheet_name="Data", index=False)
            result = read_finan_universe(path)
            self.assertEqual([row["ensembl_gene_id"] for row in result], ["ENSG1"])

    def test_lead_instrument_uses_smallest_pvalue(self):
        rows = [
            {"ensembl_gene_id": "ENSG1", "snp": "rs1", "pvalue": "1e-8", "beta": "0.1", "se": "0.01"},
            {"ensembl_gene_id": "ENSG1", "snp": "rs2", "pvalue": "1e-12", "beta": "0.2", "se": "0.01"},
        ]
        self.assertEqual(choose_lead_instruments(rows)["ENSG1"]["snp"], "rs2")

    def test_bh_adjustment_is_monotone_by_rank(self):
        rows = [{"pvalue": value} for value in (0.01, 0.02, 0.2, "")]
        apply_bh(rows, "pvalue", "fdr")
        self.assertAlmostEqual(rows[0]["fdr"], 0.03)
        self.assertAlmostEqual(rows[1]["fdr"], 0.03)
        self.assertAlmostEqual(rows[2]["fdr"], 0.2)
        self.assertNotIn("fdr", rows[3])

    def test_reference_membership_is_streamed(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "reference.bim"
            path.write_text("1 rs1 0 100 A G\n1 rs2 0 200 C T\n", encoding="ascii")
            self.assertEqual(read_reference_snps(path, {"rs2", "rs3"}), {"rs2"})


if __name__ == "__main__":
    unittest.main()
