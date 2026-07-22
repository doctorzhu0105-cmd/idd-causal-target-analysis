import importlib.util
import tempfile
import unittest
from pathlib import Path

import openpyxl


MODULE_PATH = Path(__file__).resolve().parents[1] / "py" / "mr_prepare_druggable_genes.py"
SPEC = importlib.util.spec_from_file_location("mr_prepare_druggable_genes", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class MrPrepareDruggableGenesTest(unittest.TestCase):
    def test_expected_columns_match_finan_table(self):
        self.assertEqual(len(MODULE.EXPECTED_COLUMNS), 12)
        self.assertEqual(MODULE.EXPECTED_COLUMNS[0], "ensembl_gene_id")
        self.assertEqual(MODULE.EXPECTED_COLUMNS[1], "druggability_tier")

    def test_sha256_file_is_deterministic(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "source.bin"
            path.write_bytes(b"druggable-gene-test\n")
            self.assertEqual(
                MODULE.sha256_file(path),
                "b5fb96fedbb91d4db620732625f1d1d7abcdc9d213af61d5c3459fc5b842db24",
            )

    def test_fixture_header_can_be_read_by_openpyxl(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "fixture.xlsx"
            workbook = openpyxl.Workbook()
            sheet = workbook.active
            sheet.title = "Data"
            sheet.append(MODULE.EXPECTED_COLUMNS)
            workbook.create_sheet("Column Descriptions")
            workbook.create_sheet("Legend")
            workbook.save(path)
            loaded = openpyxl.load_workbook(path, read_only=True, data_only=True)
            header = list(next(loaded["Data"].iter_rows(values_only=True)))
            loaded.close()
            self.assertEqual(header, MODULE.EXPECTED_COLUMNS)

    def test_druggability_flag_columns_are_fixed(self):
        self.assertEqual(MODULE.EXPECTED_COLUMNS[9:], [
            "small_mol_druggable",
            "bio_druggable",
            "adme_gene",
        ])


if __name__ == "__main__":
    unittest.main()
