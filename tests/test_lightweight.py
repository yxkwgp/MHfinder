import io
import os
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


class ImportSmokeTests(unittest.TestCase):
    def test_core_imports(self):
        import core.snp_filter  # noqa: F401
        import core.prep_MHs  # noqa: F401
        import core.scgwas_screen  # noqa: F401
        import cli  # noqa: F401

    def test_cli_help_exits_successfully(self):
        import cli

        commands = [cli.snp_filter_main, cli.prep_mhs_main, cli.screen_mhs_main]
        old_argv = sys.argv[:]
        try:
            for command in commands:
                sys.argv = ["cmd", "--help"]
                with self.assertRaises(SystemExit) as ctx, redirect_stdout(io.StringIO()):
                    command()
                self.assertEqual(ctx.exception.code, 0)
        finally:
            sys.argv = old_argv


class TableParsingTests(unittest.TestCase):
    def test_read_tsv_table(self):
        from core.utils import read_table

        path = ROOT / "examples" / "toy" / "sample_info.tsv"
        df = read_table(str(path))
        self.assertEqual(list(df.columns), ["FamilyID", "SampleID", "Population"])
        self.assertEqual(df.shape[0], 6)
        self.assertEqual(sorted(df["Population"].unique().tolist()), ["AFR", "EAS", "EUR"])


class MicrohaplotypeParserTests(unittest.TestCase):
    def test_candidate_mh_format(self):
        path = ROOT / "examples" / "toy" / "candidate_MHs.txt"
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(lines[0], "chr|SNP1-SNP2...SNPn")
        records = []
        for line in lines[1:]:
            chrom, snps = line.split("|")
            records.append((chrom, snps.split("-")))
        self.assertEqual(len(records), 3)
        self.assertTrue(all(len(snps) == 3 for _, snps in records))

    def test_selected_mh_format_has_three_line_blocks(self):
        path = ROOT / "examples" / "toy" / "selected_MHs_example.txt"
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(lines) % 3, 0)
        for i in range(0, len(lines), 3):
            self.assertTrue(lines[i].startswith(">"))
            self.assertTrue(lines[i + 1].startswith("*"))
            self.assertGreaterEqual(len(lines[i + 2].split("\t")), 3)


if __name__ == "__main__":
    unittest.main()
