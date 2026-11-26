import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from bin import subsample


class SubsampleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.input = self.directory / "simulation.tsv"
        self.columns = ["id", "parent_id", "death_date", "states", "seq", "note"]
        self.rows = [
            {
                "id": f"{cell_id:03d}",
                "parent_id": f"{parent:03d}" if parent else "",
                "death_date": str(death_date),
                "states": "edge" if cell_id % 2 else "center",
                "seq": "ACGT",
                "note": f"cell {cell_id}",
            }
            for cell_id, parent, death_date in (
                (1, None, 1), (2, 1, 2), (3, 1, 3), (4, 2, 3),
                (5, 2, 4), (6, 3, 5), (7, 3, 5), (8, 5, 5), (9, 5, 5),
            )
        ]
        with self.input.open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=self.columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)

    def run_cli(self, *, n=2, seed=None, name="sampled"):
        output = self.directory / name / "ancestry.tsv"
        endnodes = self.directory / name / "endnodes.tsv"
        command = [
            sys.executable, str(Path(subsample.__file__).resolve()),
            "--input", str(self.input), "--output", str(output),
            "--endnodes", str(endnodes), "--N", str(n),
        ]
        if seed is not None:
            command.extend(["--seed", str(seed)])
        result = subprocess.run(command, cwd=self.directory, capture_output=True, text=True, timeout=20)
        return result, output, endnodes

    def read_rows(self, path):
        with path.open(newline="") as file:
            reader = csv.DictReader(file, delimiter="\t")
            self.assertEqual(reader.fieldnames, self.columns)
            return list(reader)

    def test_cli_preserves_ancestry_schema_and_seeded_selection(self):
        original = self.input.read_bytes()
        result, output, endnodes = self.run_cli(seed=42)
        self.assertEqual(result.returncode, 0, result.stderr)
        ancestry = self.read_rows(output)
        survivors = self.read_rows(endnodes)
        self.assertEqual(len(survivors), 2)
        self.assertTrue(all(row["death_date"] == "5" for row in survivors))
        self.assertEqual(survivors, [row for row in ancestry if row["death_date"] == "5"])
        ids = {row["id"] for row in ancestry}
        self.assertIn("001", ids)
        self.assertEqual(len(ancestry), len(ids))
        self.assertTrue(all(not row["parent_id"] or row["parent_id"] in ids for row in ancestry))
        self.assertEqual(ancestry, [row for row in self.rows if row["id"] in ids])
        repeated, repeated_output, repeated_endnodes = self.run_cli(seed=42, name="repeated")
        self.assertEqual(repeated.returncode, 0, repeated.stderr)
        self.assertEqual(output.read_bytes(), repeated_output.read_bytes())
        self.assertEqual(endnodes.read_bytes(), repeated_endnodes.read_bytes())
        self.assertEqual(self.input.read_bytes(), original)

    def test_seed_is_optional_and_all_survivors_can_be_sampled(self):
        result, output, endnodes = self.run_cli(n=4)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({row["id"] for row in self.read_rows(endnodes)}, {"006", "007", "008", "009"})
        self.assertNotIn("004", {row["id"] for row in self.read_rows(output)})

    def test_invalid_counts_do_not_create_outputs(self):
        for n in (0, -1, 5, 1.5):
            with self.subTest(n=n):
                result, output, endnodes = self.run_cli(n=n)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertFalse(output.exists())
                self.assertFalse(endnodes.exists())

    def test_cli_rejects_overlapping_input_and_outputs(self):
        original = self.input.read_bytes()
        with self.assertRaises(SystemExit) as error:
            subsample.main([
                "--input", str(self.input), "--output", str(self.input),
                "--endnodes", str(self.directory / "endnodes.tsv"), "--N", "2",
            ])
        self.assertEqual(error.exception.code, 2)
        self.assertEqual(self.input.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
