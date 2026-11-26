import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from bin import sim


class SimulationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        rng_patch = patch.object(sim, "rng", np.random.default_rng(41))
        rng_patch.start()
        self.addCleanup(rng_patch.stop)

    def read_log(self, path):
        with path.open(newline="") as file:
            return list(csv.DictReader(file, delimiter="\t"))

    def test_reaches_target_and_finalizes_without_display_or_delays(self):
        for model in ("boundary", "unrestricted"):
            with self.subTest(model=model):
                output = io.StringIO()
                log_file = self.directory / f"{model}.tsv"
                with (
                    contextlib.redirect_stdout(output),
                    patch.object(sim, "visualize_world", side_effect=AssertionError("visualization")),
                    patch("time.sleep", side_effect=AssertionError("sleep")),
                ):
                    world = sim.run_simulation(
                        size=6, target_cells=8, growth_model=model,
                        birth_rate=1, death_rate=0, log_file=log_file,
                    )
                self.assertEqual(output.getvalue(), "")
                rows = self.read_log(log_file)
                survivors = [row for row in rows if int(row["death_date"]) == world.step + 1]
                self.assertEqual(len(survivors), 8)
                self.assertEqual(len(rows), 15)
                self.assertTrue(all(row["seq"] and set(row["seq"]) <= set("ACGT") for row in rows))
                self.assertTrue(world.update_log_file.is_file())
                self.assertEqual(world.cell_count, 0)

    def test_restarts_after_extinction_and_keeps_only_successful_attempt(self):
        log_file = self.directory / "extinction.tsv"
        with patch.object(sim.Cell, "choose_action", side_effect=["birth", "death", "death", "birth", "birth"]):
            world = sim.run_simulation(
                size=3, target_cells=3, birth_rate=1, death_rate=0, log_file=log_file,
            )
        self.assertEqual(world.step, 2)
        rows = self.read_log(log_file)
        self.assertEqual(len(rows), 5)
        self.assertEqual({row["id"] for row in rows}, {"1", "2", "3", "4", "5"})
        self.assertEqual(sum(int(row["death_date"]) == world.step + 1 for row in rows), 3)
        self.assertTrue(all(row["seq"] for row in rows))
        updates = self.read_log(world.update_log_file)
        self.assertEqual(len(updates), 2)
        self.assertTrue(all(row["action"] == "birth" and row["divided"] == "1" for row in updates))

    def test_cli_runs_only_the_selected_config_object(self):
        config = self.directory / "runs.json"
        config.write_text(json.dumps([
            {"target_cells": 2, "log_file": "unused.tsv", "birth_rate": 1, "death_rate": 0},
            {
                "size": 5, "target_cells": 5, "growth_model": "unrestricted",
                "log_file": "selected.tsv", "birth_rate": 1, "death_rate": 0,
            },
        ]))
        result = subprocess.run(
            [sys.executable, str(Path(sim.__file__).resolve()),
             "--config", str(config), "--run-index", "1"],
            cwd=self.directory, capture_output=True, text=True, timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.directory / "unused.tsv").exists())
        rows = self.read_log(self.directory / "selected.tsv")
        last_death = max(int(row["death_date"]) for row in rows)
        self.assertEqual(sum(int(row["death_date"]) == last_death for row in rows), 5)
        self.assertTrue((self.directory / "selected_updates.tsv").is_file())
        self.assertEqual(len(result.stdout.splitlines()), 2)

    def test_cli_rejects_invalid_config_without_creating_logs(self):
        config = self.directory / "runs.json"
        for parameters, index in (
            ({}, 0),
            ([], 0),
            ([{}], -1),
            ([{}], 1),
            ([None], 0),
            ([{"log_file": str(self.directory / "invalid.tsv"), "birth_rate": 2}], 0),
        ):
            with self.subTest(parameters=parameters, index=index):
                config.write_text(json.dumps(parameters))
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as error:
                        sim.main(["--config", str(config), "--run-index", str(index)])
                self.assertEqual(error.exception.code, 2)
                self.assertFalse((self.directory / "invalid.tsv").exists())


if __name__ == "__main__":
    unittest.main()
