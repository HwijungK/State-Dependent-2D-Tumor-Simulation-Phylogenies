import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import numpy as np

from bin import sim
from bin.create_beast_xml import create_beast_xml
from bin.log_to_fasta import convert_fasta
from bin.subsample import subsample_log


class BeastXmlTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.log = self.directory / "simulation.tsv"
        with patch.object(sim, "rng", np.random.default_rng(41)):
            sim.run_simulation(size=20, target_cells=64, birth_rate=1, death_rate=0, log_file=self.log)
        sampled = subsample_log(self.log, n=12, seed=42)
        self.survivors = sampled.loc[sampled["death_date"].eq(sampled["death_date"].max())]
        survivors_file = self.directory / "survivors.tsv"
        self.survivors.to_csv(survivors_file, sep="\t", index=False)
        self.fasta = convert_fasta(survivors_file, self.directory / "sampled.fasta")

    def make_xml(self, model):
        output = create_beast_xml(
            self.fasta, simulation_log=self.log, output_path=self.directory / f"{model}.xml",
            clock_model=model, chain_length=100, log_every=20,
        )
        return ET.parse(output).getroot()

    def test_alignment_clock_calibration_and_file_names(self):
        root = self.make_xml("sdevo")
        self.assertEqual(len(root.findall("data/sequence")), len(self.survivors))
        self.assertEqual(
            {sequence.get("value") for sequence in root.findall("data/sequence")},
            set(self.survivors["seq"]),
        )
        clock = root.find(".//branchRateModel")
        self.assertEqual(clock.get("spec"), "sdevo.bdmm.BdmmStateClock")
        self.assertIsNotNone(clock.find("typedTree"))
        global_rate = clock.find("parameter[@name='clock.rate']")
        self.assertEqual(global_rate.get("estimate"), "false")
        self.assertGreater(float(global_rate.text), 0)
        self.assertEqual(root.find("run").get("chainLength"), "100")
        self.assertTrue(all(logger.get("logEvery") == "20" for logger in root.iter("logger")))
        filenames = {logger.get("fileName") for logger in root.iter("logger") if logger.get("fileName")}
        self.assertEqual(filenames, {"sdevo.log", "sdevo.trees", "sdevo.typed.trees", "sdevo.typed.node.trees", "sdevo.traj"})

    def test_strict_clock_preserves_the_bdmm_prior_and_typed_tree_logging(self):
        sdevo = self.make_xml("sdevo")
        strict = self.make_xml("strict")
        clock = strict.find(".//branchRateModel")
        self.assertEqual(clock.get("spec"), "beast.evolution.branchratemodel.StrictClockModel")
        self.assertIsNone(clock.get("relativeClockRates"))
        self.assertIsNone(clock.get("normalize"))
        self.assertIsNone(clock.find("typedTree"))
        self.assertEqual(
            clock.find("parameter[@name='clock.rate']").text,
            sdevo.find(".//branchRateModel/parameter[@name='clock.rate']").text,
        )
        self.assertIsNotNone(strict.find(".//distribution[@spec='bdmmprime.distribution.BirthDeathMigrationDistribution']"))
        self.assertIsNotNone(strict.find(".//logger[@fileName='strict.typed.node.trees']"))
        for root in (sdevo, strict):
            values = root.find(".//rhoSampling/parameter[@name='values']").text
            self.assertEqual(len(values.split()), 2)


if __name__ == "__main__":
    unittest.main()
