"""Focused conversion checks for the RFD4 pathway."""
import json
import tempfile
import unittest
from pathlib import Path

from Bio.PDB import Chain, MMCIFIO, MMCIFParser, Model, PDBParser, Structure

from analyse_rfd4 import convert
from prepare_rfd4 import prepare
from use_rfd4_sequences import run


ROOT = Path(__file__).resolve().parents[2]


class RFD4ConversionTests(unittest.TestCase):
    def test_target_hotspots_and_generated_sequence_enter_refolding(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            prepare(ROOT / "benchmarkdata/5o45_pd-l1.pdb", "60-100", "A56,A115", work)
            spec = json.loads((work / "binder.json").read_text())["binder"]
            self.assertEqual(spec["contig"], "60-100, /0, A1-115")
            self.assertEqual([entry["select"] for entry in spec["conditions"]["C_HOT"]["true"]],
                             ["sel('A/*/40')", "sel('A/*/99')"])

            target = next(MMCIFParser(QUIET=True).get_structure("target", str(work / "target.cif")).get_chains())
            sample = Structure.Structure("sample")
            model = Model.Model(0)
            sample.add(model)
            # Deliberately place the target first: converter must identify it by sequence.
            target_copy = target.copy()
            target_copy.id = "Z"
            model.add(target_copy)
            binder = Chain.Chain("Q")
            for index, residue in enumerate(list(target)[:8], 1):
                item = residue.copy()
                item.id = (" ", index, " ")
                item.resname = "ALA"
                binder.add(item)
            model.add(binder)
            writer = MMCIFIO()
            writer.set_structure(sample)
            samples = work / "inference_samples" / "binder_0"
            samples.mkdir(parents=True)
            writer.save(str(samples / "binder_0_generated.cif"))

            folds = work / "folds"
            convert(work / "inference_samples", work / "target.cif", 4, folds)
            output = PDBParser(QUIET=True).get_structure("fold", str(folds / "fold_4.pdb"))
            self.assertEqual([chain.id for chain in output[0]], ["A", "B"])
            self.assertEqual(len(list(output[0]["A"].get_residues())), 8)
            self.assertEqual(len(list(output[0]["B"].get_residues())), 115)

            sequences = work / "sequences"
            run(folds, sequences, None, None, None, None)
            self.assertTrue((sequences / "fold_4_seq_0.pdb").exists())
            record = json.loads((sequences / "rfd4_seq_metadata.jsonl").read_text())
            self.assertEqual((record["fold_id"], record["seq_id"], record["seq_length"]), (4, 0, 8))

            rejected = work / "rejected"
            run(folds, rejected, None, None, None, 2.0)
            self.assertFalse((rejected / "fold_4_seq_0.pdb").exists())
            self.assertTrue((rejected / "rfd4_seq_metadata.jsonl").exists())

    def test_missing_hotspot_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "absent"):
                prepare(ROOT / "benchmarkdata/5o45_pd-l1.pdb", "60", "A999", tmp)

    def test_flexible_residues_remove_only_coordinate_conditioning(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            # A56 is deliberately both flexible and a hotspot: interaction conditioning and
            # coordinate conditioning are independent and should coexist.
            prepare(ROOT / "benchmarkdata/5o45_pd-l1.pdb", "60", "A56", work, "A17-19,A21,A56")
            conditions = json.loads((work / "binder.json").read_text())["binder"]["conditions"]

            self.assertEqual(conditions["C_SEQ"], {
                "true": [{"select": "sel('A')"}],
            })
            self.assertEqual(
                conditions["C_CRDGRP"]["null"],
                [
                    {"select": "sel('A/*/1/[N,CA,C,O]')"},
                    {"select": "sel('A/*/2/[N,CA,C,O]')"},
                    {"select": "sel('A/*/3/[N,CA,C,O]')"},
                    {"select": "sel('A/*/5/[N,CA,C,O]')"},
                    {"select": "sel('A/*/40/[N,CA,C,O]')"},
                ],
            )
            self.assertEqual(conditions["C_HOT"]["true"], [{"select": "sel('A/*/40')"}])

    def test_whole_target_can_be_structurally_flexible(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            prepare(ROOT / "benchmarkdata/5o45_pd-l1.pdb", "60", "", work, "A")
            conditions = json.loads((work / "binder.json").read_text())["binder"]["conditions"]
            self.assertEqual(len(conditions["C_CRDGRP"]["null"]), 115)
            self.assertNotIn("C_HOT", conditions)

    def test_invalid_flexible_residues_rejected(self):
        for flexible, message in (("B17", "target chain A"), ("A999", "absent"), ("A20-17", "range")):
            with self.subTest(flexible=flexible), tempfile.TemporaryDirectory() as tmp:
                with self.assertRaisesRegex(ValueError, message):
                    prepare(ROOT / "benchmarkdata/5o45_pd-l1.pdb", "60", "", tmp, flexible)


if __name__ == "__main__":
    unittest.main()
