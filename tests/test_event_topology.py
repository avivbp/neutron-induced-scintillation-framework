import json
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None

if pd is not None:
    from scripts.event_topology import (
        classify_event_topologies,
        interaction_truth_path_for_pe,
        true_single_elastic_mask_from_truth,
    )


@unittest.skipIf(pd is None, "pandas is not installed for this Python interpreter")
class EventTopologyTest(unittest.TestCase):
    def setUp(self):
        self.interactions = pd.DataFrame(
            [
                (1, "elastic", "active_lar|fiducial_lar", 0),
                (2, "elastic", "active_lar|fiducial_lar", 0),
                (2, "elastic", "active_lar", 0),
                (3, "inelastic", "active_lar|fiducial_lar", 0),
                (4, "elastic", "active_lar|fiducial_lar", 0),
                (4, "inelastic", "active_lar|fiducial_lar", 0),
                (4, "elastic", "active_lar", 1),
                (5, "capture", "active_lar", 1),
                (6, "fission", "active_lar", 0),
                (7, "other_hadronic", "active_lar", 0),
                (8, "mystery", "active_lar", 0),
                (9, "transport_exit", "", 0),
            ],
            columns=("event_id", "channel", "volume_roles", "parent_id"),
        )

    def test_channel_topologies_and_counts(self):
        result = classify_event_topologies(
            self.interactions, event_ids=list(range(1, 11))
        ).set_index("event_id")
        expected = {
            1: "single_elastic",
            2: "multiple_elastic",
            3: "inelastic",
            4: "elastic_plus_inelastic",
            5: "capture",
            6: "fission",
            7: "other_hadronic",
            8: "unclassified",
            9: "transport_exit_only",
            10: "no_recorded_interaction",
        }
        self.assertEqual(result["event_type"].to_dict(), expected)
        self.assertEqual(result.loc[4, "elastic_count"], 2)
        self.assertEqual(result.loc[4, "inelastic_count"], 1)
        self.assertEqual(result.loc[4, "secondary_neutron_interaction_count"], 1)
        self.assertEqual(result.loc[1, "fiducial_elastic_count"], 1)

    def test_truth_file_resolution_and_true_single_mask(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            pe_path = run_dir / "numPE_2_topPMTs_1_botPMTs_0_SiPMRows.csv"
            pe_path.touch()
            metadata = [
                {
                    "config_id": "C 01/test",
                    "top_pmts": 2,
                    "bottom_pmts": 1,
                    "sipm_rows": 0,
                }
            ]
            (run_dir / "scan_metadata.json").write_text(json.dumps(metadata))
            truth_path = run_dir / "neutron_interactions_C_01_test.csv"
            self.interactions.to_csv(truth_path, index=False)

            self.assertEqual(interaction_truth_path_for_pe(pe_path), truth_path)
            events = pd.DataFrame({"eventID": [1, 2, 3, 10]})
            self.assertEqual(
                true_single_elastic_mask_from_truth(events, truth_path).tolist(),
                [True, False, False, False],
            )


if __name__ == "__main__":
    unittest.main()
