#!/usr/bin/env python3

import unittest
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.neutron_scintillation_correction import apply_correction


CONFIG = {
    "mode": "simulation_native",
    "gamma_calibration": {
        "slope_photons_per_MeV": 100.0,
        "intercept_photons": 0.0,
    },
}


class NeutronScintillationCorrectionTests(unittest.TestCase):
    def test_dune_response_corrects_only_nuclear_recoil_pe(self):
        frame = pd.DataFrame(
            {
                "total_num_pe": [20.0, 4.0],
                "electronic_recoil_num_pe": [7.0, 4.0],
                "nuclear_recoil_num_pe": [10.0, 0.0],
                "other_num_pe": [3.0, 0.0],
                "nuclear_recoil_scintillation_photons": [50.0, 0.0],
                "nuclear_recoil_energy_deposit_MeV": [1.0, 0.0],
            }
        )

        result, _ = apply_correction(frame, CONFIG, "total_num_pe")

        self.assertEqual(result["numPE_electronic_recoil_raw"].tolist(), [7.0, 4.0])
        self.assertEqual(result["numPE_other_raw"].tolist(), [3.0, 0.0])
        self.assertEqual(result["numPE_nuclear_recoil_corrected"].tolist(), [20.0, 0.0])
        self.assertEqual(result["numPE_corrected"].tolist(), [30.0, 4.0])
        self.assertTrue(
            np.isnan(result.loc[1, "neutron_scintillation_correction_factor"])
        )

    def test_legacy_event_csv_keeps_whole_event_correction(self):
        frame = pd.DataFrame(
            {"numPE": [10.0], "numPhotons": [50.0], "eDep": [1.0]}
        )

        result, _ = apply_correction(frame, CONFIG, "numPE")

        self.assertEqual(result.loc[0, "numPE_corrected"], 20.0)


if __name__ == "__main__":
    unittest.main()
