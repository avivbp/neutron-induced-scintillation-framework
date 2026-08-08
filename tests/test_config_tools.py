#!/usr/bin/env python3
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.config_tools import (
    box_sensor_layout_macro_value,
    flatten_for_template,
    render_template,
)


class ConfigToolsTest(unittest.TestCase):
    def setUp(self):
        self.config = {
            "detector": {"model": "box_cryostat"},
            "geometry": {
                "active_lar": {"dimensions_cm": [100, 80, 120]},
                "fiducial": {"margin_cm": [5, 5, 5]},
                "cryostat": {"thickness_cm": 1.5},
                "neutron_detectors": [],
            },
            "sensors": {
                "sipm": {"tile_size_cm": 10},
                "pmt": {"side_length_cm": 20},
                "layouts": [
                    {
                        "id": "wall_sipms",
                        "type": "sipm",
                        "faces": ["+x", "-x"],
                        "pattern": "grid",
                        "active_size_cm": [10, 10],
                        "thickness_cm": 0.1,
                        "pitch_cm": [20, 20],
                        "edge_clearance_cm": 5,
                        "inset_cm": 0.01,
                    },
                    {
                        "id": "end_pmts",
                        "type": "pmt",
                        "faces": ["+z", "-z"],
                        "pattern": "grid",
                        "active_size_cm": [20, 20],
                        "thickness_cm": 0.1,
                        "pitch_cm": [30, 30],
                        "edge_clearance_cm": 10,
                        "inset_cm": 0.01,
                    },
                ],
            },
        }

    def test_box_layouts_render_into_run_macro(self):
        expected = (
            "wall_sipms:sipm:+x|-x:10:10:0.1:20:20:5:0.01;"
            "end_pmts:pmt:+z|-z:20:20:0.1:30:30:10:0.01"
        )
        self.assertEqual(box_sensor_layout_macro_value(self.config), expected)
        values = flatten_for_template(
            self.config,
            {"config_id": "box", "top_pmts": 0, "bottom_pmts": 0, "sipm_rows": 0},
        )
        macro = render_template(Path("macros/template_run.mac").read_text(), values)
        self.assertIn("/det/setRunLabel box", macro)
        self.assertIn("/det/setDetectorModel box_cryostat", macro)
        self.assertIn(f"/det/setBoxSensorLayouts {expected}", macro)
        self.assertNotIn("/det/setTOFWindowNs", macro)
        self.assertNotIn("{{", macro)
        self.assertAlmostEqual(values["active_area_cm2"], 9600.0)
        self.assertAlmostEqual(values["total_inner_surface_cm2"], 35200.0)
        self.assertAlmostEqual(values["coverage_percent"], 27.2727272727)

    def test_overlapping_pitch_is_rejected(self):
        self.config["sensors"]["layouts"][0]["pitch_cm"] = [5, 20]
        with self.assertRaisesRegex(ValueError, "pitch"):
            box_sensor_layout_macro_value(self.config)

    def test_legacy_config_disables_box_layouts(self):
        values = flatten_for_template({}, {})
        self.assertEqual(values["detector_model"], "nested_cell")
        self.assertEqual(values["box_sensor_layouts"], "none")
        self.assertEqual(
            values["tof_window_command"], "/det/setTOFWindowNs 40 50"
        )

    def test_empty_neutron_detector_list_renders_none(self):
        values = flatten_for_template(self.config, {})
        self.assertEqual(values["neutron_detector_config"], "none")


if __name__ == "__main__":
    unittest.main()
