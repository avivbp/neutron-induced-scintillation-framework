#!/usr/bin/env python3
"""Small helpers shared by the LAr scintillation scan scripts.

The scripts intentionally avoid non-standard dependencies except PyYAML,
matplotlib, and pandas for analysis.  If PyYAML is unavailable, install with:
    python -m pip install pyyaml
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML is required. Install it with: python -m pip install pyyaml"
    ) from exc


Number = int | float


def read_yaml(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration file {path} did not contain a YAML mapping.")
    return data


def read_scan_csv(path: str | Path) -> List[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"Sensor scan file {path} is empty.")
    for row in rows:
        for key in ("top_pmts", "bottom_pmts", "sipm_rows"):
            if key in row:
                row[key] = int(row[key])
    return rows


def deep_get(data: Mapping[str, Any], dotted_key: str, default: Any = None) -> Any:
    current: Any = data
    for part in dotted_key.split("."):
        if not isinstance(current, Mapping) or part not in current:
            return default
        current = current[part]
    return current


def neutron_detector_macro_value(config: Mapping[str, Any]) -> str:
    detectors = deep_get(config, "geometry.neutron_detectors", None)
    if detectors is None:
        detectors = [
            {"label": "A0", "angle_deg": 25.0, "distance_cm": 100.0},
            {"label": "A1", "angle_deg": 40.0, "distance_cm": 100.0},
            {"label": "A2", "angle_deg": 50.0, "distance_cm": 100.0},
            {"label": "A3", "angle_deg": 60.0, "distance_cm": 100.0},
            {"label": "A4", "angle_deg": 90.0, "distance_cm": 100.0},
        ]
    if not isinstance(detectors, list):
        raise ValueError("geometry.neutron_detectors must be a list.")

    entries: list[str] = []
    for index, detector in enumerate(detectors):
        if not isinstance(detector, Mapping):
            raise ValueError("Each neutron detector must be a mapping.")
        label = str(detector.get("label", f"A{index}"))
        angle = float(detector["angle_deg"])
        distance = float(detector["distance_cm"])
        entries.append(f"{label}:{angle:g}:{distance:g}")
    return ";".join(entries)


def tof_window_macro_value(config: Mapping[str, Any]) -> str:
    window = deep_get(config, "selection.tof_window_ns", [40.0, 50.0])
    if not isinstance(window, (list, tuple)) or len(window) != 2:
        raise ValueError("selection.tof_window_ns must contain exactly two values.")
    lo, hi = float(window[0]), float(window[1])
    return f"{lo:g} {hi:g}"


def flatten_for_template(config: Mapping[str, Any], row: Mapping[str, Any]) -> Dict[str, Any]:
    """Return all token values needed by macros/template_run.mac."""
    inner_r = float(deep_get(config, "geometry.inner_lar.radius_cm", 5.0))
    inner_h = float(deep_get(config, "geometry.inner_lar.height_cm", 10.0))

    outer_r = deep_get(config, "geometry.outer_lar.radius_cm", None)
    outer_d = deep_get(config, "geometry.outer_lar.diameter_cm", None)
    if outer_d is not None:
        outer_d = float(outer_d)
        outer_r = 0.5 * outer_d
    else:
        outer_r = float(outer_r if outer_r is not None else 17.5)
        outer_d = 2.0 * outer_r
    outer_h = float(deep_get(config, "geometry.outer_lar.height_cm", 70.0))

    pmt_pde = float(deep_get(config, "sensors.pmt.pde", 0.266))
    pmt_top_pde = float(deep_get(config, "sensors.pmt.top_pde", pmt_pde))
    pmt_bottom_pde = float(deep_get(config, "sensors.pmt.bottom_pde", pmt_top_pde))

    values: Dict[str, Any] = {
        "config_id": row.get("config_id", "config"),
        "label": row.get("label", row.get("config_id", "config")),
        "top_pmts": int(row.get("top_pmts", 0)),
        "bottom_pmts": int(row.get("bottom_pmts", 0)),
        "sipm_rows": int(row.get("sipm_rows", 0)),
        "n_events": int(deep_get(config, "run.n_events", 1000)),

        # Geometry
        "inner_lar_radius_cm": inner_r,
        "inner_lar_height_cm": inner_h,
        "inner_diameter_cm": 2.0 * inner_r,
        "inner_height_cm": inner_h,
        "outer_lar_radius_cm": outer_r,
        "outer_lar_diameter_cm": outer_d,
        "outer_lar_height_cm": outer_h,
        "outer_diameter_cm": outer_d,
        "outer_height_cm": outer_h,
        "neutron_detector_config": neutron_detector_macro_value(config),

        # Beam and event-selection cuts
        "primary_energy_MeV": float(deep_get(config, "beam.primary_energy_MeV", 2.5)),
        "tof_window_ns": tof_window_macro_value(config),

        # LAr optics/scintillation
        "lar_absorption_length_cm": float(deep_get(config, "optics.lar.absorption_length_cm", 150.0)),
        "lar_scintillation_yield_per_MeV": float(deep_get(config, "optics.lar.scintillation_yield_per_MeV", 51300.0)),
        "lar_scintillation_yield_scale": float(deep_get(config, "optics.lar.scintillation_yield_scale", 1.0)),
        "lar_birks_constant_mm_per_MeV": float(deep_get(config, "optics.lar.birks_constant_mm_per_MeV", 0.03)),
        "lar_fast_fraction": float(deep_get(config, "optics.lar.fast_fraction", 0.75)),
        "lar_fast_time_ns": float(deep_get(config, "optics.lar.fast_time_ns", 7.0)),
        "lar_slow_time_ns": float(deep_get(config, "optics.lar.slow_time_ns", 1500.0)),

        # TPB/reflector
        "tpb_efficiency": float(deep_get(config, "optics.tpb.efficiency", deep_get(config, "optics.tpb.conversion_efficiency", 1.0))),
        "tpb_conversion_efficiency": float(deep_get(config, "optics.tpb.efficiency", deep_get(config, "optics.tpb.conversion_efficiency", 1.0))),
        "reflector_reflectivity": float(deep_get(config, "optics.reflector.reflectivity", 0.95)),

        # Sensor geometry
        "pmt_side_length_cm": float(deep_get(config, "sensors.pmt.side_length_cm", 2.54)),
        "sipm_tile_size_cm": float(deep_get(config, "sensors.sipm.tile_size_cm", 0.6)),

        # Sensor response modes
        "pmt_response_mode": deep_get(config, "sensors.pmt.response_mode", "builtin_curve"),
        "sipm_response_mode": deep_get(config, "sensors.sipm.response_mode", "builtin_curve"),
        "pmt_qe_curve_file": deep_get(config, "sensors.pmt.qe_curve_file", "none"),
        "sipm_qe_curve_file": deep_get(config, "sensors.sipm.qe_curve_file", "none"),

        # Sensor optical settings
        "pmt_pde": pmt_pde,
        "pmt_top_pde": pmt_top_pde,
        "pmt_bottom_pde": pmt_bottom_pde,
        "sipm_pde": float(deep_get(config, "sensors.sipm.pde", 0.35)),
        "pmt_reflectivity": float(deep_get(config, "sensors.pmt.reflectivity", 0.10)),
        "sipm_reflectivity": float(deep_get(config, "sensors.sipm.reflectivity", 0.10)),
    }
    values.update(coverage_metrics(config, values))
    return values

def coverage_metrics(config: Mapping[str, Any], values: Mapping[str, Any]) -> Dict[str, float]:
    """Compute approximate optical coverage from the same geometric convention used
    in the thesis analysis: top+bottom circular faces plus cylindrical side area.
    """
    radius = float(deep_get(config, "geometry.inner_lar.radius_cm", 5.0))
    height = float(deep_get(config, "geometry.inner_lar.height_cm", 10.0))
    pmt_side = float(deep_get(config, "sensors.pmt.side_length_cm", 2.54))
    pmt_area = float(deep_get(config, "sensors.pmt.active_area_cm2", pmt_side * pmt_side))
    tile = float(deep_get(config, "sensors.sipm.tile_size_cm", 0.6))

    total_area = 2.0 * math.pi * radius * radius + 2.0 * math.pi * radius * height
    pmt_coverage_area = (int(values["top_pmts"]) + int(values["bottom_pmts"])) * pmt_area
    tiles_per_ring = max(1, int(round(2.0 * math.pi * radius / tile)))
    sipm_coverage_area = int(values["sipm_rows"]) * tiles_per_ring * tile * tile
    active_area = pmt_coverage_area + sipm_coverage_area
    return {
        "tiles_per_sipm_ring": float(tiles_per_ring),
        "active_area_cm2": active_area,
        "total_inner_surface_cm2": total_area,
        "coverage_fraction": active_area / total_area if total_area > 0 else 0.0,
        "coverage_percent": 100.0 * active_area / total_area if total_area > 0 else 0.0,
    }


_TOKEN_RE = re.compile(r"{{\s*([A-Za-z0-9_]+)\s*}}")


def render_template(template: str, values: Mapping[str, Any]) -> str:
    def repl(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in values:
            raise KeyError(f"Template token {{{{{key}}}}} has no value.")
        return str(values[key])
    return _TOKEN_RE.sub(repl, template)


def write_json(path: str | Path, data: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
