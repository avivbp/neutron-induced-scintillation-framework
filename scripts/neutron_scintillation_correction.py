#!/usr/bin/env python3
"""Event-level neutron scintillation correction used by the analysis.

The implemented prescription is

    numPE_corrected = numPE * Nph_expected / numPhotons
    Nph_expected = (a * eDep_MeV + b) * Leff(eDep_keV)

where ``a`` and ``b`` either come directly from the analysis configuration or
are fitted from configured gamma-calibration points.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CorrectionInfo:
    slope_photons_per_MeV: float
    intercept_photons: float
    leff_min_keV: float
    leff_max_keV: float
    extrapolation: str


def _points(value, name: str) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(f"{name}.points must contain at least two [x, y] pairs")

    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}.points must be numeric [x, y] pairs") from exc

    if array.ndim != 2 or array.shape[1] != 2 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name}.points must be finite numeric [x, y] pairs")

    order = np.argsort(array[:, 0])
    x = array[order, 0]
    y = array[order, 1]
    if np.any(np.diff(x) <= 0):
        raise ValueError(f"{name}.points energies must be unique")
    return x, y


def _gamma_line(config: dict) -> tuple[float, float]:
    gamma = config.get("gamma_calibration", {})
    if not isinstance(gamma, dict):
        raise ValueError("gamma_calibration must be a mapping")

    if "slope_photons_per_MeV" in gamma:
        slope = float(gamma["slope_photons_per_MeV"])
        intercept = float(gamma.get("intercept_photons", 0.0))
    else:
        energies, photons = _points(gamma.get("points"), "gamma_calibration")
        unit = str(gamma.get("energy_unit", "MeV")).lower()
        if unit == "kev":
            energies = energies / 1000.0
        elif unit != "mev":
            raise ValueError("gamma_calibration.energy_unit must be keV or MeV")
        slope, intercept = np.polyfit(energies, photons, deg=1)

    if not np.isfinite(slope) or not np.isfinite(intercept) or slope <= 0:
        raise ValueError("the gamma calibration line must have a positive finite slope")
    return float(slope), float(intercept)


def _interpolate_leff(
    energy_keV: np.ndarray,
    point_energy_keV: np.ndarray,
    point_leff: np.ndarray,
    extrapolation: str,
) -> np.ndarray:
    if extrapolation == "clip":
        return np.interp(energy_keV, point_energy_keV, point_leff)

    result = np.interp(energy_keV, point_energy_keV, point_leff)
    outside = (energy_keV < point_energy_keV[0]) | (energy_keV > point_energy_keV[-1])

    if extrapolation == "nan":
        result[outside] = np.nan
        return result

    if extrapolation == "error":
        if np.any(outside):
            observed = energy_keV[np.isfinite(energy_keV)]
            raise ValueError(
                "eDep falls outside the configured Leff range "
                f"[{point_energy_keV[0]:g}, {point_energy_keV[-1]:g}] keV; "
                f"observed range is [{observed.min():g}, {observed.max():g}] keV"
            )
        return result

    if extrapolation == "linear":
        below = energy_keV < point_energy_keV[0]
        above = energy_keV > point_energy_keV[-1]
        low_slope = (
            (point_leff[1] - point_leff[0])
            / (point_energy_keV[1] - point_energy_keV[0])
        )
        high_slope = (
            (point_leff[-1] - point_leff[-2])
            / (point_energy_keV[-1] - point_energy_keV[-2])
        )
        result[below] = point_leff[0] + low_slope * (
            energy_keV[below] - point_energy_keV[0]
        )
        result[above] = point_leff[-1] + high_slope * (
            energy_keV[above] - point_energy_keV[-1]
        )
        return result

    raise ValueError("leff.extrapolation must be error, nan, clip, or linear")


def apply_correction(
    frame: pd.DataFrame,
    config: dict,
    pe_column: str,
) -> tuple[pd.DataFrame, CorrectionInfo]:
    """Return a copy of *frame* with raw, expected, and corrected PE columns."""
    for column in ("numPhotons", "eDep"):
        if column not in frame.columns:
            raise ValueError(
                f"neutron scintillation correction requires CSV column {column!r}"
            )

    slope, intercept = _gamma_line(config)
    leff_config = config.get("leff", {})
    if not isinstance(leff_config, dict):
        raise ValueError("leff must be a mapping")
    result = frame.copy()
    raw_pe = pd.to_numeric(result[pe_column], errors="coerce").to_numpy(dtype=float)
    raw_photons = pd.to_numeric(result["numPhotons"], errors="coerce").to_numpy(dtype=float)
    edep_MeV = pd.to_numeric(result["eDep"], errors="coerce").to_numpy(dtype=float)
    edep_keV = 1000.0 * edep_MeV

    if "constant" in leff_config:
        constant_leff = float(leff_config["constant"])
        if not np.isfinite(constant_leff) or constant_leff < 0:
            raise ValueError("leff.constant must be a finite non-negative number")
        leff = np.full(len(result), constant_leff, dtype=float)
        leff_min_keV = 0.0
        leff_max_keV = float("inf")
        extrapolation = "constant"
    else:
        leff_energy, leff_values = _points(leff_config.get("points"), "leff")
        unit = str(leff_config.get("energy_unit", "keV")).lower()
        if unit == "mev":
            leff_energy = leff_energy * 1000.0
        elif unit != "kev":
            raise ValueError("leff.energy_unit must be keV or MeV")
        if np.any(leff_values < 0):
            raise ValueError("Leff values must be non-negative")
        extrapolation = str(leff_config.get("extrapolation", "error")).lower()
        leff = _interpolate_leff(
            edep_keV, leff_energy, leff_values, extrapolation
        )
        leff_min_keV = float(leff_energy[0])
        leff_max_keV = float(leff_energy[-1])
    expected_gamma_photons = slope * edep_MeV + intercept
    expected_nr_photons = expected_gamma_photons * leff

    valid = (
        np.isfinite(raw_pe)
        & np.isfinite(raw_photons)
        & (raw_photons > 0)
        & np.isfinite(edep_MeV)
        & (edep_MeV >= 0)
        & np.isfinite(expected_nr_photons)
        & (expected_nr_photons >= 0)
    )
    factor = np.full(len(result), np.nan, dtype=float)
    corrected_pe = np.full(len(result), np.nan, dtype=float)
    factor[valid] = expected_nr_photons[valid] / raw_photons[valid]
    corrected_pe[valid] = raw_pe[valid] * factor[valid]

    result["numPE_raw"] = raw_pe
    result["Leff"] = leff
    result["numPhotons_expected_gamma"] = expected_gamma_photons
    result["numPhotons_expected_NR"] = expected_nr_photons
    result["neutron_scintillation_correction_factor"] = factor
    result["numPE_corrected"] = corrected_pe

    return result, CorrectionInfo(
        slope_photons_per_MeV=slope,
        intercept_photons=intercept,
        leff_min_keV=leff_min_keV,
        leff_max_keV=leff_max_keV,
        extrapolation=extrapolation,
    )
