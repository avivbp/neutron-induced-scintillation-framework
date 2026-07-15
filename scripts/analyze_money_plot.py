#!/usr/bin/env python3
"""
Create mean-photoelectron versus optical-coverage plots.

For every optical configuration, the script computes:
  * the arithmetic mean number of photoelectrons;
  * the event-to-event sample standard deviation;
  * the standard error on the mean (SEM);
  * an optional percentile-bootstrap confidence interval for the mean.

Plots produced:
  1. One combined plot using all detector labels together.
  2. One plot for each configured or observed neutron detector.

All plots include the same hard-coded, detector-independent PE threshold bands.

The uncertainty shown on the plots is selected with:
  --uncertainty sem
or
  --uncertainty bootstrap

The bootstrap option replaces the plotted SEM error bars; both SEM and
bootstrap quantities are still written to the output summary CSV files.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from neutron_scintillation_correction import apply_correction


PE_COLUMN_CANDIDATES = ("numPE", "nPE", "pe", "PE")

# Hard-coded literature-motivated PE threshold bands.
# These are deliberately detector-independent.
THRESHOLD_BANDS = (
    (2.0, 4.0, "red", "Noise-dominated region"),
    (20.0, 30.0, "purple", "Trigger-threshold region"),
    (40.0, 80.0, "green", "PSD-useful region"),
)



def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream) or {}


def nested(config: dict, keys: tuple[str, ...], default):
    value = config
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def configured_neutron_detectors(config: dict) -> list[dict]:
    detectors = nested(config, ("geometry", "neutron_detectors"), None)
    if detectors is None:
        detectors = [
            {"label": "A0", "angle_deg": 25.0, "distance_cm": 100.0},
            {"label": "A1", "angle_deg": 40.0, "distance_cm": 100.0},
            {"label": "A2", "angle_deg": 50.0, "distance_cm": 100.0},
            {"label": "A3", "angle_deg": 60.0, "distance_cm": 100.0},
            {"label": "A4", "angle_deg": 90.0, "distance_cm": 100.0},
        ]
    if not isinstance(detectors, list):
        return []
    return [
        detector
        for detector in detectors
        if isinstance(detector, dict) and "label" in detector
    ]


def detector_labels_for_analysis(config: dict, observed: pd.Series) -> list[str]:
    labels = [
        str(detector["label"]).strip()
        for detector in configured_neutron_detectors(config)
    ]
    for detector in sorted(observed.dropna().astype(str).str.strip().unique()):
        if detector and detector not in labels:
            labels.append(detector)
    return labels


def detector_plot_title(config: dict, detector_label: str) -> str:
    for detector in configured_neutron_detectors(config):
        if str(detector.get("label", "")).strip() == detector_label:
            angle = detector.get("angle_deg")
            distance = detector.get("distance_cm")
            if angle is not None and distance is not None:
                return (
                    f"Neutron detector {detector_label}: "
                    f"{float(angle):g} deg, {float(distance):g} cm"
                )
    return f"Neutron detector {detector_label}"


def parse_optical_configuration(path: Path) -> tuple[int, int, int] | None:
    match = re.search(
        r"numPE_(\d+)_topPMTs_(\d+)_botPMTs_(\d+)_SiPMRows",
        path.name,
    )
    if match is None:
        return None
    return tuple(int(value) for value in match.groups())


def coverage_percent(config: dict, top: int, bottom: int, sipm_rows: int) -> float:
    radius_cm = float(
        nested(config, ("geometry", "inner_lar", "radius_cm"), 5.0)
    )
    height_cm = float(
        nested(config, ("geometry", "inner_lar", "height_cm"), 10.0)
    )
    pmt_side_cm = float(
        nested(config, ("sensors", "pmt", "side_length_cm"), 2.54)
    )
    sipm_tile_cm = float(
        nested(config, ("sensors", "sipm", "tile_size_cm"), 0.6)
    )

    total_area = (
        2.0 * math.pi * radius_cm**2
        + 2.0 * math.pi * radius_cm * height_cm
    )
    pmt_area = pmt_side_cm**2

    tiles_per_row = int((2.0 * math.pi * radius_cm) // sipm_tile_cm)
    sipm_row_area = tiles_per_row * sipm_tile_cm**2

    instrumented_area = (
        pmt_area * (top + bottom)
        + sipm_row_area * sipm_rows
    )
    return 100.0 * instrumented_area / total_area


def find_pe_column(frame: pd.DataFrame) -> str | None:
    for column in PE_COLUMN_CANDIDATES:
        if column in frame.columns:
            return column
    return None


def bootstrap_mean_interval(
    values: np.ndarray,
    confidence: float,
    samples: int,
    rng: np.random.Generator,
    chunk_size: int = 1000,
) -> tuple[float, float, float]:
    """
    Return (bootstrap standard error, lower percentile, upper percentile).

    Resampling is performed in chunks to avoid allocating a very large
    samples-by-events array for high-statistics runs.
    """
    n_events = values.size
    if n_events == 0:
        return np.nan, np.nan, np.nan
    if n_events == 1:
        value = float(values[0])
        return 0.0, value, value
    if samples < 2:
        raise ValueError("--bootstrap-samples must be at least 2")

    bootstrap_means = np.empty(samples, dtype=float)

    start = 0
    while start < samples:
        stop = min(start + chunk_size, samples)
        n_chunk = stop - start
        indices = rng.integers(
            low=0,
            high=n_events,
            size=(n_chunk, n_events),
        )
        bootstrap_means[start:stop] = values[indices].mean(axis=1)
        start = stop

    alpha = 1.0 - confidence
    lower_percentile = 100.0 * alpha / 2.0
    upper_percentile = 100.0 * (1.0 - alpha / 2.0)

    lower, upper = np.percentile(
        bootstrap_means,
        [lower_percentile, upper_percentile],
    )
    bootstrap_se = float(np.std(bootstrap_means, ddof=1))

    return bootstrap_se, float(lower), float(upper)


def summarize_values(
    values: np.ndarray,
    confidence: float,
    bootstrap_samples: int,
    rng: np.random.Generator,
) -> dict[str, float | int]:
    n_events = int(values.size)
    mean_pe = float(np.mean(values))
    std_pe = float(np.std(values, ddof=1)) if n_events > 1 else 0.0
    sem_pe = std_pe / math.sqrt(n_events) if n_events > 1 else 0.0

    bootstrap_se, bootstrap_low, bootstrap_high = bootstrap_mean_interval(
        values=values,
        confidence=confidence,
        samples=bootstrap_samples,
        rng=rng,
    )

    return {
        "n_selected_events": n_events,
        "mean_pe": mean_pe,
        "std_pe": std_pe,
        "sem_pe": sem_pe,
        "bootstrap_se_pe": bootstrap_se,
        "bootstrap_ci_low_pe": bootstrap_low,
        "bootstrap_ci_high_pe": bootstrap_high,
        "bootstrap_err_low_pe": mean_pe - bootstrap_low,
        "bootstrap_err_high_pe": bootstrap_high - mean_pe,
    }


def find_input_files(input_dir: Path) -> list[Path]:
    """
    Read only files directly in input_dir.

    This deliberately avoids recursively finding duplicate numPE files in
    build directories or old output directories.
    """
    return sorted(
        input_dir.glob(
            "numPE_*_topPMTs_*_botPMTs_*_SiPMRows.csv"
        )
    )


def load_summaries(
    input_dir: Path,
    config: dict,
    confidence: float,
    bootstrap_samples: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    files = find_input_files(input_dir)
    if not files:
        raise FileNotFoundError(
            f"No numPE files were found directly in {input_dir}"
        )
    print(f"[analysis] found {len(files)} numPE files in {input_dir}")

    combined_rows: list[dict] = []
    detector_rows: list[dict] = []

    # A reproducible but different random stream for each file/subset.
    master_rng = np.random.default_rng(seed)
    correction_config = nested(
        config,
        ("analysis", "neutron_scintillation_correction"),
        {},
    )
    correction_enabled = bool(
        isinstance(correction_config, dict)
        and correction_config.get("enabled", False)
    )

    for path in files:
        parsed = parse_optical_configuration(path)
        if parsed is None:
            print(f"[analysis] warning: cannot parse configuration from {path.name}")
            continue

        top, bottom, sipm_rows = parsed
        coverage = coverage_percent(config, top, bottom, sipm_rows)

        try:
            frame = pd.read_csv(path)
        except Exception as exc:
            print(f"[analysis] warning: could not read {path}: {exc}")
            continue

        # Fix headers such as "detector    ".
        frame.columns = frame.columns.astype(str).str.strip()

        pe_column = find_pe_column(frame)
        if pe_column is None:
            print(
                f"[analysis] warning: no PE column in {path.name}; "
                f"columns={frame.columns.tolist()}"
            )
            continue

        if correction_enabled:
            try:
                frame, correction_info = apply_correction(
                    frame,
                    correction_config,
                    pe_column,
                )
            except ValueError as exc:
                raise SystemExit(
                    f"Neutron scintillation correction failed for {path.name}: {exc}"
                ) from exc
            pe_column = "numPE_corrected"
            n_valid = int(frame[pe_column].notna().sum())
            leff_description = (
                "constant"
                if correction_info.extrapolation == "constant"
                else (
                    f"{correction_info.leff_min_keV:g}-"
                    f"{correction_info.leff_max_keV:g} keV "
                    f"({correction_info.extrapolation})"
                )
            )
            print(
                "[analysis] neutron correction: "
                f"{path.name}: {n_valid}/{len(frame)} valid events; "
                f"gamma a={correction_info.slope_photons_per_MeV:g} photons/MeV, "
                f"b={correction_info.intercept_photons:g} photons; "
                f"Leff={leff_description}"
            )

        pe_values = pd.to_numeric(
            frame[pe_column],
            errors="coerce",
        )
        valid_all = pe_values[np.isfinite(pe_values)].to_numpy(dtype=float)

        if valid_all.size == 0:
            print(f"[analysis] warning: no valid PE values in {path.name}")
        else:
            rng = np.random.default_rng(
                master_rng.integers(0, np.iinfo(np.uint32).max)
            )
            combined_rows.append(
                {
                    "source_file": path.name,
                    "top_pmts": top,
                    "bottom_pmts": bottom,
                    "sipm_rows": sipm_rows,
                    "coverage_percent": coverage,
                    **summarize_values(
                        valid_all,
                        confidence,
                        bootstrap_samples,
                        rng,
                    ),
                    **correction_diagnostics(frame),
                }
            )

        if "detector" not in frame.columns:
            print(
                f"[analysis] warning: no detector column in {path.name}; "
                "detector-specific plots cannot use this file"
            )
            continue

        detector_labels = frame["detector"].astype(str).str.strip()

        for detector in detector_labels_for_analysis(config, detector_labels):
            mask = detector_labels.eq(detector)
            detector_values = pd.to_numeric(
                frame.loc[mask, pe_column],
                errors="coerce",
            )
            detector_values = detector_values[
                np.isfinite(detector_values)
            ].to_numpy(dtype=float)

            if detector_values.size == 0:
                continue

            rng = np.random.default_rng(
                master_rng.integers(0, np.iinfo(np.uint32).max)
            )
            detector_rows.append(
                {
                    "detector": detector,
                    "source_file": path.name,
                    "top_pmts": top,
                    "bottom_pmts": bottom,
                    "sipm_rows": sipm_rows,
                    "coverage_percent": coverage,
                    **summarize_values(
                        detector_values,
                        confidence,
                        bootstrap_samples,
                        rng,
                    ),
                    **correction_diagnostics(frame.loc[mask]),
                }
            )

    combined = pd.DataFrame(combined_rows)
    by_detector = pd.DataFrame(detector_rows)

    if combined.empty:
        raise SystemExit(
            f"Found {len(files)} numPE files in {input_dir}, but none "
            "contained valid PE event rows. Check whether the simulation "
            "produced only headers or wrote event CSVs to another directory."
        )

    if not combined.empty:
        combined = combined.sort_values(
            ["coverage_percent", "top_pmts", "bottom_pmts", "sipm_rows"]
        ).reset_index(drop=True)

    if not by_detector.empty:
        by_detector = by_detector.sort_values(
            [
                "detector",
                "coverage_percent",
                "top_pmts",
                "bottom_pmts",
                "sipm_rows",
            ]
        ).reset_index(drop=True)

    return combined, by_detector


def correction_diagnostics(frame: pd.DataFrame) -> dict[str, float | int | str]:
    """Add auditable raw/corrected quantities without changing event cuts."""
    if "numPE_corrected" not in frame.columns:
        return {"pe_quantity": "raw numPE"}

    valid_frame = frame.loc[frame["numPE_corrected"].notna()]

    def finite_mean(column: str) -> float:
        values = pd.to_numeric(
            valid_frame[column], errors="coerce"
        ).to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        return float(np.mean(values)) if values.size else np.nan

    return {
        "pe_quantity": "thesis-corrected numPE",
        "n_correction_valid": len(valid_frame),
        "mean_pe_raw": finite_mean("numPE_raw"),
        "mean_pe_corrected": finite_mean("numPE_corrected"),
        "mean_edep_MeV": finite_mean("eDep"),
        "mean_num_photons_raw": finite_mean("numPhotons"),
        "mean_num_photons_expected_NR": finite_mean("numPhotons_expected_NR"),
        "mean_correction_factor": finite_mean(
            "neutron_scintillation_correction_factor"
        ),
    }


def detector_angles(config: dict) -> dict[str, float]:
    return {
        str(detector["label"]).strip(): float(detector["angle_deg"])
        for detector in configured_neutron_detectors(config)
        if detector.get("angle_deg") is not None
    }


def plot_raw_and_corrected_vs_angle(
    frame: pd.DataFrame,
    config: dict,
    output_dir: Path,
) -> None:
    """Make one thesis-correction audit plot per optical configuration."""
    if frame.empty or "mean_pe_corrected" not in frame.columns:
        return

    angle_by_detector = detector_angles(config)
    output_dir.mkdir(parents=True, exist_ok=True)

    for source_file, subset in frame.groupby("source_file", sort=False):
        subset = subset.copy()
        subset["angle_deg"] = subset["detector"].map(angle_by_detector)
        subset = subset.dropna(
            subset=["angle_deg", "mean_pe_raw", "mean_pe_corrected"]
        ).sort_values("angle_deg")
        if subset.empty:
            continue

        figure, axis = plt.subplots(figsize=(8, 5.5))
        axis.plot(
            subset["angle_deg"],
            subset["mean_pe_raw"],
            "o-",
            label="Raw detected PE",
        )
        axis.plot(
            subset["angle_deg"],
            subset["mean_pe_corrected"],
            "s-",
            label="Thesis-corrected PE",
        )
        axis.set_xlabel("Neutron detector angle [deg]")
        axis.set_ylabel("Mean number of photoelectrons")
        axis.set_title(Path(source_file).stem)
        axis.grid(True, alpha=0.3)
        axis.legend()
        figure.tight_layout()
        output_path = output_dir / f"{Path(source_file).stem}_pe_vs_angle.png"
        figure.savefig(output_path, dpi=300)
        plt.close(figure)
        print(f"[analysis] wrote {output_path}")


def y_errors(frame: pd.DataFrame, uncertainty: str) -> np.ndarray:
    if uncertainty == "sem":
        return frame["sem_pe"].to_numpy(dtype=float)

    lower = frame["bootstrap_err_low_pe"].to_numpy(dtype=float)
    upper = frame["bootstrap_err_high_pe"].to_numpy(dtype=float)
    return np.vstack([lower, upper])


def add_threshold_bands(axis: plt.Axes) -> float:
    """Add detector-independent PE threshold bands to an axes.

    Returns the highest band edge so the caller can ensure that all bands
    remain visible after Matplotlib autoscales the plotted data.
    """
    highest_edge = 0.0

    for lower, upper, color, label in THRESHOLD_BANDS:
        axis.axhspan(
            lower,
            upper,
            color=color,
            alpha=0.13,
            label=label,
            zorder=0,
        )
        highest_edge = max(highest_edge, upper)

    return highest_edge


def deduplicate_legend(axis: plt.Axes, title: str | None = None) -> None:
    """Remove duplicate legend entries while preserving their order."""
    handles, labels = axis.get_legend_handles_labels()
    unique: dict[str, object] = {}

    for handle, label in zip(handles, labels):
        if label and label not in unique:
            unique[label] = handle

    axis.legend(
        unique.values(),
        unique.keys(),
        title=title,
    )


def plot_summary(
    frame: pd.DataFrame,
    output_path: Path,
    uncertainty: str,
    confidence: float,
    title: str | None = None,
    corrected: bool = False,
) -> None:
    if frame.empty:
        return

    frame = frame.sort_values("coverage_percent")

    coverage = frame["coverage_percent"].to_numpy(dtype=float)
    mean_pe = frame["mean_pe"].to_numpy(dtype=float)
    errors = y_errors(frame, uncertainty)
    sipm_rows = frame["sipm_rows"].to_numpy(dtype=int)

    pmt_only = sipm_rows == 0
    with_sipm = sipm_rows > 0

    figure, axis = plt.subplots(figsize=(9, 6))
    highest_band_edge = add_threshold_bands(axis)

    if np.any(pmt_only):
        pmt_errors = (
            errors[pmt_only]
            if errors.ndim == 1
            else errors[:, pmt_only]
        )
        axis.errorbar(
            coverage[pmt_only],
            mean_pe[pmt_only],
            yerr=pmt_errors,
            fmt="o",
            capsize=3,
            label="PMTs only",
        )

    if np.any(with_sipm):
        sipm_errors = (
            errors[with_sipm]
            if errors.ndim == 1
            else errors[:, with_sipm]
        )
        axis.errorbar(
            coverage[with_sipm],
            mean_pe[with_sipm],
            yerr=sipm_errors,
            fmt="s",
            capsize=3,
            label="PMTs + SiPM rows",
        )

    if title:
        axis.set_title(title)

    axis.set_xlabel("Optical coverage [%]")
    axis.set_ylabel(
        "Mean corrected number of photoelectrons"
        if corrected
        else "Mean number of photoelectrons"
    )
    axis.grid(True, alpha=0.3)

    # Keep all threshold bands visible without clipping higher-valued data.
    current_lower, current_upper = axis.get_ylim()
    axis.set_ylim(0.0, max(current_upper, highest_band_edge))

    uncertainty_label = (
        "SEM"
        if uncertainty == "sem"
        else f"{100.0 * confidence:g}% bootstrap CI"
    )
    deduplicate_legend(
        axis,
        title=f"Error bars: {uncertainty_label}",
    )

    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)

    print(f"[analysis] wrote {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot mean photoelectrons versus optical coverage, combined "
            "and separately for each neutron detector."
        )
    )
    parser.add_argument(
        "--config",
        default="config/example_config.yaml",
        help="YAML configuration used for geometry and sensor dimensions.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help=(
            "Directory containing numPE CSV files. By default, use "
            "run.output_dir from the selected config."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for summary CSVs and figures. By default, use "
            "run.output_dir from the selected config."
        ),
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="Accepted for compatibility with older run_scan.py versions.",
    )
    parser.add_argument(
        "--uncertainty",
        choices=("sem", "bootstrap"),
        default="sem",
        help=(
            "Error bars to plot: SEM, or a percentile-bootstrap confidence "
            "interval for the mean. Default: sem"
        ),
    )
    parser.add_argument(
        "--bootstrap-confidence",
        type=float,
        default=0.68,
        help="Bootstrap confidence level as a fraction. Default: 0.68",
    )
    parser.add_argument(
        "--bootstrap-samples",
        type=int,
        default=10000,
        help="Number of bootstrap resamples. Default: 10000",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=12345,
        help="Random seed for reproducible bootstrap intervals.",
    )

    args = parser.parse_args()

    if not 0.0 < args.bootstrap_confidence < 1.0:
        parser.error("--bootstrap-confidence must be between 0 and 1")
    if args.bootstrap_samples < 2:
        parser.error("--bootstrap-samples must be at least 2")

    config_path = Path(args.config)
    config = load_yaml(config_path)
    configured_output_dir = nested(
        config, ("run", "output_dir"), "output/tau_lab_example"
    )
    input_dir = Path(args.input_dir or configured_output_dir)
    output_dir = Path(args.output_dir or configured_output_dir)
    correction_enabled = bool(
        nested(
            config,
            ("analysis", "neutron_scintillation_correction", "enabled"),
            False,
        )
    )
    print(
        "[analysis] neutron scintillation correction: "
        + ("enabled" if correction_enabled else "disabled")
    )

    print(f"[analysis] uncertainty mode: {args.uncertainty}")
    if args.uncertainty == "bootstrap":
        print(
            "[analysis] bootstrap settings: "
            f"confidence={args.bootstrap_confidence}, "
            f"samples={args.bootstrap_samples}, "
            f"seed={args.bootstrap_seed}"
        )

    combined, by_detector = load_summaries(
        input_dir=input_dir,
        config=config,
        confidence=args.bootstrap_confidence,
        bootstrap_samples=args.bootstrap_samples,
        seed=args.bootstrap_seed,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    detector_figures_dir = figures_dir / "by_neutron_detector"

    combined.to_csv(
        output_dir / "run_summary.csv",
        index=False,
    )
    by_detector.to_csv(
        output_dir / "run_summary_by_detector.csv",
        index=False,
    )

    if correction_enabled:
        plot_raw_and_corrected_vs_angle(
            by_detector,
            config,
            figures_dir / "correction_vs_detector_angle",
        )

    plot_summary(
        combined,
        figures_dir / "coverage_vs_mean_pe.png",
        uncertainty=args.uncertainty,
        confidence=args.bootstrap_confidence,
        corrected=correction_enabled,
    )

    if by_detector.empty or "detector" not in by_detector.columns:
        print(
            "[analysis] warning: no detector-specific data were found; "
            "no neutron-detector plots were created"
        )
        return

    for detector in by_detector["detector"].dropna().astype(str).unique():
        subset = by_detector[by_detector["detector"] == detector].copy()
        if subset.empty:
            print(f"[analysis] warning: no data found for detector {detector}")
            continue

        plot_summary(
            subset,
            detector_figures_dir
            / f"coverage_vs_mean_pe_detector_{detector}.png",
            uncertainty=args.uncertainty,
            confidence=args.bootstrap_confidence,
            title=detector_plot_title(config, detector),
            corrected=correction_enabled,
        )


if __name__ == "__main__":
    main()
