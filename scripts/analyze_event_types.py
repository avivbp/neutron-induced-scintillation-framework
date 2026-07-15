#!/usr/bin/env python3
"""Summarize neutron event types among rows that passed simulation cuts.

Each neutron ``numPE_*.csv`` contains only events written after the detected,
positive-PE, and TOF requirements.  This script classifies those rows using
``numElasticSensitive`` and ``scatteredNotSensitive`` into mutually exclusive
event types and reports counts and relative frequencies.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


FILE_PATTERN = "numPE_*_topPMTs_*_botPMTs_*_SiPMRows.csv"
CATEGORIES = (
    "true_single_elastic",
    "single_elastic_external_scatter",
    "multiple_elastic",
    "other_or_inelastic",
)
LABELS = {
    "true_single_elastic": "True single elastic",
    "single_elastic_external_scatter": "Single elastic + external scatter",
    "multiple_elastic": "Multiple elastic",
    "other_or_inelastic": "Other / inelastic / capture",
}
COLORS = {
    "true_single_elastic": "#2ca02c",
    "single_elastic_external_scatter": "#ff7f0e",
    "multiple_elastic": "#d62728",
    "other_or_inelastic": "#7f7f7f",
}


def find_run_dirs(paths: list[Path]) -> list[Path]:
    """Resolve explicit run dirs or discover direct children containing CSVs."""
    run_dirs: list[Path] = []
    for path in paths:
        if list(path.glob(FILE_PATTERN)):
            run_dirs.append(path)
            continue
        if not path.is_dir():
            print(f"[event types] warning: not a directory: {path}")
            continue
        run_dirs.extend(
            child
            for child in sorted(path.iterdir())
            if child.is_dir() and list(child.glob(FILE_PATTERN))
        )
    # Preserve discovery order while removing duplicate resolved paths.
    return list(dict.fromkeys(path.resolve() for path in run_dirs))


def as_boolean(values: pd.Series) -> pd.Series:
    """Interpret the numeric/bool spellings expected in simulation CSVs."""
    numeric = pd.to_numeric(values, errors="coerce")
    result = numeric.fillna(0).ne(0)
    text = values.astype(str).str.strip().str.lower()
    result.loc[text.isin(("true", "yes", "y"))] = True
    return result


def load_and_classify(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame.columns = frame.columns.astype(str).str.strip()
    required = {"numElasticSensitive", "scatteredNotSensitive"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")

    elastic = pd.to_numeric(frame["numElasticSensitive"], errors="coerce")
    external = as_boolean(frame["scatteredNotSensitive"])
    clean = pd.Series(True, index=frame.index)
    for column in ("numInelastic", "nCapture", "InelasticSensitive"):
        if column in frame.columns:
            clean &= pd.to_numeric(frame[column], errors="coerce").fillna(0).eq(0)

    category = np.full(len(frame), "other_or_inelastic", dtype=object)
    category[(elastic == 1) & ~external & clean] = "true_single_elastic"
    category[(elastic == 1) & external & clean] = "single_elastic_external_scatter"
    # Multiple fiducial elastics take precedence regardless of external scatter.
    category[(elastic >= 2) & clean] = "multiple_elastic"

    result = frame.copy()
    result["event_type"] = category
    result["source_file"] = path.name
    return result


def summary_rows(
    frame: pd.DataFrame,
    run_dir: Path,
    grouping: str,
    group_value: str,
    scattering_angle_deg: float | None = None,
) -> list[dict]:
    counts = frame["event_type"].value_counts()
    total = len(frame)
    rows = []
    for category in CATEGORIES:
        count = int(counts.get(category, 0))
        rows.append(
            {
                "run_directory": str(run_dir),
                "run_name": run_dir.name,
                "grouping": grouping,
                "group": group_value,
                "scattering_angle_deg": scattering_angle_deg,
                "event_type": category,
                "event_type_label": LABELS[category],
                "count": count,
                "total_passed_events": total,
                "relative_frequency": count / total if total else np.nan,
                "percent": 100.0 * count / total if total else np.nan,
            }
        )
    return rows


def configured_detector_angles(run_dir: Path) -> dict[str, float]:
    config_path = run_dir / "config_used.yaml"
    if not config_path.exists():
        return {}
    try:
        with config_path.open("r", encoding="utf-8") as stream:
            config = yaml.safe_load(stream) or {}
        detectors = config.get("geometry", {}).get("neutron_detectors", [])
        return {
            str(detector["label"]).strip(): float(detector["angle_deg"])
            for detector in detectors
            if isinstance(detector, dict)
            and "label" in detector
            and "angle_deg" in detector
        }
    except (OSError, TypeError, ValueError, yaml.YAMLError) as exc:
        print(f"[event types] warning: could not read {config_path}: {exc}")
        return {}


def summarize_run(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    classified_frames: list[pd.DataFrame] = []
    rows: list[dict] = []
    for path in sorted(run_dir.glob(FILE_PATTERN)):
        try:
            frame = load_and_classify(path)
        except (OSError, pd.errors.ParserError, ValueError) as exc:
            print(f"[event types] warning: {exc}")
            continue
        if frame.empty:
            print(f"[event types] warning: no passed events in {path}")
        classified_frames.append(frame)
        rows.extend(summary_rows(frame, run_dir, "file", path.name))

    if not classified_frames:
        return pd.DataFrame(), pd.DataFrame()

    combined = pd.concat(classified_frames, ignore_index=True)
    rows.extend(summary_rows(combined, run_dir, "run", "all_files"))
    if "detector" in combined.columns:
        angle_by_detector = configured_detector_angles(run_dir)
        detector = combined["detector"].astype(str).str.strip()
        for label in sorted(detector[detector.ne("")].unique()):
            rows.extend(
                summary_rows(
                    combined.loc[detector == label],
                    run_dir,
                    "detector",
                    label,
                    angle_by_detector.get(label),
                )
            )
    return pd.DataFrame(rows), combined


def plot_run_frequencies(summary: pd.DataFrame, output_path: Path) -> None:
    run_rows = summary[summary["grouping"] == "run"].copy()
    if run_rows.empty:
        return
    pivot = run_rows.pivot(
        index="run_name", columns="event_type", values="percent"
    ).reindex(columns=CATEGORIES, fill_value=0.0)
    axis = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(max(8.0, 1.4 * len(pivot)), 6.0),
        color=[COLORS[value] for value in pivot.columns],
    )
    axis.set_xlabel("Neutron run directory")
    axis.set_ylabel("Fraction of passed events [%]")
    axis.set_ylim(0.0, 100.0)
    axis.grid(axis="y", alpha=0.3)
    axis.legend(
        [LABELS[value] for value in pivot.columns],
        title="Event type",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
    )
    axis.figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(output_path, dpi=300)
    plt.close(axis.figure)


def natural_key(value: str) -> list[object]:
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    ]


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "run"


def plot_detector_frequencies(
    summary: pd.DataFrame,
    run_name: str,
    output_path: Path,
) -> None:
    detector_rows = summary[
        (summary["grouping"] == "detector")
        & (summary["run_name"] == run_name)
    ].copy()
    if detector_rows.empty:
        return

    order_frame = detector_rows[["group", "scattering_angle_deg"]].drop_duplicates()
    angle_by_detector = order_frame.set_index("group")["scattering_angle_deg"]
    detector_order = sorted(
        detector_rows["group"].unique(),
        key=lambda label: (
            pd.isna(angle_by_detector.get(label)),
            angle_by_detector.get(label)
            if not pd.isna(angle_by_detector.get(label))
            else float("inf"),
            natural_key(label),
        ),
    )
    pivot = detector_rows.pivot(
        index="group", columns="event_type", values="percent"
    ).reindex(index=detector_order, columns=CATEGORIES, fill_value=0.0)
    totals = (
        detector_rows.groupby("group")["total_passed_events"]
        .first()
        .reindex(detector_order)
    )

    axis = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(max(8.0, 1.25 * len(pivot)), 6.3),
        color=[COLORS[value] for value in pivot.columns],
    )
    for index, total in enumerate(totals):
        axis.text(
            index,
            101.5,
            f"N={int(total)}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    tick_labels = []
    for label in detector_order:
        angle = angle_by_detector.get(label)
        tick_labels.append(
            f"{label}\n{float(angle):g}°" if not pd.isna(angle) else label
        )
    axis.set_xticklabels(tick_labels, rotation=0)
    axis.set_xlabel("Neutron detector and scattering angle")
    axis.set_ylabel("Fraction of passed events [%]")
    axis.set_title(f"Event types by scattering detector: {run_name}")
    axis.set_ylim(0.0, 112.0)
    axis.grid(axis="y", alpha=0.3)
    axis.legend(
        [LABELS[value] for value in pivot.columns],
        title="Event type",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
    )
    axis.figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(output_path, dpi=300)
    plt.close(axis.figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_dirs",
        nargs="*",
        type=Path,
        default=[Path("output")],
        help=(
            "Neutron run directories, or parent directories whose direct "
            "children are neutron runs. Default: output"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/event_type_analysis"),
        help="Summary CSV and plot directory.",
    )
    args = parser.parse_args()

    run_dirs = find_run_dirs(args.run_dirs)
    if not run_dirs:
        raise SystemExit("No neutron run directories containing numPE CSVs found")

    summaries: list[pd.DataFrame] = []
    for run_dir in run_dirs:
        summary, combined = summarize_run(run_dir)
        if summary.empty:
            continue
        summaries.append(summary)
        run_total = summary[
            (summary["grouping"] == "run")
            & (summary["event_type"] == CATEGORIES[0])
        ]["total_passed_events"].iloc[0]
        print(
            f"[event types] {run_dir}: {run_total} passed rows across "
            f"{combined['source_file'].nunique()} CSV files"
        )

    if not summaries:
        raise SystemExit("No valid neutron event rows were found")

    result = pd.concat(summaries, ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "event_type_frequencies.csv"
    result.to_csv(summary_path, index=False)
    plot_path = args.output_dir / "event_type_frequencies_by_run.png"
    plot_run_frequencies(result, plot_path)
    detector_plot_paths: list[Path] = []
    for run_name in result["run_name"].unique():
        detector_path = (
            args.output_dir
            / f"event_type_frequencies_by_detector_{safe_filename(run_name)}.png"
        )
        plot_detector_frequencies(result, run_name, detector_path)
        if detector_path.exists():
            detector_plot_paths.append(detector_path)
    print(f"[event types] wrote {summary_path}")
    print(f"[event types] wrote {plot_path}")
    for detector_path in detector_plot_paths:
        print(f"[event types] wrote {detector_path}")


if __name__ == "__main__":
    main()
