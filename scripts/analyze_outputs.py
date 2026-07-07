#!/usr/bin/env python3
"""Aggregate Geant4 CSV outputs into package-level summary tables and figures.

This script is deliberately tolerant of Aviv's older output naming scheme
(numPE_<top>_topPMTs_<bottom>_botPMTs_<rows>_SiPMRows.csv) and of newer
package-style event_summary.csv files.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List

import matplotlib.pyplot as plt
import pandas as pd

from config_tools import deep_get, read_yaml


OLD_NUMPE_RE = re.compile(
    r"numPE_(?P<top>\d+)_topPMTs_(?P<bottom>\d+)_botPMTs_(?P<rows>\d+)_SiPMRows\.csv"
)


def read_metadata(path: Path) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return pd.DataFrame(data)


def find_old_numpe_files(input_dir: Path) -> list[Path]:
    return [Path(p) for p in glob.glob(str(input_dir / "**" / "numPE_*_topPMTs_*_botPMTs_*_SiPMRows.csv"), recursive=True)]


def load_event_file(path: Path) -> pd.DataFrame:
    # Older files often have no header.  This list follows the content written in EventAction.cc.
    old_cols = [
        "event_id", "num_pe", "num_photons", "edep_MeV",
        "num_elastic_sensitive", "num_inelastic_sensitive", "num_inelastic_total",
        "n_capture", "tof_ns", "ls_detector", "recoil_energy_MeV",
    ]
    try:
        df = pd.read_csv(path)
        # If the first column name looks numeric, pandas treated the first data row as a header.
        if all(str(c).strip().replace(".", "", 1).replace("-", "", 1).isdigit() for c in df.columns[:2]):
            df = pd.read_csv(path, header=None)
    except pd.errors.ParserError:
        df = pd.read_csv(path, header=None)
    if list(df.columns) == list(range(len(df.columns))):
        df.columns = old_cols[: len(df.columns)]
    return df


def attach_config_from_filename(df: pd.DataFrame, path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    match = OLD_NUMPE_RE.search(path.name)
    if match:
        top = int(match.group("top"))
        bottom = int(match.group("bottom"))
        rows = int(match.group("rows"))
        sel = metadata[
            (metadata["top_pmts"].astype(int) == top)
            & (metadata["bottom_pmts"].astype(int) == bottom)
            & (metadata["sipm_rows"].astype(int) == rows)
        ]
        if len(sel) > 0:
            info = sel.iloc[0].to_dict()
            for key, val in info.items():
                if key not in df.columns:
                    df[key] = val
        else:
            df["top_pmts"] = top
            df["bottom_pmts"] = bottom
            df["sipm_rows"] = rows
            df["config_id"] = f"top{top}_bot{bottom}_rows{rows}"
    else:
        if "config_id" not in df.columns:
            df["config_id"] = path.stem
    return df


def apply_selection(df: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
    sel = pd.Series(True, index=df.index)
    cuts = config.get("selection", {})
    if cuts.get("require_num_pe_positive", True) and "num_pe" in df:
        sel &= pd.to_numeric(df["num_pe"], errors="coerce").fillna(0) > 0
    if "tof_ns" in df and "tof_window_ns" in cuts and cuts["tof_window_ns"]:
        lo, hi = cuts["tof_window_ns"]
        tof = pd.to_numeric(df["tof_ns"], errors="coerce")
        sel &= tof.between(float(lo), float(hi), inclusive="both")
    if cuts.get("reject_inelastic_sensitive", False) and "num_inelastic_sensitive" in df:
        sel &= pd.to_numeric(df["num_inelastic_sensitive"], errors="coerce").fillna(0) == 0
    return df[sel].copy()


def summarize(df: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cid, group in df.groupby("config_id", dropna=False):
        sg = selected[selected["config_id"] == cid]
        meta = group.iloc[0].to_dict()
        num_pe = pd.to_numeric(sg.get("num_pe", pd.Series(dtype=float)), errors="coerce")
        rows.append({
            "config_id": cid,
            "label": meta.get("label", cid),
            "top_pmts": meta.get("top_pmts", None),
            "bottom_pmts": meta.get("bottom_pmts", None),
            "sipm_rows": meta.get("sipm_rows", None),
            "coverage_fraction": meta.get("coverage_fraction", None),
            "coverage_percent": meta.get("coverage_percent", None),
            "n_raw_events_in_file": len(group),
            "n_selected_events": len(sg),
            "mean_pe": float(num_pe.mean()) if len(num_pe) else float("nan"),
            "std_pe": float(num_pe.std(ddof=1)) if len(num_pe) > 1 else float("nan"),
            "sem_pe": float(num_pe.std(ddof=1) / math.sqrt(len(num_pe))) if len(num_pe) > 1 else float("nan"),
            "median_pe": float(num_pe.median()) if len(num_pe) else float("nan"),
        })
    out = pd.DataFrame(rows)
    if "coverage_fraction" in out:
        out = out.sort_values(["coverage_fraction", "config_id"], na_position="last")
    return out


def plot_coverage(summary: pd.DataFrame, output_dir: Path, fmt: str = "png") -> None:
    figdir = output_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    clean = summary.dropna(subset=["coverage_percent", "mean_pe"])
    if clean.empty:
        print("No coverage/PE data available for plotting.")
        return
    plt.figure(figsize=(7, 4.5))
    yerr = clean["sem_pe"] if "sem_pe" in clean else None
    plt.errorbar(clean["coverage_percent"], clean["mean_pe"], yerr=yerr, marker="o", linestyle="none", capsize=3)
    plt.xlabel("Optical coverage [%]")
    plt.ylabel("Mean detected photoelectrons")
    plt.title("Detected PE yield vs optical coverage")
    plt.tight_layout()
    plt.savefig(figdir / f"coverage_vs_mean_pe.{fmt}", dpi=200)
    plt.close()


def plot_pe_histograms(selected: pd.DataFrame, output_dir: Path, fmt: str = "png") -> None:
    figdir = output_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    if "num_pe" not in selected:
        return
    plt.figure(figsize=(7, 4.5))
    for cid, group in selected.groupby("config_id"):
        vals = pd.to_numeric(group["num_pe"], errors="coerce").dropna()
        if len(vals):
            plt.hist(vals, bins=40, histtype="step", label=str(cid))
    plt.xlabel("Detected PE")
    plt.ylabel("Events")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / f"pe_distributions.{fmt}", dpi=200)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/example_config.yaml")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--input-dir", default=".")
    parser.add_argument("--output-dir", default="output/scan")
    args = parser.parse_args()

    config = read_yaml(args.config)
    metadata = read_metadata(Path(args.metadata))
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    for path in find_old_numpe_files(input_dir):
        df = load_event_file(path)
        df = attach_config_from_filename(df, path, metadata)
        frames.append(df)

    if not frames:
        raise SystemExit("No Geant4 event CSV files found. Check analysis.expected_event_file_patterns or output paths.")

    raw = pd.concat(frames, ignore_index=True)
    selected = apply_selection(raw, config)
    summary = summarize(raw, selected)

    raw.to_csv(output_dir / "event_summary_raw.csv", index=False)
    selected.to_csv(output_dir / "event_summary_selected.csv", index=False)
    summary.to_csv(output_dir / "run_summary.csv", index=False)

    fmt = config.get("analysis", {}).get("figure_format", "png")
    plot_coverage(summary, output_dir, fmt)
    plot_pe_histograms(selected, output_dir, fmt)
    print(f"wrote {output_dir / 'run_summary.csv'}")


if __name__ == "__main__":
    main()
