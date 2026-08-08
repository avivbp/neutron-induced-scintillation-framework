#!/usr/bin/env python3
"""Summarize uncut DUNE-like optical response and apply NR-only correction."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from neutron_scintillation_correction import apply_correction


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream) or {}


def finite_mean(frame: pd.DataFrame, column: str) -> float:
    if column not in frame:
        return float("nan")
    return float(pd.to_numeric(frame[column], errors="coerce").mean())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    config = load_yaml(Path(args.config))
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    correction = config.get("analysis", {}).get(
        "neutron_scintillation_correction", {}
    )
    correction_enabled = bool(correction.get("enabled", False))
    paths = sorted(
        path
        for path in input_dir.glob("dune_event_response_*.csv")
        if not path.stem.endswith("_corrected")
    )
    if not paths:
        raise SystemExit(f"No DUNE event-response CSVs found in {input_dir}")

    summaries: list[dict] = []
    for path in paths:
        frame = pd.read_csv(path)
        if "total_num_pe" not in frame:
            raise ValueError(f"{path} is missing total_num_pe")
        pe_components = [
            "electronic_recoil_num_pe",
            "nuclear_recoil_num_pe",
            "other_num_pe",
        ]
        missing = [column for column in pe_components if column not in frame]
        if missing:
            raise ValueError(f"{path} is missing PE components: {missing}")
        total_pe = pd.to_numeric(frame["total_num_pe"], errors="coerce")
        component_pe = sum(
            pd.to_numeric(frame[column], errors="coerce")
            for column in pe_components
        )
        if not total_pe.equals(component_pe):
            raise ValueError(f"{path} has PE components that do not sum to total_num_pe")

        analyzed = frame
        corrected_path = None
        if correction_enabled:
            analyzed, _ = apply_correction(
                frame, correction, pe_column="total_num_pe"
            )
            corrected_path = output_dir / f"{path.stem}_corrected.csv"
            analyzed.to_csv(corrected_path, index=False)

        total_pe = pd.to_numeric(
            analyzed["total_num_pe"], errors="coerce"
        )
        summary = {
            "config_id": path.stem.removeprefix("dune_event_response_"),
            "n_events": len(analyzed),
            "n_events_with_pe": int((total_pe > 0).sum()),
            "mean_total_num_pe_raw": finite_mean(analyzed, "total_num_pe"),
            "mean_electronic_recoil_num_pe_raw": finite_mean(
                analyzed, "electronic_recoil_num_pe"
            ),
            "mean_nuclear_recoil_num_pe_raw": finite_mean(
                analyzed, "nuclear_recoil_num_pe"
            ),
            "mean_other_num_pe_raw": finite_mean(analyzed, "other_num_pe"),
            "correction_enabled": correction_enabled,
            "mean_total_num_pe_corrected": finite_mean(
                analyzed, "numPE_corrected"
            ),
            "mean_nuclear_recoil_num_pe_corrected": finite_mean(
                analyzed, "numPE_nuclear_recoil_corrected"
            ),
            "corrected_csv": str(corrected_path) if corrected_path else "",
        }
        summaries.append(summary)

    summary_path = output_dir / "dune_response_summary.csv"
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    print(f"[DUNE response] wrote {summary_path}")


if __name__ == "__main__":
    main()
