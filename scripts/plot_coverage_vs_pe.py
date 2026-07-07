#!/usr/bin/env python3
"""Plot mean PE versus optical coverage from run_summary.csv."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", nargs="?", default="output/tau_lab_example/run_summary.csv")
    parser.add_argument("--out", default="output/tau_lab_example/figures/coverage_vs_mean_pe.png")
    args = parser.parse_args()

    df = pd.read_csv(args.summary).dropna(subset=["coverage_percent", "mean_pe"])
    if df.empty:
        raise SystemExit("No plottable rows found in run_summary.csv")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 4.5))
    yerr = df["sem_pe"] if "sem_pe" in df.columns else None
    plt.errorbar(df["coverage_percent"], df["mean_pe"], yerr=yerr, marker="o", linestyle="none", capsize=3)
    for _, row in df.iterrows():
        plt.annotate(str(row["config_id"]), (row["coverage_percent"], row["mean_pe"]), xytext=(3, 3), textcoords="offset points", fontsize=8)
    plt.xlabel("Optical coverage [%]")
    plt.ylabel("Mean detected photoelectrons")
    plt.title("Detected PE yield vs optical coverage")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
