#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def build_scan_rows():
    configs = []

    # 1-4 top PMTs, no bottom PMTs or SiPM rows
    configs.extend((top, 0, 0) for top in range(1, 5))

    # Four top PMTs, then 1-4 bottom PMTs
    configs.extend((4, bottom, 0) for bottom in range(1, 5))

    # Four top + four bottom PMTs, then 1-10 SiPM rows
    configs.extend((4, 4, rows) for rows in range(1, 11))

    rows = []
    for i, (top, bottom, sipm_rows) in enumerate(configs):
        rows.append({
            "config_id": f"C{i:02d}",
            "top_pmts": top,
            "bottom_pmts": bottom,
            "sipm_rows": sipm_rows,
            "label": f"top{top}_bot{bottom}_sipm{sipm_rows}",
        })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Generate the thesis-style optical-coverage scan CSV."
    )
    parser.add_argument(
        "--output",
        default="config/sensor_scan.csv",
        help="Output CSV path.",
    )
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = build_scan_rows()
    fieldnames = [
        "config_id",
        "top_pmts",
        "bottom_pmts",
        "sipm_rows",
        "label",
    ]

    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} configurations to {output}")


if __name__ == "__main__":
    main()
