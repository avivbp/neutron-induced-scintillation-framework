#!/usr/bin/env python3
"""Generate one Geant4 macro per detector-coverage configuration."""
from __future__ import annotations

import argparse
from pathlib import Path

from config_tools import flatten_for_template, read_scan_csv, read_yaml, render_template, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/example_config.yaml")
    parser.add_argument("--scan", default="config/sensor_scan.csv")
    parser.add_argument("--template", default=None)
    parser.add_argument("--outdir", default=None)
    args = parser.parse_args()

    config = read_yaml(args.config)
    rows = read_scan_csv(args.scan)
    template_path = Path(args.template or config.get("run", {}).get("macro_template", "macros/template_run.mac"))
    outdir = Path(args.outdir or config.get("run", {}).get("generated_macro_dir", "macros/generated"))
    outdir.mkdir(parents=True, exist_ok=True)

    template = template_path.read_text(encoding="utf-8")
    metadata = []
    for row in rows:
        values = flatten_for_template(config, row)
        macro_text = render_template(template, values)
        macro_path = outdir / f"{values['config_id']}.mac"
        macro_path.write_text(macro_text, encoding="utf-8")
        metadata.append({"macro": str(macro_path), **values})
        print(f"wrote {macro_path}")

    write_json(outdir / "scan_metadata.json", metadata)


if __name__ == "__main__":
    main()
