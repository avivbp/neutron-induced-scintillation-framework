#!/usr/bin/env python3
"""Build/run a coverage scan and then aggregate the results.

Typical use from Code/LabSimulation:
    python scripts/run_scan.py --config config/example_config.yaml --scan config/sensor_scan.csv --build
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from config_tools import flatten_for_template, read_scan_csv, read_yaml, render_template, write_json


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def maybe_build(source_dir: Path, build_dir: Path, jobs: int) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    run(["cmake", "-S", str(source_dir), "-B", str(build_dir)])
    run(["cmake", "--build", str(build_dir), f"-j{jobs}"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/example_config.yaml")
    parser.add_argument("--scan", default="config/sensor_scan.csv")
    parser.add_argument("--build", action="store_true", help="Run CMake before the scan.")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--dry-run", action="store_true", help="Generate macros/metadata but do not run Geant4.")
    args = parser.parse_args()

    root = Path.cwd()
    config = read_yaml(args.config)
    rows = read_scan_csv(args.scan)

    outdir = root / config.get("run", {}).get("output_dir", "output/scan")
    macro_dir = root / config.get("run", {}).get("generated_macro_dir", "macros/generated")
    macro_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    template_path = root / config.get("run", {}).get("macro_template", "macros/template_run.mac")
    template = template_path.read_text(encoding="utf-8")

    # Keep a copy of the inputs used for reproducibility.
    shutil.copy2(args.config, outdir / "config_used.yaml")
    shutil.copy2(args.scan, outdir / "sensor_scan_used.csv")

    build_dir = root / config.get("run", {}).get("build_dir", "build")
    if args.build:
        maybe_build(root, build_dir, args.jobs)

    executable = config.get("run", {}).get("executable", "./LArLightSim")
    exe_path = Path(executable)
    if not exe_path.is_absolute():
        # Prefer executable inside build dir, then current dir.
        build_candidate = build_dir / exe_path.name
        exe_path = build_candidate if build_candidate.exists() else root / exe_path

    metadata = []
    for row in rows:
        values = flatten_for_template(config, row)
        macro_path = macro_dir / f"{values['config_id']}.mac"
        macro_path.write_text(render_template(template, values), encoding="utf-8")
        metadata.append({"macro": str(macro_path), **values})
        if args.dry_run:
            continue
        n_threads = str(config.get("run", {}).get("n_threads", 1))
        run([str(exe_path), str(macro_path), n_threads])

    write_json(outdir / "scan_metadata.json", metadata)

    if not args.dry_run and config.get("analysis", {}).get("make_plots", True):
        run([
            "python3", "scripts/analyze_outputs.py",
            "--config", args.config,
            "--metadata", str(outdir / "scan_metadata.json"),
            "--input-dir", str(root),
            "--output-dir", str(outdir),
        ])


if __name__ == "__main__":
    main()
