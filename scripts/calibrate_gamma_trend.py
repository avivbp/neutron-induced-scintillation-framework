#!/usr/bin/env python3
"""Run full-absorption gamma calibrations and fit N_ph = a E_dep + b.

The generated ``config_with_gamma_trend.yaml`` is a copy of the input config
with the fitted line inserted and neutron scintillation correction enabled.
"""

from __future__ import annotations

import argparse
import copy
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from config_tools import flatten_for_template, read_yaml, render_template


def run(command: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=str(cwd) if cwd else None, check=True)


def run_logged(
    command: list[str],
    cwd: Path,
    log_path: Path,
    requested_events: int,
) -> None:
    print("+", " ".join(command), f"(log: {log_path})")
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        started = time.monotonic()
        while process.poll() is None:
            time.sleep(5.0)
            try:
                log.flush()
                log_text = log_path.read_text(encoding="utf-8", errors="replace")
                matches = re.findall(
                    r"\[gamma progress\] primaries completed: (\d+)",
                    log_text,
                )
                completed = int(matches[-1]) if matches else 0
            except OSError:
                completed = 0
            elapsed = time.monotonic() - started
            percent = 100.0 * completed / requested_events
            print(
                "[gamma calibration] Geant4 running: "
                f"{elapsed:.0f} s elapsed; {completed}/{requested_events} "
                f"primaries completed ({percent:.1f}%)",
                flush=True,
            )
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)
        elapsed = time.monotonic() - started
        print(
            "[gamma calibration] Geant4 complete: "
            f"{requested_events}/{requested_events} primaries completed "
            f"in {elapsed:.0f} s",
            flush=True,
        )


def executable_path(root: Path, config: dict) -> Path:
    build_dir = root / config.get("run", {}).get("build_dir", "build")
    executable = Path(config.get("run", {}).get("executable", "./LArLightSim"))
    if executable.is_absolute():
        return executable
    build_candidate = build_dir / executable.name
    return build_candidate if build_candidate.exists() else root / executable


def calibration_settings(config: dict, args: argparse.Namespace) -> dict:
    configured = config.get("analysis", {}).get("gamma_trend_calibration", {})
    if not isinstance(configured, dict):
        raise ValueError("analysis.gamma_trend_calibration must be a mapping")

    energies = args.energies_keV or configured.get(
        "energies_keV", [122.0, 511.0, 662.0, 1274.0]
    )
    energies = sorted({float(value) for value in energies})
    if len(energies) < 2 or any(value <= 0 for value in energies):
        raise ValueError("at least two positive gamma energies are required")

    source = args.source_position_cm or configured.get(
        "source_position_cm", [0.0, 0.0, 0.0]
    )
    if len(source) != 3:
        raise ValueError("source_position_cm must contain x, y, z")

    return {
        "energies_keV": energies,
        "n_events": int(
            args.events
            if args.events is not None
            else configured.get("n_events_per_energy", 1000)
        ),
        "n_threads": int(
            args.threads
            if args.threads is not None
            else configured.get(
                "n_threads", config.get("run", {}).get("n_threads", 1)
            )
        ),
        "source_position_cm": [float(value) for value in source],
        "minimum_deposit_fraction": float(
            args.minimum_deposit_fraction
            if args.minimum_deposit_fraction is not None
            else configured.get("minimum_deposit_fraction", 0.99)
        ),
    }


def energy_tag(energy_keV: float) -> str:
    return f"{energy_keV:g}".replace(".", "p") + "keV"


def generate_and_run(
    root: Path,
    config: dict,
    settings: dict,
    output_dir: Path,
    dry_run: bool,
) -> None:
    template_path = root / "macros/template_gamma_calibration.mac"
    template = template_path.read_text(encoding="utf-8")
    macro_dir = output_dir / "macros"
    macro_dir.mkdir(parents=True, exist_ok=True)
    executable = executable_path(root, config)

    row = {
        "config_id": "gamma_calibration",
        "label": "gamma_calibration",
        # Required by the shared template flattener but intentionally unused
        # by the gamma calibration macro: this calibration is sensor-agnostic.
        "top_pmts": 0,
        "bottom_pmts": 0,
        "sipm_rows": 0,
    }

    for energy_keV in settings["energies_keV"]:
        tag = energy_tag(energy_keV)
        run_dir = output_dir / "runs" / tag
        run_dir.mkdir(parents=True, exist_ok=True)
        event_path = run_dir / "gamma_events.csv"
        event_path.write_text(
            "eventID,primaryEnergy_MeV,totalLArEDep_MeV,"
            "totalLArNumPhotons\n",
            encoding="utf-8",
        )

        values = flatten_for_template(config, row)
        values.update(
            {
                "gamma_energy_MeV": energy_keV / 1000.0,
                "gamma_source_x_cm": settings["source_position_cm"][0],
                "gamma_source_y_cm": settings["source_position_cm"][1],
                "gamma_source_z_cm": settings["source_position_cm"][2],
                "n_events": settings["n_events"],
            }
        )
        macro_path = macro_dir / f"gamma_{tag}.mac"
        macro_path.write_text(render_template(template, values), encoding="utf-8")
        if not dry_run:
            run_logged(
                [str(executable), str(macro_path), str(settings["n_threads"])],
                cwd=run_dir,
                log_path=run_dir / "geant4.log",
                requested_events=settings["n_events"],
            )


def fit_points(settings: dict, output_dir: Path) -> tuple[pd.DataFrame, float, float, float]:
    rows: list[dict] = []
    fraction = settings["minimum_deposit_fraction"]

    for energy_keV in settings["energies_keV"]:
        event_path = output_dir / "runs" / energy_tag(energy_keV) / "gamma_events.csv"
        if not event_path.exists():
            print(f"[gamma calibration] warning: missing {event_path}")
            continue
        frame = pd.read_csv(event_path)
        for column in ("totalLArEDep_MeV", "totalLArNumPhotons"):
            if column not in frame.columns:
                raise ValueError(f"{event_path} is missing column {column}")
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

        valid = frame[
            np.isfinite(frame["totalLArEDep_MeV"])
            & np.isfinite(frame["totalLArNumPhotons"])
            & (frame["totalLArEDep_MeV"] > 0)
            & (frame["totalLArNumPhotons"] > 0)
        ].copy()
        selected = valid[
            valid["totalLArEDep_MeV"] >= fraction * energy_keV / 1000.0
        ]
        if selected.empty:
            print(
                "[gamma calibration] warning: no selected deposits for "
                f"{energy_keV:g} keV"
            )
            continue

        n = len(selected)
        photon_std = (
            float(selected["totalLArNumPhotons"].std(ddof=1))
            if n > 1 else 0.0
        )
        rows.append(
            {
                "primary_energy_keV": energy_keV,
                "n_events_requested": settings["n_events"],
                "n_events_with_deposit": len(valid),
                "n_events_selected": n,
                "mean_totalLArEDep_MeV": float(
                    selected["totalLArEDep_MeV"].mean()
                ),
                "std_totalLArEDep_MeV": float(
                    selected["totalLArEDep_MeV"].std(ddof=1)
                ) if n > 1 else 0.0,
                "mean_totalLArNumPhotons": float(
                    selected["totalLArNumPhotons"].mean()
                ),
                "std_totalLArNumPhotons": photon_std,
                "sem_totalLArNumPhotons": photon_std / math.sqrt(n),
            }
        )

    points = pd.DataFrame(rows)
    if len(points) < 2:
        raise ValueError(
            "fewer than two gamma energies had selected events; increase the "
            "event count or lower minimum_deposit_fraction"
        )

    x = points["mean_totalLArEDep_MeV"].to_numpy(dtype=float)
    y = points["mean_totalLArNumPhotons"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, deg=1)
    predicted = slope * x + intercept
    residual_sum = float(np.sum((y - predicted) ** 2))
    total_sum = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0 else 1.0
    points["fit_numPhotons"] = predicted
    points["residual_numPhotons"] = y - predicted
    return points, float(slope), float(intercept), r_squared


def write_results(
    config: dict,
    settings: dict,
    output_dir: Path,
    points: pd.DataFrame,
    slope: float,
    intercept: float,
    r_squared: float,
) -> None:
    points.to_csv(output_dir / "gamma_calibration_points.csv", index=False)
    fit = {
        "gamma_calibration": {
            "slope_photons_per_MeV": slope,
            "intercept_photons": intercept,
            "r_squared": r_squared,
            "fit_x_quantity": "mean total inner+outer LAr deposited energy",
            "event_selection": "full absorption in total LAr",
            "minimum_deposit_fraction": settings["minimum_deposit_fraction"],
        }
    }
    with (output_dir / "gamma_trend.yaml").open("w", encoding="utf-8") as stream:
        yaml.safe_dump(fit, stream, sort_keys=False)

    # Preserve the user's YAML mapping order in the reusable config. A
    # dump/load round-trip sorts keys by default and needlessly rearranges
    # detector fields such as label, angle_deg, and distance_cm.
    updated = copy.deepcopy(config)
    correction = updated.setdefault("analysis", {}).setdefault(
        "neutron_scintillation_correction", {}
    )
    correction["enabled"] = True
    correction["gamma_calibration"] = {
        "slope_photons_per_MeV": slope,
        "intercept_photons": intercept,
    }
    with (output_dir / "config_with_gamma_trend.yaml").open(
        "w", encoding="utf-8"
    ) as stream:
        yaml.safe_dump(updated, stream, sort_keys=False)

    figure, axis = plt.subplots(figsize=(7.5, 5.5))
    axis.errorbar(
        points["mean_totalLArEDep_MeV"],
        points["mean_totalLArNumPhotons"],
        yerr=points["sem_totalLArNumPhotons"],
        fmt="o",
        capsize=3,
        label="Gamma simulations",
    )
    x_line = np.linspace(
        0.0, 1.05 * points["mean_totalLArEDep_MeV"].max(), 200
    )
    axis.plot(
        x_line,
        slope * x_line + intercept,
        label=f"Nph = {slope:.6g} Edep + {intercept:.6g}\n$R^2$ = {r_squared:.6f}",
    )
    axis.set_xlabel("Mean total deposited energy in inner + outer LAr [MeV]")
    axis.set_ylabel("Mean photons for full-absorption events")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "gamma_trend.png", dpi=300)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/example_config.yaml")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--energies-keV", nargs="+", type=float)
    parser.add_argument("--events", type=int)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--source-position-cm", nargs=3, type=float)
    parser.add_argument("--minimum-deposit-fraction", type=float)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--fit-only", action="store_true")
    args = parser.parse_args()

    if args.dry_run and args.fit_only:
        parser.error("--dry-run and --fit-only cannot be combined")
    if args.events is not None and args.events <= 0:
        parser.error("--events must be positive")
    if args.threads is not None and args.threads <= 0:
        parser.error("--threads must be positive")
    if args.minimum_deposit_fraction is not None and not (
        0.0 <= args.minimum_deposit_fraction <= 1.0
    ):
        parser.error("--minimum-deposit-fraction must be between 0 and 1")

    root = Path.cwd()
    config = read_yaml(args.config)
    settings = calibration_settings(config, args)
    if not 0.0 <= settings["minimum_deposit_fraction"] <= 1.0:
        parser.error("configured minimum_deposit_fraction must be between 0 and 1")
    if settings["n_threads"] <= 0:
        parser.error("configured n_threads must be positive")
    output_dir = Path(
        args.output_dir
        or config.get("analysis", {})
        .get("gamma_trend_calibration", {})
        .get("output_dir", "output/gamma_calibration")
    )
    if not output_dir.is_absolute():
        output_dir = root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.config, output_dir / "config_used.yaml")

    if args.build:
        build_dir = root / config.get("run", {}).get("build_dir", "build")
        run(["cmake", "-S", str(root), "-B", str(build_dir)])
        run(["cmake", "--build", str(build_dir), f"-j{args.jobs}"])

    if not args.fit_only:
        generate_and_run(root, config, settings, output_dir, args.dry_run)
    if args.dry_run:
        print(f"[gamma calibration] generated macros in {output_dir / 'macros'}")
        return

    points, slope, intercept, r_squared = fit_points(settings, output_dir)
    write_results(
        config, settings, output_dir, points, slope, intercept, r_squared
    )
    print(
        "[gamma calibration] fitted trend: "
        f"N_ph = ({slope:.10g} photons/MeV) E_dep + "
        f"({intercept:.10g} photons); R^2={r_squared:.8f}"
    )
    print(
        "[gamma calibration] reusable config: "
        f"{output_dir / 'config_with_gamma_trend.yaml'}"
    )


if __name__ == "__main__":
    main()
