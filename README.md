# Neutron-Induced Scintillation Framework

A Geant4 and Python workflow for simulating neutron-induced scintillation in
liquid argon (LAr), scanning optical sensor coverage, calibrating the simulated
gamma/electron-recoil light yield, correcting neutron-recoil photoelectron
production, and auditing the topology of selected neutron events.

The framework has two distinct simulation modes driven by two Python entry
points:

- `scripts/calibrate_gamma_trend.py` always runs gamma calibration macros.
- `scripts/run_scan.py` always runs neutron sensor-scan macros.

The YAML file supplies parameters; it does **not** select gamma versus neutron
mode by itself.

## Workflow at a glance

```text
config/example_config.yaml
        |
        |  calibrate_gamma_trend.py
        v
output/gamma_calibration/config_with_gamma_trend.yaml
        |
        |  calcCuts (optional: choose neutron energy and TOF window)
        |  run_scan.py
        v
output/<neutron-run>/numPE_*.csv
        |
        +--> analyze_money_plot.py       corrected PE summaries and plots
        +--> analyze_event_types.py      event-topology frequencies
```

The normal sequence is:

1. Configure detector geometry, LAr optics, gamma energies, and run settings in
   `config/example_config.yaml`.
2. Build the Geant4 executable.
3. Run the gamma calibration.
4. Use the generated `config_with_gamma_trend.yaml` for neutron simulation and
   post-analysis.
5. For each neutron energy, calculate the correct TOF window and use a distinct
   output directory.
6. Generate the sensor scan, run the neutron simulation, and inspect corrected
   PE and event-topology results.

## Requirements

- Linux or another environment supported by Geant4
- CMake 3.16 or newer
- A C++ compiler compatible with the installed Geant4 version
- Geant4 with multithreading and optical physics
- Python 3.10 or newer
- Python packages:
  - `numpy`
  - `pandas`
  - `PyYAML`
  - `matplotlib`

Create a local Python environment if needed:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy pandas pyyaml matplotlib
```

## Build

With Geant4 already available in the shell:

```bash
cmake -S . -B build-mt -DCMAKE_PREFIX_PATH=/path/to/geant4-install
cmake --build build-mt -j"$(nproc)"
```

The project builds two equivalent executable names:

- `build-mt/TestEm5`
- `build-mt/LArLightSim`

The example configuration currently uses `./TestEm5`; the scripts first look
for that executable inside `run.build_dir`.

## Configuration model

### Base configuration

`config/example_config.yaml` is the editable source configuration. It defines:

- executable, build directory, event count, thread count, and output directory;
- neutron primary energy;
- inner and outer LAr geometry;
- neutron detector labels, scattering angles, and distances;
- LAr scintillation yield, yield scale, Birks constant, timing, and absorption;
- TPB and reflector properties;
- PMT and SiPM response settings;
- TOF and event-selection settings;
- gamma calibration energies and statistics;
- neutron scintillation-correction settings.

### Generated calibrated configuration

Gamma calibration writes:

```text
output/gamma_calibration/config_with_gamma_trend.yaml
```

This is a copy of the base configuration with the fitted gamma light-yield
coefficients inserted and neutron correction enabled. It is the natural input
for subsequent neutron runs and post-analysis.

Do not copy the fitted coefficients manually unless deliberately maintaining a
separate configuration. The generated file is already reusable.

### Selecting a calibrated configuration

The Bash helpers described below use:

```bash
${LARSIM_CALIBRATED_CONFIG:-$LARSIM_DIR/output/gamma_calibration/config_with_gamma_trend.yaml}
```

Thus, the standard location works without setting another variable. Set
`LARSIM_CALIBRATED_CONFIG` only to use a different calibration or a derived
per-energy config:

```bash
export LARSIM_CALIBRATED_CONFIG="$LARSIM_DIR/output/gamma_calibration/neutron_5MeV.yaml"
```

Return to the default with:

```bash
unset LARSIM_CALIBRATED_CONFIG
```

## Sensor coverage scan

The sensor scan CSV controls which optical configurations are simulated:

```text
config/sensor_scan.csv
```

Regenerate it after editing `scripts/generate_sensor_scan.py`:

```bash
python3 scripts/generate_sensor_scan.py
```

Inspect it with:

```bash
column -s, -t config/sensor_scan.csv
```

Generate neutron macros without running them:

```bash
python3 scripts/run_scan.py \
    --config output/gamma_calibration/config_with_gamma_trend.yaml \
    --dry-run
```

`run_scan.py` regenerates every macro before execution, so previously generated
macros do not control a later run. The current `run.n_events` value is rendered
as `/run/beamOn` each time.

## Gamma calibration

### Purpose

Gamma calibration establishes the simulated electron-recoil light-yield line:

```text
N_ph,gamma = a * E_dep + b
```

The gamma simulation is sensor-independent. It records produced LAr
scintillation photons, not detected photoelectrons.

During calibration:

- the primary is changed to `gamma`;
- the source position comes from
  `analysis.gamma_trend_calibration.source_position_cm`;
- scintillation is enabled in both inner and outer LAr;
- energy deposition is accumulated across inner and outer LAr;
- original LAr scintillation photons are counted once at creation;
- TPB wavelength-shifted photons are excluded from the production count;
- `numPE` is not used by the calibration;
- only full-absorption events enter the fitted points by default.

The normal neutron source remains unchanged at `(-300, 0, 0) cm`, directed
along `+x`. The gamma macro overrides only the gamma source position.

### Relevant YAML

```yaml
analysis:
  gamma_trend_calibration:
    energies_keV: [122.0, 511.0, 662.0, 1274.0]
    n_events_per_energy: 10000
    n_threads: 11
    source_position_cm: [0.0, 0.0, 0.0]
    minimum_deposit_fraction: 0.99
    output_dir: output/gamma_calibration
```

Gamma calibration inherits the LAr optical/material values from the same YAML,
including `birks_constant_mm_per_MeV`. If Birks or another light-production
setting changes, rerun the gamma calibration and use the newly generated config
for neutron work.

The global electron-recoil yield should remain at full scale. Nuclear recoils
are reduced separately during photon creation:

```yaml
optics:
  lar:
    scintillation_yield_per_MeV: 51300.0
    scintillation_yield_scale: 1.0
    birks_constant_mm_per_MeV: 0.03
    ion_scintillation_yield_scale: 0.3
    outer_scintillation: false
```

Geant4 11 makes its built-in particle-dependent scintillation mode mutually
exclusive with Birks saturation. This framework therefore retains the normal
Birks calculation and statistically keeps each LAr scintillation photon from
an ion recoil with probability `ion_scintillation_yield_scale`. This occurs
before photon transport and sensor detection. Gamma-induced electron-recoil
photons are not thinned.

Outer-LAr scintillation is disabled for neutron runs, so neutron `numPhotons`,
`eDep`, and `numPE` describe inner-cell scintillation only. Gamma calibration
enables both identical LAr volumes and fits total created photons against total
deposited energy; that material light-yield relation is then evaluated at the
neutron event's inner-cell `eDep`.

### Run

```bash
python3 scripts/calibrate_gamma_trend.py \
    --config config/example_config.yaml
```

Useful overrides:

```bash
# Short test
python3 scripts/calibrate_gamma_trend.py \
    --config config/example_config.yaml \
    --events 100

# Explicit energies and threads
python3 scripts/calibrate_gamma_trend.py \
    --config config/example_config.yaml \
    --energies-keV 122 511 662 1274 \
    --events 10000 \
    --threads 8

# Generate macros only
python3 scripts/calibrate_gamma_trend.py \
    --config config/example_config.yaml \
    --dry-run

# Refit existing gamma CSVs without rerunning Geant4
python3 scripts/calibrate_gamma_trend.py \
    --config config/example_config.yaml \
    --fit-only
```

Use `--build` when C++ has changed. Runtime-only YAML changes such as Birks do
not require compilation.

### Thread safety and progress

Gamma events run in parallel. A Geant4 mutex protects the complete shared CSV
append operation. An atomic counter tracks completed primary gamma events.
Geant4 writes progress checkpoints every 100 completed primaries, while the
Python driver reports the latest checkpoint every five seconds and prints the
exact requested total after successful completion.

Detailed Geant4 output is stored per energy in `geant4.log`; calibration runs do
not open a visualization window.

### Outputs

```text
output/gamma_calibration/
├── config_used.yaml
├── config_with_gamma_trend.yaml
├── gamma_calibration_points.csv
├── gamma_trend.png
├── gamma_trend.yaml
├── macros/
└── runs/
    ├── 122keV/
    │   ├── gamma_events.csv
    │   └── geant4.log
    └── ...
```

Each gamma event CSV contains:

```text
eventID,primaryEnergy_MeV,totalLArEDep_MeV,totalLArNumPhotons
```

At each energy, valid events must have positive total-LAr energy deposition and
photon production. The full-absorption selection additionally requires:

```text
totalLArEDep_MeV >= minimum_deposit_fraction * primaryEnergy_MeV
```

The calibration point is the mean deposited energy and mean photon count among
those selected events. A mean is required because Geant4 fluctuates
scintillation production event by event.

## Neutron energy and TOF cuts

`calcCuts` calculates the kinematic TOF range for single elastic neutron
scattering and can update a YAML configuration:

```bash
./calcCuts ENERGY_MEV DISTANCE_CM \
    --update-config path/to/config.yaml \
    --output-dir output/unique_run_name
```

Example for 5 MeV neutrons and 100 cm detector distance:

```bash
./calcCuts 5 100 \
    --update-config output/gamma_calibration/config_with_gamma_trend.yaml \
    --output-dir output/neutron_5MeV
```

This updates:

- `beam.primary_energy_MeV`;
- `selection.tof_window_ns`;
- `run.output_dir` when `--output-dir` is supplied.

Inline and block-list YAML forms of `tof_window_ns` are equivalent and both are
supported.

For reproducibility across several neutron energies, derive one config per
energy instead of repeatedly modifying the calibrated master:

```bash
cp output/gamma_calibration/config_with_gamma_trend.yaml \
   output/gamma_calibration/neutron_5MeV.yaml

./calcCuts 5 100 \
    --update-config output/gamma_calibration/neutron_5MeV.yaml \
    --output-dir output/neutron_5MeV
```

Use a distinct output directory for every neutron energy to avoid replacing
earlier CSVs and summaries.

## Neutron simulation

Run a calibrated sensor scan:

```bash
python3 scripts/run_scan.py \
    --config output/gamma_calibration/config_with_gamma_trend.yaml
```

`run_scan.py`:

1. reads `config/sensor_scan.csv`;
2. regenerates one macro per optical configuration;
3. copies the YAML and scan CSV into the run output directory;
4. runs every macro with `run.n_threads`;
5. writes neutron `numPE_*.csv` files into `run.output_dir`;
6. automatically invokes corrected PE analysis when
   `analysis.make_plots: true`.

The neutron CSV rows are written after the detected-neutron, positive-PE, and
TOF-window conditions pass. Important columns include:

- `numPE`
- `numPhotons`
- `eDep`
- `numElasticSensitive`
- `numInelastic`
- `nCapture`
- `tof`
- `detector`
- `nucleusRecoilEnergy`
- `scatteredNotSensitive`
- `InelasticSensitive`

## Neutron scintillation correction

For newly simulated data, gamma calibration enables this normalization:

```text
E_dep               = innerCell energy deposit from the event CSV
N_ph,expected,gamma = a * E_dep,MeV + b
correction_factor   = N_ph,expected,gamma / N_ph,simulated
numPE_corrected     = numPE_raw * correction_factor
```

This is evaluated separately for every event using inner-cell `eDep`,
`numPhotons`, and `numPE`. Only after each event is corrected are values grouped
and averaged for plots. The formula does not multiply by Leff because the `0.3`
ion efficiency was already sampled during photon production.

### Legacy post-analysis correction

The previous deterministic correction is retained only for reproducing older
runs:

```text
E_dep               = innerCell energy deposit from the event CSV
N_ph,expected,gamma = a * E_dep,MeV + b
N_ph,expected,NR    = N_ph,expected,gamma * Leff(E_dep)
correction_factor   = N_ph,expected,NR / N_ph,simulated
numPE_corrected     = numPE_raw * correction_factor
```

Do not select `legacy_post` for newly simulated data using
`ion_scintillation_yield_scale`; doing so would apply Leff twice. Gamma
calibration writes the fitted line into the generated config, selects
`mode: simulation_native`, and enables the gamma-only normalization.

For reproducing an older unquenched simulation, the legacy constant efficiency
is still accepted:

```yaml
leff:
  constant: 0.3
```

An energy-dependent form is also supported by removing `constant` and providing
at least two `[energy, Leff]` points plus an extrapolation policy (`error`,
`nan`, `clip`, or `linear`).

The gamma relation may be supplied either as fitted coefficients:

```yaml
gamma_calibration:
  slope_photons_per_MeV: 51300.0
  intercept_photons: 0.0
```

or as at least two energy/photon points. If a slope is present, coefficients
take precedence.

## PE post-analysis

Standalone analysis:

```bash
python3 scripts/analyze_money_plot.py \
    --config output/gamma_calibration/config_with_gamma_trend.yaml
```

When `--input-dir` and `--output-dir` are omitted, both default to the selected
config's `run.output_dir`.

Bootstrap example:

```bash
python3 scripts/analyze_money_plot.py \
    --config output/gamma_calibration/config_with_gamma_trend.yaml \
    --uncertainty bootstrap \
    --bootstrap-samples 10000 \
    --bootstrap-confidence 0.68
```

Outputs include:

- `run_summary.csv`
- `run_summary_by_detector.csv`
- combined coverage-versus-PE plots
- one coverage plot per neutron detector
- a separate true-single-elastic mean-PE series on every money plot
- raw-versus-corrected event-level PE histograms under one subdirectory per
  neutron detector and one PNG per optical configuration
- raw-versus-corrected PE plots by detector angle
- mean-PE plots with event-to-event sample standard deviation error bars
- correction diagnostics such as raw/corrected mean PE and mean correction
  factor.

## Event-topology frequency analysis

Analyze event types among rows already written after the simulation cuts:

```bash
python3 scripts/analyze_event_types.py \
    output/neutron_2MeV \
    output/neutron_14MeV
```

Categories are mutually exclusive:

- **True single elastic:** exactly one sensitive-LAr elastic scatter, no
  `scatteredNotSensitive` flag, no inelastic interaction, and no capture.
- **Single elastic + external scatter:** exactly one sensitive-LAr elastic
  scatter with `scatteredNotSensitive`, no inelastic interaction, and no
  capture.
- **Multiple elastic:** at least two sensitive-LAr elastic scatters, with or
  without an external scatter, and no inelastic interaction or capture.
- **Other / inelastic / capture:** every remaining passed row.

`scatteredNotSensitive` is the combined available flag for scattering outside
the sensitive LAr, including external, cryostat, inner-layer, and external
inelastic activity.

Outputs under `output/event_type_analysis/` include:

- `event_type_frequencies.csv` with file-, run-, and detector-level counts;
- `event_type_frequencies_by_run.png`;
- one detector/angle frequency plot per neutron run.

Detector plots are ordered and labeled by the scattering angles archived in
each run's `config_used.yaml`. Every bar is annotated with `N=<passed events>`.

## Recommended Bash helpers

Add the following to `~/.bashrc`, adjusting `LARSIM_DIR` and `G4INSTALL` for the
local installation:

```bash
# ============================================================
# LAr / Geant4 simulation helpers
# ============================================================

export LARSIM_DIR="$HOME/projects/neutron-induced-scintillation-framework"
export G4INSTALL="$HOME/Downloads/geant4-v11.0.2-install"

# Optional override. If unset, calibrated helpers use:
# $LARSIM_DIR/output/gamma_calibration/config_with_gamma_trend.yaml
# export LARSIM_CALIBRATED_CONFIG="$HOME/path/to/config_with_gamma_trend.yaml"

larsim-env() {
    source "$G4INSTALL/bin/geant4.sh"
    cd "$LARSIM_DIR" || return
    if [ -d "$LARSIM_DIR/.venv" ]; then
        source "$LARSIM_DIR/.venv/bin/activate"
    fi
}

larsim-build() {
    larsim-env || return
    cmake -S . -B build-mt -DCMAKE_PREFIX_PATH="$G4INSTALL"
    cmake --build build-mt -j"$(nproc)"
}

larsim-clean-build() {
    larsim-env || return
    rm -rf build-mt
    cmake -S . -B build-mt -DCMAKE_PREFIX_PATH="$G4INSTALL"
    cmake --build build-mt -j"$(nproc)"
}

larsim-dry() {
    larsim-env || return
    python3 scripts/run_scan.py \
        --config config/example_config.yaml \
        --dry-run "$@"
}

larsim-run() {
    larsim-env || return
    python3 scripts/run_scan.py \
        --config config/example_config.yaml "$@"
}

larsim-macros() {
    larsim-env || return
    python3 scripts/generate_sensor_scan.py "$@"
    python3 scripts/make_macros.py \
        --config config/example_config.yaml
}

larsim-gamma-macros() {
    larsim-env || return
    python3 scripts/calibrate_gamma_trend.py \
        --config config/example_config.yaml \
        --dry-run "$@"
}

larsim-gamma-run() {
    larsim-env || return
    python3 scripts/calibrate_gamma_trend.py \
        --config config/example_config.yaml "$@"
}

larsim-gamma-fit() {
    larsim-env || return
    python3 scripts/calibrate_gamma_trend.py \
        --config config/example_config.yaml \
        --fit-only "$@"
}

larsim-neutron-run() {
    larsim-env || return
    local config="${LARSIM_CALIBRATED_CONFIG:-$LARSIM_DIR/output/gamma_calibration/config_with_gamma_trend.yaml}"
    if [ ! -f "$config" ]; then
        echo "Missing calibrated config: $config" >&2
        echo "Run larsim-gamma-run or set LARSIM_CALIBRATED_CONFIG." >&2
        return 1
    fi
    python3 scripts/run_scan.py --config "$config" "$@"
}

larsim-neutron-dry() {
    larsim-env || return
    local config="${LARSIM_CALIBRATED_CONFIG:-$LARSIM_DIR/output/gamma_calibration/config_with_gamma_trend.yaml}"
    if [ ! -f "$config" ]; then
        echo "Missing calibrated config: $config" >&2
        return 1
    fi
    python3 scripts/run_scan.py --config "$config" --dry-run "$@"
}

larsim-post() {
    larsim-env || return
    local config="${LARSIM_CALIBRATED_CONFIG:-$LARSIM_DIR/output/gamma_calibration/config_with_gamma_trend.yaml}"
    if [ ! -f "$config" ]; then
        echo "Missing calibrated config: $config" >&2
        echo "Run larsim-gamma-run or set LARSIM_CALIBRATED_CONFIG." >&2
        return 1
    fi
    python3 scripts/analyze_money_plot.py --config "$config" "$@"
}

larsim-event-types() {
    larsim-env || return
    python3 scripts/analyze_event_types.py "$@"
}

larsim-test() {
    larsim-env || return
    ./build-mt/LArLightSim vis.mac 1
}
```

Reload after editing:

```bash
source ~/.bashrc
```

### Helper summary

| Helper | Purpose |
|---|---|
| `larsim-env` | Source Geant4, enter the repository, and activate `.venv`. |
| `larsim-build` | Configure and incrementally build the Geant4 executables. |
| `larsim-clean-build` | Remove `build-mt` and perform a clean build. |
| `larsim-dry` | Generate base-config neutron macros without running them. |
| `larsim-run` | Run a neutron scan directly from `example_config.yaml` without requiring gamma calibration. |
| `larsim-macros` | Regenerate `sensor_scan.csv` and base-config macros. |
| `larsim-gamma-macros` | Generate gamma calibration macros only. |
| `larsim-gamma-run` | Run gamma simulations, fit the trend, and create the calibrated YAML. |
| `larsim-gamma-fit` | Refit existing gamma event CSVs without rerunning Geant4. |
| `larsim-neutron-dry` | Generate neutron macros using the selected calibrated YAML. |
| `larsim-neutron-run` | Run the calibrated neutron scan and automatic corrected post-analysis. |
| `larsim-post` | Rerun corrected PE analysis on existing neutron CSVs. |
| `larsim-event-types` | Calculate event-topology frequencies for one or more run directories. |
| `larsim-test` | Start the visualization/test macro with one thread. |

## Practical end-to-end example

```bash
# 1. Enter environment and build after source-code changes
larsim-build

# 2. Configure gamma settings and Birks in example_config.yaml, then calibrate
larsim-gamma-run

# 3. Regenerate the optical coverage scan after editing its generator
python3 scripts/generate_sensor_scan.py

# 4. Make a reproducible config for a 5 MeV neutron run
cp output/gamma_calibration/config_with_gamma_trend.yaml \
   output/gamma_calibration/neutron_5MeV.yaml

./calcCuts 5 100 \
    --update-config output/gamma_calibration/neutron_5MeV.yaml \
    --output-dir output/neutron_5MeV

# 5. Select that calibrated config and inspect generated macros
export LARSIM_CALIBRATED_CONFIG="$LARSIM_DIR/output/gamma_calibration/neutron_5MeV.yaml"
larsim-neutron-dry

# 6. Run all sensor configurations; corrected PE analysis runs afterward
larsim-neutron-run

# 7. Rerun corrected PE analysis if desired
larsim-post --uncertainty bootstrap --bootstrap-samples 10000

# 8. Audit event topology by detector/scattering angle
larsim-event-types output/neutron_5MeV
```

## Reproducibility notes

- Every neutron run archives `config_used.yaml` and `sensor_scan_used.csv` in
  its output directory.
- Every gamma calibration archives `config_used.yaml` and produces its own
  macros, raw event CSVs, logs, fitted points, and reusable calibrated config.
- Keep one output directory and preferably one derived YAML per neutron energy.
- Gamma calibration coefficients are tied to the configured electron-recoil
  LAr scintillation model, including yield scale and Birks constant. The ion
  scale affects neutron photon production but not gamma calibration.
- YAML inline and block sequence styles are semantically equivalent.
- Generated files under `output/`, CSVs, logs, build directories, and generated
  macros are ignored by default through `.gitignore`.

## Troubleshooting

### Gamma run appears silent

Gamma Geant4 output is redirected to an energy-specific `geant4.log`. The main
driver prints periodic primary-completion checkpoints. Follow the full log with:

```bash
tail -f output/gamma_calibration/runs/122keV/geant4.log
```

### Calibrated config is missing

Run `larsim-gamma-run`, or point helpers to another file:

```bash
export LARSIM_CALIBRATED_CONFIG=/absolute/path/to/config_with_gamma_trend.yaml
```

### Changed event count is not reflected

`larsim-neutron-run` uses the calibrated config, not `example_config.yaml`.
Change `run.n_events` in the selected calibrated config, then use
`larsim-neutron-dry` to verify `/run/beamOn` in `macros/generated/*.mac`.

### Matplotlib cache warning

If the default cache is not writable:

```bash
export MPLCONFIGDIR="$LARSIM_DIR/.matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
```
