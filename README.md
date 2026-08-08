# Neutron-Induced Scintillation Framework

A Geant4 and Python workflow for simulating neutron-induced scintillation in
liquid argon (LAr). It supports both the original nested cylindrical cell and a
generic DUNE-like box cryostat, configurable optical sensors and external
neutron taggers, gamma/electron-recoil light-yield calibration, recoil-resolved
photoelectron analysis, uncut neutron-interaction truth, particle-resolved
energy/light/PE bookkeeping, and channel-aware event topology analysis.

## Capabilities

| Area | Supported capabilities |
| --- | --- |
| Detector geometry | Legacy concentric cylindrical LAr cells (`nested_cell`) or a configurable rectangular active-LAr volume, inset fiducial volume, and steel shell (`box_cryostat`). |
| Volume semantics | Logical-volume roles for active LAr, fiducial LAr, inactive LAr, cryostat, optical detector, and external detector. Analysis does not depend on legacy volume names. |
| Optical sensors | Legacy top/bottom PMT patches and cylindrical SiPM rows; box-face PMT or SiPM grids on any combination of `+x`, `-x`, `+y`, `-y`, `+z`, and `-z`. |
| Sensor response | Flat PDE, scalable built-in spectral curves, raw built-in curves, scalable CSV curves, or raw CSV curves; independent PMT/SiPM reflectivity and PDE. |
| Physics | `QGSP_BERT_HP` neutron transport, EM option 4, optical scintillation/WLS/transport, Birks saturation, and statistical neutron/alpha/nuclear-recoil light scaling. Distributed batch templates disable Cerenkov. |
| Neutron tagging | Any number of labeled external neutron detectors with configurable angle and distance, plus a configurable TOF window. |
| Interaction transport | The legacy nested cell intentionally stops an event at its first neutron-inelastic interaction. The box cryostat retains elastic, inelastic, capture, fission, secondary-neutron, and mixed cascades. |
| Truth output | Uncut interaction records for primary and secondary neutrons, including volume roles, process/channel, kinematics, local deposit, and produced secondaries. |
| Signal bookkeeping | Per-event active-LAr energy deposition, produced scintillation photons, and detected PE split by primary neutron, secondary neutron, gamma, electron/positron, proton, alpha, nuclear recoil, and other. Optical origin survives wavelength shifting. |
| Calibration and analysis | Multithreaded gamma calibration, event-level neutron light normalization, coverage-versus-PE plots, detector-angle plots, bootstrap/SEM uncertainties, and truth-defined topology frequencies. |
| Reproducibility | Generated macros, archived YAML/scan inputs, per-configuration run labels, deterministic analysis seeds, and thread-safe CSV writers. |

Independent of the chosen detector model, the framework has two execution
modes driven by two Python entry points:

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
        +--> neutron_interactions_<config_id>.csv   uncut interaction truth
        +--> particle_class_summary_<config_id>.csv uncut energy/light/PE truth
        +--> dune_event_response_<config_id>.csv    uncut ER/NR optical signal
        +--> analyze_money_plot.py                  corrected PE summaries
        +--> analyze_event_types.py                 topology frequencies
```

The normal sequence is:

1. Start from the complete nested-cell or DUNE-like example below and save it
   as a YAML file.
2. Build the Geant4 executable.
3. Run the gamma calibration.
4. Use the generated `config_with_gamma_trend.yaml` for neutron simulation and
   post-analysis.
5. For tagged `nested_cell` runs, calculate the correct TOF window. Standalone
   `box_cryostat` runs with no external taggers skip this step.
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
- `nested_cell` or `box_cryostat` detector selection;
- neutron primary energy;
- nested cylindrical dimensions or box active/fiducial/cryostat dimensions;
- neutron detector labels, scattering angles, and distances;
- LAr scintillation yield, yield scale, Birks constant, timing, and absorption;
- TPB and reflector properties;
- PMT and SiPM response settings and optional box-face layouts;
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

## Detector models

| Behavior | `nested_cell` | `box_cryostat` |
| --- | --- | --- |
| Active geometry | Inner LAr cylinder; outer concentric LAr cylinder can optionally scintillate | Rectangular active-LAr box with an optional inset fiducial box |
| Cryostat | Existing concentric cylindrical layers | Configurable stainless-steel box shell |
| Optical placement | Four possible PMTs on each endcap plus selectable cylindrical SiPM rows | Fixed face-grid layouts declared in YAML |
| Coverage variation | Rows of `config/sensor_scan.csv` activate progressively more prebuilt sensors | Change `sensors.layouts` in YAML; use one scan row per fixed box layout |
| Neutron inelastic behavior | Record the triggering interaction, then abort the event to preserve the elastic-only legacy selection | Continue the complete cascade, including secondary neutron tracks |
| Primary use | Small tagged-neutron cell and regression workflow | General/DUNE-like active-LAr studies and mixed-channel truth |

Detector construction assigns logical-volume roles used uniformly by stepping,
scintillation, truth, and analysis code: `active_lar` contributes energy and
light, `fiducial_lar` is eligible for fiducial selections, `inactive_lar` is
transport-only LAr, `cryostat` marks structure, `optical_detector` marks sensor
surfaces, and `external_detector` marks neutron taggers. A logical volume can
carry more than one role; for example, the nested inner cell is both active and
fiducial.

`box_cryostat` is a configurable DUNE-like rectangular model, not a hard-coded
official DUNE module or imported engineering geometry. Dimensions and sensor
layouts are user inputs. Very large LAr volumes with full optical transport can
be computationally expensive; start with low event counts and sparse sensor
layouts.

The default neutron gun starts at `(-300, 0, 0) cm` and points along `+x`.
Choose box dimensions that leave this source in the intended location. The
validated example below spans `x = -200...+200 cm`, placing the source 100 cm
outside its upstream face.

### Complete nested-cell setup

Save this as `config/nested_cell.yaml`. It is a complete configuration for the
legacy two-cylinder detector, a three-angle external neutron-detector array,
gamma calibration, the standard optical-coverage scan, and post-analysis.

```yaml
run:
  executable: ./LArLightSim
  build_dir: build-mt
  output_dir: output/nested_cell_5MeV
  n_events: 10000
  n_threads: 4
  macro_template: macros/template_run.mac
  generated_macro_dir: macros/generated/nested_cell

detector:
  model: nested_cell

beam:
  primary_energy_MeV: 5.0

geometry:
  inner_lar:
    radius_cm: 5.0
    height_cm: 10.0
  outer_lar:
    diameter_cm: 35.0
    height_cm: 70.0
  neutron_detectors:
    - label: A0
      angle_deg: 25.0
      distance_cm: 100.0
    - label: A1
      angle_deg: 60.0
      distance_cm: 100.0
    - label: A2
      angle_deg: 90.0
      distance_cm: 100.0

optics:
  lar:
    scintillation_yield_per_MeV: 51300.0
    scintillation_yield_scale: 1.0
    ion_scintillation_yield_scale: 0.3
    outer_scintillation: false
    absorption_length_cm: 150.0
    fast_time_ns: 7.0
    slow_time_ns: 1500.0
    fast_fraction: 0.75
    birks_constant_mm_per_MeV: 0.03
  tpb:
    efficiency: 1.0
  reflector:
    reflectivity: 0.98

sensors:
  pmt:
    side_length_cm: 2.54
    reflectivity: 0.20
    response_mode: builtin_curve
    pde: 0.266
    top_pde: 0.266
    bottom_pde: 0.266
    qe_curve_file: none
  sipm:
    tile_size_cm: 0.60
    reflectivity: 0.40
    response_mode: builtin_curve
    pde: 0.402
    qe_curve_file: none

selection:
  require_detected: true
  require_num_pe_positive: true
  tof_window_ns: [30.0, 40.0]
  reject_inelastic_sensitive: true
  reject_external_scatter: false

analysis:
  make_plots: true
  gamma_trend_calibration:
    energies_keV: [122.0, 511.0, 662.0, 1274.0]
    n_events_per_energy: 2000
    n_threads: 4
    source_position_cm: [0.0, 0.0, 0.0]
    minimum_deposit_fraction: 0.99
    output_dir: output/gamma_nested_cell
  neutron_scintillation_correction:
    enabled: false
    mode: simulation_native
    gamma_calibration:
      slope_photons_per_MeV: 51300.0
      intercept_photons: 0.0
```

Run the complete nested-cell workflow:

```bash
# Fit the gamma/electron-recoil light-yield line.
python3 scripts/calibrate_gamma_trend.py \
    --config config/nested_cell.yaml

# Preserve the calibrated master and derive a neutron-energy-specific config.
cp output/gamma_nested_cell/config_with_gamma_trend.yaml \
   output/gamma_nested_cell/nested_cell_5MeV.yaml

# All taggers in this example are 100 cm from the LAr target.
./calcCuts 5 100 \
    --update-config output/gamma_nested_cell/nested_cell_5MeV.yaml \
    --output-dir output/nested_cell_5MeV

# Use the standard progressive PMT/SiPM scan.
python3 scripts/run_scan.py \
    --config output/gamma_nested_cell/nested_cell_5MeV.yaml \
    --scan config/sensor_scan.csv

python3 scripts/analyze_event_types.py output/nested_cell_5MeV
```

In this model, `innerCell` has the `active_lar` and `fiducial_lar` roles.
`outerCell` is `inactive_lar` when `outer_scintillation: false`. Inelastic
neutron events are recorded in interaction truth and then intentionally
aborted, so selected-event analysis remains elastic-only.

### Complete DUNE-like box setup

Save this as `config/dune_like_box.yaml`. The example uses a
`400 x 400 x 1000 cm` active-LAr box, a 10 cm inset fiducial region, a 1 cm
steel shell, PMT grids on the endcaps, and SiPM grids on the four long faces.
It intentionally has no external neutron taggers or TOF selection: the optical
response is the detector signal. A ready-to-run copy is checked in at
`config/dune_like_box.yaml`.

```yaml
run:
  executable: ./LArLightSim
  build_dir: build-mt
  output_dir: output/dune_like_14MeV
  n_events: 10
  n_threads: 1
  macro_template: macros/template_run.mac
  generated_macro_dir: macros/generated/dune_like

detector:
  model: box_cryostat

beam:
  primary_energy_MeV: 14.0

geometry:
  active_lar:
    dimensions_cm: [400.0, 400.0, 1000.0]
  fiducial:
    margin_cm: [10.0, 10.0, 10.0]
  cryostat:
    thickness_cm: 1.0
  neutron_detectors: []

optics:
  lar:
    scintillation_yield_per_MeV: 51300.0
    scintillation_yield_scale: 1.0
    ion_scintillation_yield_scale: 0.3
    outer_scintillation: false
    absorption_length_cm: 150.0
    fast_time_ns: 7.0
    slow_time_ns: 1500.0
    fast_fraction: 0.75
    birks_constant_mm_per_MeV: 0.03
  tpb:
    efficiency: 1.0
  reflector:
    reflectivity: 0.98

sensors:
  pmt:
    side_length_cm: 20.0
    active_size_cm: [20.0, 20.0]
    reflectivity: 0.20
    response_mode: builtin_curve
    pde: 0.266
    top_pde: 0.266
    bottom_pde: 0.266
    qe_curve_file: none
  sipm:
    tile_size_cm: 10.0
    active_size_cm: [10.0, 10.0]
    reflectivity: 0.40
    response_mode: builtin_curve
    pde: 0.402
    qe_curve_file: none
  layouts:
    - id: long_wall_sipms
      type: sipm
      faces: [+x, -x, +y, -y]
      pattern: grid
      active_size_cm: [10.0, 10.0]
      thickness_cm: 0.10
      pitch_cm: [100.0, 100.0]
      edge_clearance_cm: 25.0
      inset_cm: 0.10
    - id: endcap_pmts
      type: pmt
      faces: [+z, -z]
      pattern: grid
      active_size_cm: [20.0, 20.0]
      thickness_cm: 0.10
      pitch_cm: [100.0, 100.0]
      edge_clearance_cm: 25.0
      inset_cm: 0.10

analysis:
  make_plots: true
  gamma_trend_calibration:
    energies_keV: [122.0, 511.0, 662.0, 1274.0]
    n_events_per_energy: 2000
    n_threads: 4
    source_position_cm: [0.0, 0.0, 0.0]
    minimum_deposit_fraction: 0.99
    output_dir: output/gamma_dune_like
  neutron_scintillation_correction:
    enabled: false
    mode: simulation_native
    gamma_calibration:
      slope_photons_per_MeV: 51300.0
      intercept_photons: 0.0
```

Box sensor layouts are constructed before the run and all declared tiles are
active. The legacy `top_pmts`, `bottom_pmts`, and `sipm_rows` scan columns do
not enable or disable box tiles. For this fixed layout, save a one-row scan as
`config/dune_like_scan.csv`:

```csv
config_id,top_pmts,bottom_pmts,sipm_rows,label
DUNE00,0,0,0,dune_like_fixed_layout
```

Run the complete DUNE-like workflow:

```bash
# Fast uncalibrated smoke run: 10 events, one thread, no tagger or TOF.
python3 scripts/run_scan.py \
    --config config/dune_like_box.yaml \
    --scan config/dune_like_scan.csv \
    --build

# Optional gamma normalization for a larger analysis run.
python3 scripts/calibrate_gamma_trend.py \
    --config config/dune_like_box.yaml

python3 scripts/run_scan.py \
    --config output/gamma_dune_like/config_with_gamma_trend.yaml \
    --scan config/dune_like_scan.csv
```

The box model never applies the legacy inelastic-event abort. Its interaction
truth and particle-class summary therefore include the complete inelastic
cascade and secondary-neutron transport. `dune_event_response_DUNE00.csv`
contains one row for every event regardless of PE, interaction type, external
tagger, or TOF. `run_scan.py` summarizes it with
`scripts/analyze_dune_response.py`.

### Box layout field reference

| Field | Meaning |
| --- | --- |
| `id` | Stable layout name containing letters, numbers, `_`, or `-`. |
| `type` | `pmt` or `sipm`; selects the corresponding response model. |
| `faces` | One or more of `+x`, `-x`, `+y`, `-y`, `+z`, `-z`. |
| `pattern` | Currently `grid`. |
| `active_size_cm` | Two-dimensional active footprint of one tile. |
| `thickness_cm` | Tile depth inside the active-LAr host volume. |
| `pitch_cm` | Center-to-center spacing along the two local face axes; each pitch must be at least the corresponding active size. |
| `edge_clearance_cm` | Minimum unused border around each selected face. |
| `inset_cm` | Offset inward from the active-LAr boundary. With a separate fiducial box, its margin on that face must be at least `thickness_cm + inset_cm`. |

### Sensor response modes

Both `sensors.pmt` and `sensors.sipm` accept the same response modes:

| Mode | Behavior |
| --- | --- |
| `flat` | Constant detection efficiency set by `pde` (or PMT `top_pde`/`bottom_pde`). |
| `builtin_curve` | Built-in wavelength response rescaled so its peak equals the requested PDE. |
| `builtin_curve_raw` | Built-in wavelength response used without peak rescaling. |
| `csv_curve` | User CSV response rescaled to the requested peak PDE. |
| `csv_curve_raw` | User CSV response used as an absolute efficiency curve. |

For CSV modes, set `qe_curve_file` to a file containing either
`wavelength_nm,qe` or `energy_eV,qe`; `qe` may be a fraction or percentage.
Reflectivity remains independently configurable. The physical PDE is converted
to Geant4 boundary efficiency so reflectivity and detection retain the intended
probabilities.

At simulation time, the legacy `numPE_*.csv` writer always requires an
external neutron detection, positive PE, and a TOF inside
`selection.tof_window_ns`. The other `selection` booleans are consumed by the
older `scripts/analyze_outputs.py` workflow; they do not loosen or tighten the
Geant4 CSV writer. The box-specific `dune_event_response_<config_id>.csv` and
the uncut interaction and particle summaries are independent of all these
selected-event conditions.

## Sensor coverage scan

For `nested_cell`, the sensor scan CSV controls which of the prebuilt top PMTs,
bottom PMTs, and cylindrical SiPM rows are active:

```text
config/sensor_scan.csv
```

Regenerate the standard nested-cell scan after editing
`scripts/generate_sensor_scan.py`:

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

For `box_cryostat`, sensor geometry comes entirely from `sensors.layouts` and
all declared tiles are active. Use a one-row scan CSV for each fixed YAML
layout, as shown in the DUNE-like example. The row's `config_id` names
`neutron_interactions_<config_id>.csv` and
`particle_class_summary_<config_id>.csv`; its three legacy sensor-count fields
remain in the selected-event filename but do not control box tiles. To compare
several box coverages, use separate YAML/layout and output-directory pairs.

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
- all relevant LAr is active (`nested_cell` enables both cylinders; the box
  uses its active-LAr volume);
- energy deposition is accumulated across all active LAr;
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
Birks calculation and statistically keeps each LAr scintillation photon whose
parent is a neutron, alpha, or heavier nuclear recoil with probability
`ion_scintillation_yield_scale`. This preserves the legacy proton behavior and
occurs before photon transport and sensor detection. Gamma-induced
electron-recoil photons are not thinned.

For the nested neutron model, outer-LAr scintillation is normally disabled, so
`numPhotons`, `eDep`, and `numPE` describe the inner cell. Gamma calibration
temporarily enables both nested LAr cylinders. For the box model, these values
describe all volumes carrying the `active_lar` role. In either geometry, the
calibration fits total created photons against active-LAr deposited energy and
the neutron correction evaluates that relation using the event's active-LAr
`eDep`.

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

Every neutron interaction is also written independently of those event cuts to
`neutron_interactions_<config_id>.csv`. This includes interactions by primary
and secondary neutron tracks. The channel column uses `elastic`, `inelastic`,
`capture`, `fission`, `other_hadronic`, `transport_exit`, or `unclassified`;
the original Geant4 process name is retained alongside it.

Each interaction row contains the event, interaction, track, parent, and step
identifiers; physical volume and configured volume roles; position and time;
neutron kinetic energy before and after the step; local energy deposit; and all
secondaries created by the interaction. The `secondaries` field is a
pipe-separated list of `particle:PDG-code:kinetic-energy-keV` entries. In
multithreaded runs, rows from one event remain together, but event blocks can
appear in completion order rather than numeric event order.

Interaction transport follows the detector model. The legacy `nested_cell`
model preserves its elastic-only workflow by aborting the current event after
the first neutron-inelastic interaction (the triggering interaction itself is
still recorded). The general `box_cryostat` model allows the complete
inelastic cascade and secondary-neutron tracks to continue.

Every event also receives a particle-resolved energy, scintillation, and PE summary
in `particle_class_summary_<config_id>.csv`, without the detected-neutron,
positive-PE, or TOF cuts used for the legacy event CSV. Energy deposited in
active LAr, scintillation photons produced there, and detected PE are split into primary
neutron, secondary neutron, gamma, electron/positron, proton, alpha, nuclear
recoil, and other classes. The total columns are the sums of those classes.
In multithreaded runs, rows can appear in completion order rather than event-ID
order.

Scintillation photons are attributed to their parent particle after the
configured ion-yield thinning. That origin is propagated from the original VUV
photon to its wavelength-shifted descendants, allowing sensor detections to be
assigned to the same particle class. This preserves the meaning of the legacy
`numPhotons` field. For `nested_cell`, a rejected inelastic event's summary contains only
energy and light produced before its intentional abort; `box_cryostat`
summaries cover the complete cascade.

### Neutron-run output reference

| Output | Scope and purpose |
| --- | --- |
| `numPE_<top>_topPMTs_<bottom>_botPMTs_<rows>_SiPMRows.csv` | Legacy selected-event stream: external neutron detected, positive PE, and TOF inside the configured window. Contains PE, produced photons, active-LAr deposit, detector label, and legacy interaction counters. |
| `dune_event_response_<config_id>.csv` | One uncut box-model row per event with total, electronic-recoil, nuclear-recoil, and other PE, produced photons, deposited energy, and interaction counts. Requires neither a tagger nor TOF. |
| `neutron_interactions_<config_id>.csv` | Uncut interaction-level truth for primary and secondary neutron tracks. Events with no recordable neutron interaction naturally have no row. |
| `particle_class_summary_<config_id>.csv` | One uncut row per event with total and per-particle-class active-LAr energy, produced light, and detected PE. |
| `scan_metadata.json` | Exact macro/config values for every scan row; also maps a selected-event filename back to its interaction-truth `config_id`. |
| `config_used.yaml` | Archived YAML used by the run. |
| `sensor_scan_used.csv` | Archived sensor scan used by the run. |

The event ID is the join key within one configuration. Multithreaded output can
be in completion order, so join or sort by event ID rather than assuming row
order. Across several scan configurations, use both `config_id` (or source
filename) and event ID because Geant4 event IDs restart for every macro.

## Neutron scintillation correction

For newly simulated box-model data, gamma calibration normalizes only the
nuclear-recoil optical component:

```text
N_ph,expected,gamma = a * E_dep,NR,MeV + b
NR_correction       = N_ph,expected,gamma / N_ph,NR,simulated
numPE_NR,corrected  = numPE_NR,raw * NR_correction
numPE_total,corr    = numPE_ER,raw + numPE_other,raw + numPE_NR,corrected
```

Electronic-recoil and other-origin PE remain untouched. This is evaluated
separately for every event using the recoil-resolved columns in
`dune_event_response_<config_id>.csv`. The formula does not multiply by Leff
because the configured ion efficiency was already sampled during photon
production. Events without NR light retain their raw total PE.

The nuclear-recoil group matches the particle classes using the configured ion
scintillation-yield scale: primary/secondary neutrons, alpha particles, and
nuclear recoil ions. In practice, scintillation is created by charged tracks.
Gamma and electron/positron origins form the electronic-recoil group. Proton
and unidentified origins are retained explicitly in `other` to preserve the
historical proton-yield policy and make the component sum auditable.

Legacy selected-event CSVs do not contain recoil-split PE, so their correction
continues to normalize the whole event as before.

### Legacy post-analysis correction

The previous deterministic correction is retained only for reproducing older
runs:

```text
E_dep               = active-LAr energy deposit from the event CSV
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

For current runs, each selected `eventID` is joined to its matching
`neutron_interactions_<config_id>.csv`; this includes interactions by primary
and secondary neutron tracks. Categories are mutually exclusive:

- **No recorded interaction** and **transport exit only**;
- **single elastic** and **multiple elastic**;
- **inelastic** and **elastic + inelastic**;
- **capture**, **fission**, **other hadronic**, and **unclassified**.

Fission, capture, other-hadronic, and unclassified records take precedence over
the elastic/inelastic categories; the event-level output retains every channel
count. A true-single-elastic event used by the PE money plots must contain
exactly one elastic interaction in fiducial LAr and no other hadronic
interaction anywhere. When interaction-truth files are absent, older run
directories remain analyzable through a documented approximation using the
legacy event columns.

Outputs under `output/event_type_analysis/` include:

- `event_type_frequencies.csv` with file-, run-, and detector-level counts;
- `passed_event_topologies.csv` with each selected event's topology and, when
  interaction truth is available, its channel counts, fiducial-elastic count,
  and primary/secondary-neutron counts;
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

## Helper-based nested-cell example

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
