# Configurable detector-response patch

Copy the files into your `Code/LabSimulation` directory:

```bash
cp include/DetectorConstruction.hh include/
cp include/DetectorMessenger.hh include/
cp src/DetectorConstruction.cc src/
cp src/DetectorMessenger.cc src/
cp scripts/config_tools.py scripts/
cp scripts/run_scan.py scripts/
cp macros/template_run.mac macros/
cp config/example_config.yaml config/
```

Then rebuild and regenerate macros:

```bash
cmake --build build-mt -j$(nproc)
python3 scripts/run_scan.py --config config/example_config.yaml --scan config/sensor_scan.csv --dry-run
cat macros/generated/A0.mac
python3 scripts/run_scan.py --config config/example_config.yaml --scan config/sensor_scan.csv
```

The generated macro now includes pre-initialization geometry/material commands
and idle-state detector-response commands. See `docs/DETECTOR_RESPONSE_MODES.md`.


## Raw built-in/CSV curves

Use `builtin_curve_raw` to use the built-in wavelength-dependent PMT/SiPM curve exactly as it is stored in the C++ code. In this mode the configured `pde`, `top_pde`, and `bottom_pde` values are ignored for that sensor type.

Use `csv_curve_raw` to use a user-provided CSV curve as an absolute QE/PDE curve, without scaling its peak to the configured PDE.

Use `builtin_curve` or `csv_curve` when you want the curve shape to be preserved but the peak to be scaled to the configured PDE.
