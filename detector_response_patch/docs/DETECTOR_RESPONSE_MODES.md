# Detector response modes

The patched files add three optical-detector response modes for PMTs and SiPMs.
They are selected in `config/example_config.yaml` and inserted into the generated
Geant4 macro.

## Modes

```yaml
sensors:
  pmt:
    response_mode: builtin_curve   # flat, builtin_curve, csv_curve
    pde: 0.266
    top_pde: 0.266
    bottom_pde: 0.266
    reflectivity: 0.10
    qe_curve_file: none
  sipm:
    response_mode: builtin_curve
    pde: 0.35
    reflectivity: 0.10
    qe_curve_file: none
```

### `flat`
Uses a wavelength-independent detector response:

```text
physical detection probability = pde
Geant4 surface EFFICIENCY = pde / (1 - reflectivity)
```

### `builtin_curve`
Uses the hard-coded PMT/SiPM quantum-efficiency curves already in
`DetectorConstruction.cc`. The `pde` value scales the peak of the curve.
For example, if the built-in SiPM curve peaks at 0.402 and `pde: 0.35`,
the entire curve is multiplied by `0.35 / 0.402` before the usual Geant4
surface correction.

Set `pde: 0` to use the built-in curve exactly as written, with no peak scaling.

### `csv_curve`
Reads a user-provided CSV file. The CSV may be either:

```csv
wavelength_nm,qe
350,0.12
400,0.31
425,0.36
```

or:

```csv
energy_eV,qe
2.5,0.20
2.8,0.32
3.1,0.36
```

`qe` may be written as a fraction or a percent. For example, both `0.35` and
`35` are accepted.

As in `builtin_curve`, `pde` scales the curve peak. Set `pde: 0` to use the
CSV curve as an absolute curve.

## Macro commands added

```text
/det/setPMTResponseMode flat|builtin_curve|csv_curve
/det/setSiPMResponseMode flat|builtin_curve|csv_curve
/det/setPMTCurveFile path/to/pmt.csv
/det/setSiPMCurveFile path/to/sipm.csv
```

## Important convention

The config value `pde` should always be the intended physical detection
probability, not the internally corrected Geant4 `EFFICIENCY` value. The C++
code performs the correction internally:

```text
EFFICIENCY = PDE / (1 - reflectivity)
```


## Raw built-in/CSV curves

Use `builtin_curve_raw` to use the built-in wavelength-dependent PMT/SiPM curve exactly as it is stored in the C++ code. In this mode the configured `pde`, `top_pde`, and `bottom_pde` values are ignored for that sensor type.

Use `csv_curve_raw` to use a user-provided CSV curve as an absolute QE/PDE curve, without scaling its peak to the configured PDE.

Use `builtin_curve` or `csv_curve` when you want the curve shape to be preserved but the peak to be scaled to the configured PDE.
