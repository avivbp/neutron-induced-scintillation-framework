# General Detector and DUNE Extension Design

## Goal

Extend the framework from its nested cylindrical LAr-cell geometry to detectors
with one or more active LAr volumes, including a DUNE-like box-shaped cryostat,
without breaking the existing small-cell workflow.

The extension has two independent axes:

1. physics analysis must use configured volume roles instead of physical-volume
   names such as `innerCell`; and
2. optical sensors must be placed on detector surfaces using shape-aware layouts
   instead of assuming top/bottom PMT groups and cylindrical SiPM rings.

## Compatibility baseline

The existing nested-cell configuration remains supported and is the regression
baseline. Before a new detector shape is introduced, the volume-role and sensor
layout abstractions must reproduce the current geometry, event columns, photon
counts, and sensor activation behavior for the existing example configuration.

Legacy output column names may be retained temporarily, but new code must not
assign physics meaning based solely on `innerCell` or `outerCell` names.

## Detector volume roles

Logical volumes can have one or more analysis roles:

- `active_lar`: energy deposition and scintillation are counted;
- `fiducial_lar`: interactions are eligible for fiducial selections;
- `inactive_lar`: LAr present for transport but excluded from the active signal;
- `cryostat`: detector structure surrounding the LAr;
- `optical_detector`: a photon-detection surface; and
- `external_detector`: an optional neutron or auxiliary detector.

Detector construction owns the role registry and exposes queries using logical
volume identity. Stepping and event actions query these roles rather than
comparing strings. A volume may be both `active_lar` and `fiducial_lar`.

Initial mappings are:

| Detector model | Logical volume | Roles |
| --- | --- | --- |
| Nested cell | `innerCell` | `active_lar`, `fiducial_lar` |
| Nested cell, optional | `outerCell` | `active_lar` or `inactive_lar` |
| DUNE-like cryostat | configured LAr volume(s) | `active_lar` |
| DUNE-like cryostat | configured inset volume(s) | `fiducial_lar` |

## Interaction truth

Event summaries will count all neutron channels, not only clean elastic
scatters. The initial channel vocabulary is:

- `elastic` (`hadElastic`);
- `inelastic` (`neutronInelastic`);
- `capture` (`nCapture`);
- `fission` when provided by the selected physics list;
- `other_hadronic`;
- `transport_exit`; and
- `unclassified`, retaining the original Geant4 process name.

An interaction record should retain event, track and parent identifiers,
position, time, volume role, process and channel, neutron energy before and
after the interaction, local energy deposit, and produced secondaries. Event
summaries are derived from these records rather than replacing them.

The implementation writes these records to
`neutron_interactions_<config_id>.csv` without applying the optical-detection
or TOF cuts used for the legacy event CSV. Secondary lists are encoded as
`particle:PDG-code:kinetic-energy-keV` entries. Ordinary geometry crossings
are omitted, while transportation out of the simulation world is retained as
`transport_exit`.

The detector exposes this behavior through the generic
`AllowsAllNeutronInteractions` capability. It is enabled for `box_cryostat`;
the legacy `nested_cell` model retains its elastic-only selection by aborting
an event immediately after recording its first neutron-inelastic interaction.

Energy deposition and scintillation must also be attributable to depositing
particle classes. Inelastic events can mix nuclear recoils, secondary neutrons,
gammas, electrons, protons, alphas, and heavier ions; total inelastic deposited
energy must not be treated as one nuclear recoil.

The implementation writes one uncut event row to
`particle_class_summary_<config_id>.csv`. Active-LAr energy and produced
scintillation photons are split into primary neutron, secondary neutron, gamma,
electron/positron, proton, alpha, nuclear recoil, and other classes. Totals are
computed from the class values. Photon attribution occurs after the existing
ion-yield thinning and before optical transport, wavelength shifting, and
sensor response. Consequently, an intentionally aborted `nested_cell`
inelastic event contains only pre-abort deposits, while `box_cryostat` records
the complete inelastic cascade.

## Shape-aware optical sensor placement

Sensor response (PDE, reflectivity, wavelength dependence, and active/inactive
state) remains separate from sensor placement.

A sensor layout produces a list of placements, each containing:

- sensor type and stable identifier;
- host active-LAr volume;
- position and orientation;
- face or surface identifier;
- active dimensions; and
- grouping metadata used by coverage scans.

Placement strategies are shape-specific:

### Cylinder

- planar grids on the top and bottom caps;
- azimuthal rings or staggered tiles on the barrel; and
- configurable edge clearance and tile pitch.

The first cylinder strategy must reproduce the existing four-top, four-bottom,
and SiPM-row convention.

### Box

- independent two-dimensional grids on the `+x`, `-x`, `+y`, `-y`, `+z`, and
  `-z` faces;
- face selection and per-face sensor type;
- configurable pitch, edge clearance, offsets, and optional staggering; and
- correct inward-facing orientation for every face.

Tiles are placed only when their full active footprint fits inside the selected
face after clearance is applied. Coverage is computed from placed active area
divided by the selected instrumentable surface area, rather than from a
cylinder-specific row formula.

Future imported or irregular geometries can supply explicit sensor transforms
without changing photon detection or analysis code.

## Configuration direction

The legacy sensor-scan columns remain supported during migration. New layouts
will use named placement groups, for example:

```yaml
detector:
  model: box_cryostat

geometry:
  active_lar:
    shape: box
    dimensions_cm: [1400.0, 1200.0, 6000.0]

sensors:
  layouts:
    - id: long_wall_sipms
      type: sipm
      faces: [+x, -x]
      pattern: grid
      pitch_cm: [20.0, 20.0]
      edge_clearance_cm: 10.0
    - id: endcap_sipms
      type: sipm
      faces: [+z, -z]
      pattern: grid
      pitch_cm: [20.0, 20.0]
      edge_clearance_cm: 10.0
```

Exact DUNE dimensions and photon-detector layouts are inputs to this model, not
hard-coded defaults.

## Implementation checkpoints

Each checkpoint is implemented, tested, committed, and pushed independently:

1. add the logical-volume role registry while preserving nested-cell results;
2. replace active/fiducial `innerCell` checks with role queries;
3. introduce shape-independent sensor placement records;
4. reproduce the current cylindrical sensor layout through the new interface;
5. add box-face grid placement and geometry/coverage tests;
6. add the configurable box-cryostat detector model;
7. write interaction-level truth for every neutron channel;
8. add particle-class energy and scintillation bookkeeping; and
9. update analyses with channel-aware topology categories.

Checkpoint 9 joins cut-selected neutron event rows to the matching per-config
interaction-truth file. The analysis distinguishes no-record/transport-only,
single/multiple elastic, inelastic, elastic-plus-inelastic, capture, fission,
other-hadronic, and unclassified topologies. It writes both frequency summaries
and event-level channel counts. The PE analysis uses truth-defined single
fiducial elastics when available and falls back to legacy event columns for
pre-truth run directories.

Large DUNE-scale optical transport may require photon libraries, fast optical
simulation, or parameterized response. That performance choice is deliberately
kept outside the geometry, interaction-truth, and placement interfaces.
