#include "SensorPlacement.hh"

#include <array>
#include <cmath>
#include <stdexcept>

namespace {

struct FaceFrame {
  G4String label;
  G4ThreeVector u;
  G4ThreeVector v;
  G4ThreeVector inward;
  G4ThreeVector center;
  G4double extentU;
  G4double extentV;
};

G4bool IsBoxFace(SensorSurface face)
{
  return face >= SensorSurface::BoxPositiveX &&
         face <= SensorSurface::BoxNegativeZ;
}

std::size_t BoxFaceIndex(SensorSurface face)
{
  if (!IsBoxFace(face)) {
    throw std::invalid_argument("surface is not a box face");
  }
  return static_cast<std::size_t>(face) -
         static_cast<std::size_t>(SensorSurface::BoxPositiveX);
}

FaceFrame MakeFaceFrame(SensorSurface face, const G4ThreeVector& dimensions)
{
  const G4double hx = 0.5 * dimensions.x();
  const G4double hy = 0.5 * dimensions.y();
  const G4double hz = 0.5 * dimensions.z();

  switch (face) {
    case SensorSurface::BoxPositiveX:
      return {"+x", {0, 1, 0}, {0, 0, -1}, {-1, 0, 0}, {hx, 0, 0},
              dimensions.y(), dimensions.z()};
    case SensorSurface::BoxNegativeX:
      return {"-x", {0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {-hx, 0, 0},
              dimensions.y(), dimensions.z()};
    case SensorSurface::BoxPositiveY:
      return {"+y", {1, 0, 0}, {0, 0, 1}, {0, -1, 0}, {0, hy, 0},
              dimensions.x(), dimensions.z()};
    case SensorSurface::BoxNegativeY:
      return {"-y", {1, 0, 0}, {0, 0, -1}, {0, 1, 0}, {0, -hy, 0},
              dimensions.x(), dimensions.z()};
    case SensorSurface::BoxPositiveZ:
      return {"+z", {1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {0, 0, hz},
              dimensions.x(), dimensions.y()};
    case SensorSurface::BoxNegativeZ:
      return {"-z", {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -hz},
              dimensions.x(), dimensions.y()};
    default:
      throw std::invalid_argument("sensor grid face is not a box face");
  }
}

G4int GridCount(G4double extent, G4double activeSize, G4double pitch,
                G4double clearance)
{
  const G4double usable = extent - 2.0 * clearance;
  if (activeSize > usable) {
    return 0;
  }
  return static_cast<G4int>(std::floor((usable - activeSize) / pitch)) + 1;
}

}  // namespace

std::vector<SensorPlacement>
GenerateBoxFaceGridPlacements(const BoxFaceGridConfig& config)
{
  if (!IsBoxFace(config.face)) {
    throw std::invalid_argument("box sensor grid requires one of the six box faces");
  }
  if (config.group.empty()) {
    throw std::invalid_argument("box sensor grid group must not be empty");
  }
  if (config.boxDimensions.x() <= 0.0 || config.boxDimensions.y() <= 0.0 ||
      config.boxDimensions.z() <= 0.0) {
    throw std::invalid_argument("box dimensions must be positive");
  }
  if (config.activeWidth <= 0.0 || config.activeHeight <= 0.0 ||
      config.thickness <= 0.0) {
    throw std::invalid_argument("sensor dimensions must be positive");
  }
  if (config.pitchU < config.activeWidth ||
      config.pitchV < config.activeHeight) {
    throw std::invalid_argument("sensor pitch must be at least the active size");
  }
  if (config.edgeClearance < 0.0 || config.inset < 0.0) {
    throw std::invalid_argument("sensor clearance and inset must be non-negative");
  }

  const auto frame = MakeFaceFrame(config.face, config.boxDimensions);
  const G4int countU = GridCount(frame.extentU, config.activeWidth,
                                config.pitchU, config.edgeClearance);
  const G4int countV = GridCount(frame.extentV, config.activeHeight,
                                config.pitchV, config.edgeClearance);
  if (countU == 0 || countV == 0) {
    return {};
  }

  const G4double startU = -0.5 * (countU - 1) * config.pitchU;
  const G4double startV = -0.5 * (countV - 1) * config.pitchV;
  const G4ThreeVector surfaceOffset =
      frame.inward * (config.inset + 0.5 * config.thickness);
  const G4RotationMatrix rotation(frame.u, frame.v, frame.inward);

  std::vector<SensorPlacement> placements;
  placements.reserve(countU * countV);
  for (G4int row = 0; row < countV; ++row) {
    for (G4int column = 0; column < countU; ++column) {
      const G4int index = row * countU + column;
      SensorPlacement placement;
      placement.id = config.group + "/" + frame.label + "/" +
                     std::to_string(row) + "/" + std::to_string(column);
      placement.group = config.group;
      placement.physicalVolumeName = config.physicalVolumeName;
      placement.type = config.type;
      placement.surface = config.face;
      placement.hostLogicalVolume = config.hostLogicalVolume;
      placement.hostPhysicalVolume = config.hostPhysicalVolume;
      placement.position = frame.center + surfaceOffset +
          frame.u * (startU + column * config.pitchU) +
          frame.v * (startV + row * config.pitchV);
      placement.rotation = rotation;
      placement.hasRotation = true;
      placement.activeWidth = config.activeWidth;
      placement.activeHeight = config.activeHeight;
      placement.thickness = config.thickness;
      placement.groupIndex = row;
      placement.copyNumber = index;
      placements.push_back(placement);
    }
  }

  return placements;
}

std::vector<BoxTPBSlabPlacement>
GenerateBoxTPBSlabPlacements(const BoxTPBSlabConfig& config)
{
  if (config.boxDimensions.x() <= 0.0 ||
      config.boxDimensions.y() <= 0.0 ||
      config.boxDimensions.z() <= 0.0) {
    throw std::invalid_argument("box dimensions must be positive");
  }
  if (config.thickness <= 0.0) {
    throw std::invalid_argument("TPB thickness must be positive");
  }
  for (const auto inset : config.wallInsets) {
    if (inset < 0.0) {
      throw std::invalid_argument("TPB wall inset must be non-negative");
    }
  }

  const auto depth = [&](SensorSurface face) {
    return config.wallInsets[BoxFaceIndex(face)] + config.thickness;
  };
  const auto centerBetween = [](G4double low, G4double high) {
    return 0.5 * (low + high);
  };

  const G4double hx = 0.5 * config.boxDimensions.x();
  const G4double hy = 0.5 * config.boxDimensions.y();
  const G4double hz = 0.5 * config.boxDimensions.z();
  const G4double t = config.thickness;

  const G4double xLow = -hx + depth(SensorSurface::BoxNegativeX);
  const G4double xHigh = hx - depth(SensorSurface::BoxPositiveX);
  const G4double yLow = -hy + depth(SensorSurface::BoxNegativeY);
  const G4double yHigh = hy - depth(SensorSurface::BoxPositiveY);
  const G4double zLow = -hz + depth(SensorSurface::BoxNegativeZ);
  const G4double zHigh = hz - depth(SensorSurface::BoxPositiveZ);
  if (xHigh <= xLow || yHigh <= yLow || zHigh <= zLow) {
    throw std::invalid_argument("TPB slab depths leave no box interior");
  }

  // Partition the inner coating without overlaps: Z slabs span complete end
  // faces, Y slabs stop at the Z slabs, and X slabs stop at both Y and Z.
  // Adjacent pieces therefore meet at their edges while every ray normal to a
  // coated box face crosses TPB before reaching a wall-mounted sensor.
  return {
      {SensorSurface::BoxPositiveX,
       {hx - config.wallInsets[BoxFaceIndex(SensorSurface::BoxPositiveX)] -
            0.5 * t,
        centerBetween(yLow, yHigh), centerBetween(zLow, zHigh)},
       {t, yHigh - yLow, zHigh - zLow}},
      {SensorSurface::BoxNegativeX,
       {-hx + config.wallInsets[BoxFaceIndex(SensorSurface::BoxNegativeX)] +
             0.5 * t,
        centerBetween(yLow, yHigh), centerBetween(zLow, zHigh)},
       {t, yHigh - yLow, zHigh - zLow}},
      {SensorSurface::BoxPositiveY,
       {0.0,
        hy - config.wallInsets[BoxFaceIndex(SensorSurface::BoxPositiveY)] -
            0.5 * t,
        centerBetween(zLow, zHigh)},
       {config.boxDimensions.x(), t, zHigh - zLow}},
      {SensorSurface::BoxNegativeY,
       {0.0,
        -hy + config.wallInsets[BoxFaceIndex(SensorSurface::BoxNegativeY)] +
              0.5 * t,
        centerBetween(zLow, zHigh)},
       {config.boxDimensions.x(), t, zHigh - zLow}},
      {SensorSurface::BoxPositiveZ,
       {0.0, 0.0,
        hz - config.wallInsets[BoxFaceIndex(SensorSurface::BoxPositiveZ)] -
            0.5 * t},
       {config.boxDimensions.x(), config.boxDimensions.y(), t}},
      {SensorSurface::BoxNegativeZ,
       {0.0, 0.0,
        -hz + config.wallInsets[BoxFaceIndex(SensorSurface::BoxNegativeZ)] +
              0.5 * t},
       {config.boxDimensions.x(), config.boxDimensions.y(), t}},
  };
}
