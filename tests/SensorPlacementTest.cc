#include "SensorPlacement.hh"

#include <array>
#include <cmath>
#include <stdexcept>

namespace {

constexpr G4double kTolerance = 1e-9;

G4bool Near(G4double left, G4double right)
{
  return std::abs(left - right) < kTolerance;
}

void Require(G4bool condition, const char* message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void CheckNormal(const SensorPlacement& placement,
                 const G4ThreeVector& expected)
{
  const auto normal = placement.rotation * G4ThreeVector(0, 0, 1);
  Require(Near(normal.x(), expected.x()), "incorrect sensor normal x");
  Require(Near(normal.y(), expected.y()), "incorrect sensor normal y");
  Require(Near(normal.z(), expected.z()), "incorrect sensor normal z");
}

BoxFaceGridConfig BaseConfig(SensorSurface face)
{
  BoxFaceGridConfig config;
  config.group = "wall_sipms";
  config.type = SensorType::SiPM;
  config.face = face;
  config.boxDimensions = {100.0, 60.0, 40.0};
  config.activeWidth = 10.0;
  config.activeHeight = 10.0;
  config.thickness = 2.0;
  config.pitchU = 20.0;
  config.pitchV = 15.0;
  config.edgeClearance = 5.0;
  config.inset = 1.0;
  return config;
}

}  // namespace

int main()
{
  const std::array<std::pair<SensorSurface, G4ThreeVector>, 6> faces = {{
      {SensorSurface::BoxPositiveX, {-1, 0, 0}},
      {SensorSurface::BoxNegativeX, {1, 0, 0}},
      {SensorSurface::BoxPositiveY, {0, -1, 0}},
      {SensorSurface::BoxNegativeY, {0, 1, 0}},
      {SensorSurface::BoxPositiveZ, {0, 0, -1}},
      {SensorSurface::BoxNegativeZ, {0, 0, 1}},
  }};

  for (const auto& [face, normal] : faces) {
    const auto placements = GenerateBoxFaceGridPlacements(BaseConfig(face));
    Require(!placements.empty(), "box face unexpectedly has no placements");
    for (const auto& placement : placements) {
      CheckNormal(placement, normal);
      Require(placement.id.find("wall_sipms/") == 0,
              "sensor ID does not start with group");
      Require(placement.hasRotation, "box sensor has no rotation");
    }
  }

  const auto positiveX = GenerateBoxFaceGridPlacements(
      BaseConfig(SensorSurface::BoxPositiveX));
  Require(positiveX.size() == 6, "unexpected +x grid size");
  Require(Near(positiveX.front().position.x(), 48.0), "incorrect face inset");
  Require(Near(positiveX.front().position.y(), -20.0), "incorrect first u");
  Require(Near(positiveX.front().position.z(), 7.5), "incorrect first v");
  Require(Near(positiveX.back().position.y(), 20.0), "incorrect last u");
  Require(Near(positiveX.back().position.z(), -7.5), "incorrect last v");

  auto tooLarge = BaseConfig(SensorSurface::BoxPositiveX);
  tooLarge.activeWidth = 100.0;
  tooLarge.pitchU = 100.0;
  Require(GenerateBoxFaceGridPlacements(tooLarge).empty(),
          "oversize sensors should produce no placements");

  auto overlapping = BaseConfig(SensorSurface::BoxPositiveX);
  overlapping.pitchU = 5.0;
  G4bool rejected = false;
  try {
    GenerateBoxFaceGridPlacements(overlapping);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  Require(rejected, "overlapping sensor pitch was not rejected");

  return 0;
}
