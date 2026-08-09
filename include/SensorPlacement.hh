#ifndef SensorPlacement_h
#define SensorPlacement_h 1

#include "G4RotationMatrix.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;

enum class SensorType {
  PMT,
  SiPM,
};

enum class SensorSurface {
  CylinderTop,
  CylinderBottom,
  CylinderBarrel,
  BoxPositiveX,
  BoxNegativeX,
  BoxPositiveY,
  BoxNegativeY,
  BoxPositiveZ,
  BoxNegativeZ,
  Explicit,
};

// Shape-independent description of one photosensor transform. Layout
// generators create these records; detector construction consumes them.
struct SensorPlacement {
  G4String id;
  G4String group;
  G4String physicalVolumeName;
  SensorType type = SensorType::SiPM;
  SensorSurface surface = SensorSurface::Explicit;
  G4LogicalVolume* hostLogicalVolume = nullptr;
  G4VPhysicalVolume* hostPhysicalVolume = nullptr;
  G4ThreeVector position;
  G4RotationMatrix rotation;
  G4bool hasRotation = false;
  G4double activeWidth = 0.0;
  G4double activeHeight = 0.0;
  G4double thickness = 0.0;
  G4int groupIndex = 0;
  G4int copyNumber = 0;
};

struct BoxFaceGridConfig {
  G4String group;
  G4String physicalVolumeName = "SensorTilePV";
  SensorType type = SensorType::SiPM;
  SensorSurface face = SensorSurface::BoxPositiveX;
  G4ThreeVector boxDimensions;
  G4double activeWidth = 0.0;
  G4double activeHeight = 0.0;
  G4double thickness = 0.0;
  G4double pitchU = 0.0;
  G4double pitchV = 0.0;
  G4double edgeClearance = 0.0;
  G4double inset = 0.0;
  G4LogicalVolume* hostLogicalVolume = nullptr;
  G4VPhysicalVolume* hostPhysicalVolume = nullptr;
};

std::vector<SensorPlacement>
GenerateBoxFaceGridPlacements(const BoxFaceGridConfig& config);

// Six non-overlapping TPB slabs forming an inner coating for a box.  Each
// wallInset is the distance from the corresponding box wall to the slab face
// nearest that wall.  Different values are supported so the coating can sit
// in front of configurable sensor depths on each face.
struct BoxTPBSlabConfig {
  G4ThreeVector boxDimensions;
  std::array<G4double, 6> wallInsets{};
  G4double thickness = 0.0;
};

struct BoxTPBSlabPlacement {
  SensorSurface face = SensorSurface::BoxPositiveX;
  G4ThreeVector position;
  G4ThreeVector dimensions;
};

std::vector<BoxTPBSlabPlacement>
GenerateBoxTPBSlabPlacements(const BoxTPBSlabConfig& config);

#endif
