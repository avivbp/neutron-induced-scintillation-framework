#ifndef SensorPlacement_h
#define SensorPlacement_h 1

#include "G4RotationMatrix.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

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

#endif
