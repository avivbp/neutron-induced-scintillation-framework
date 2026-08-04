//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "EventAction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4Tubs.hh"
#include "G4OpticalSurface.hh"
#include "G4Orb.hh"
#include "G4Colour.hh"
#include "G4Scintillation.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <unordered_map>
#include <vector>

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  enum class VolumeRole : std::uint32_t {
    ActiveLAr       = 1u << 0,
    FiducialLAr     = 1u << 1,
    InactiveLAr     = 1u << 2,
    Cryostat        = 1u << 3,
    OpticalDetector = 1u << 4,
    ExternalDetector = 1u << 5,
  };

private:
  std::unordered_map<const G4LogicalVolume*, std::uint32_t> fVolumeRoles;

  // ------------------ Config knobs (set by messenger) ------------------
  G4int    fTopPMTs   = 0;     // 0..4 active tiles
  G4int    fBotPMTs   = 0;     // 0..4 active tiles
  G4int    fSiPMRows  = 0;     // 0..16 active rows

  G4double fRefl_PMT  = 0.35;  // reflectivity of active PMT tile
  G4double fRefl_SiPM = 0.4;  // reflectivity of active SiPM tile

  // User-facing PDE values are the intended physical detection probabilities.
  // The conversion to Geant4 optical-surface EFFICIENCY = PDE/(1-reflectivity)
  // is done only when building the material properties table.
  G4double fPDE_TopPMT = 0.266; // requested peak/flat PDE
  G4double fPDE_BotPMT = 0.266; // requested peak/flat PDE
  G4double fPDE_SiPM  = 0.35;  // requested peak/flat PDE

  enum class SensorResponseMode {
    Flat,
    BuiltinCurve,     // built-in spectral shape, scaled to requested peak PDE
    BuiltinCurveRaw,  // built-in curve used exactly as coded
    CsvCurve,         // CSV spectral shape, scaled to requested peak PDE
    CsvCurveRaw       // CSV values used exactly as supplied
  };

  SensorResponseMode fPMTResponseMode  = SensorResponseMode::BuiltinCurve;
  SensorResponseMode fSiPMResponseMode = SensorResponseMode::BuiltinCurve;
  G4String fPMTCurveFile  = "";
  G4String fSiPMCurveFile = "";

  G4MaterialPropertiesTable* fMPT_ESR = nullptr;
  G4MaterialPropertiesTable* fMPT_TopPMT_Active = nullptr;
  G4MaterialPropertiesTable* fMPT_BotPMT_Active = nullptr;
  G4MaterialPropertiesTable* fMPT_SiPM_Active = nullptr;

  G4double fCached_PDE_TopPMT  = -1.0;
  G4double fCached_PDE_BotPMT  = -1.0;
  G4double fCached_Refl_PMT = -1.0;
  G4double fCached_PDE_SiPM = -1.0;
  G4double fCached_Refl_SiPM = -1.0;


  // ------------------ Patch bookkeeping ------------------
  std::vector<G4VPhysicalVolume*> fTopPMT_PV;
  std::vector<G4VPhysicalVolume*> fBotPMT_PV;
  std::vector<G4OpticalSurface*>  fTopPMT_Surf;
  std::vector<G4OpticalSurface*>  fBotPMT_Surf;

  std::vector<G4VPhysicalVolume*> fSiPM_PV;
  std::vector<G4OpticalSurface*>  fSiPM_Surf;
  std::vector<int>                fSiPM_RowOfTile;

  G4int fSiPM_TilesPerRow = 0; // computed

  // ------------------ Geometry config ------------------
  G4double fInnerRadius = 5.0 * CLHEP::cm;
  G4double fInnerHeight = 10.0 * CLHEP::cm;
  G4double fOuterDiameter = 35.0 * CLHEP::cm;
  G4double fOuterHeight = 70.0 * CLHEP::cm;

  struct NeutronDetectorConfig {
    G4String name;
    G4double angle;
    G4double distance;
  };

  std::vector<NeutronDetectorConfig> fNeutronDetectors = {
    {"A0", 25.0 * CLHEP::deg, 100.0 * CLHEP::cm},
    {"A1", 40.0 * CLHEP::deg, 100.0 * CLHEP::cm},
    {"A2", 50.0 * CLHEP::deg, 100.0 * CLHEP::cm},
    {"A3", 60.0 * CLHEP::deg, 100.0 * CLHEP::cm},
    {"A4", 90.0 * CLHEP::deg, 100.0 * CLHEP::cm},
  };

  // ------------------ Sensor geometry config ------------------
  G4double fPMTSide = 2.54 * CLHEP::cm;
  G4double fSiPMTile = 0.6 * CLHEP::cm;

  // ------------------ Event-selection config ------------------
  G4double fTOFMinNs = 40.0;
  G4double fTOFMaxNs = 50.0;
  G4double fPrimaryEnergy = 2.5 * CLHEP::MeV;

  // ------------------ Optical/scintillation config ------------------
  G4double fLArScintYieldPerMeV = 51300.0;
  G4double fLArScintYieldScale = 1.0;
  // Nuclear-recoil scintillation efficiency applied to LAr scintillation
  // photons after Geant4's Birks calculation and before optical transport.
  G4double fLArIonScintYieldScale = 0.3;
  G4double fLArAbsLength = 150.0 * CLHEP::cm;
  G4double fLArFastTime = 7.0 * CLHEP::ns;
  G4double fLArSlowTime = 1500.0 * CLHEP::ns;
  G4double fLArFastFraction = 1.0;
  G4double fTPBEfficiency = 1.0;
  G4double fReflectorReflectivity = 0.98;

  // ------------------ Build + apply methods ------------------

public:
  void RegisterVolumeRole(const G4LogicalVolume* volume, VolumeRole role);
  G4bool HasVolumeRole(const G4LogicalVolume* volume, VolumeRole role) const;
  void ClearVolumeRoles();

  void BuildPMTPatches();   // builds 4 top + 4 bottom (max)
  void BuildSiPMPatches();  // builds max rows and tiles
  void ApplySensorConfig(); // updates EFFICIENCY/REFLECTIVITY per patch
  void UpdateConfigurableOpticalProperties();

  // ------------------ Setters used by messenger ------------------
  public:
  void SetTopPMTs(G4int v)   { fTopPMTs  = v; }
  void SetBotPMTs(G4int v)   { fBotPMTs  = v; }
  void SetSiPMRows(G4int v)  { fSiPMRows = v; }

  G4int GetTopPMTs()  const  {return fTopPMTs;}
  G4int GetBotPMTs()  const  {return fBotPMTs;}
  G4int GetSiPMRows() const  {return fSiPMRows;}


  void SetPDETopPMT(G4double v) { fPDE_TopPMT = v; fMPT_TopPMT_Active = nullptr; }
  void SetPDEBotPMT(G4double v) { fPDE_BotPMT = v; fMPT_BotPMT_Active = nullptr; }
  void SetPDESiPM(G4double v)  { fPDE_SiPM  = v; fMPT_SiPM_Active = nullptr; }

  //void SetReflESR(G4double v)  { fRefl_ESR  = v; }
  void SetReflPMT(G4double v)  { fRefl_PMT  = v; fMPT_TopPMT_Active = nullptr; fMPT_BotPMT_Active = nullptr; }
  void SetReflSiPM(G4double v) { fRefl_SiPM = v; fMPT_SiPM_Active = nullptr; }

  void SetPMTResponseMode(const G4String& mode) {
      const auto m = ToLower(mode);
      if (m == "flat") {
        fPMTResponseMode = SensorResponseMode::Flat;
      } else if (m == "builtin_curve_raw" || m == "builtin_raw" || m == "raw_builtin") {
        fPMTResponseMode = SensorResponseMode::BuiltinCurveRaw;
      } else if (m == "csv_curve_raw" || m == "csv_raw" || m == "raw_csv") {
        fPMTResponseMode = SensorResponseMode::CsvCurveRaw;
      } else if (m == "csv_curve" || m == "csv" || m == "file") {
        fPMTResponseMode = SensorResponseMode::CsvCurve;
      } else {
        fPMTResponseMode = SensorResponseMode::BuiltinCurve;
      }
      fMPT_TopPMT_Active = nullptr;
      fMPT_BotPMT_Active = nullptr;
  }

  void SetSiPMResponseMode(const G4String& mode) {
      const auto m = ToLower(mode);
      if (m == "flat") {
        fSiPMResponseMode = SensorResponseMode::Flat;
      } else if (m == "builtin_curve_raw" || m == "builtin_raw" || m == "raw_builtin") {
        fSiPMResponseMode = SensorResponseMode::BuiltinCurveRaw;
      } else if (m == "csv_curve_raw" || m == "csv_raw" || m == "raw_csv") {
        fSiPMResponseMode = SensorResponseMode::CsvCurveRaw;
      } else if (m == "csv_curve" || m == "csv" || m == "file") {
        fSiPMResponseMode = SensorResponseMode::CsvCurve;
      } else {
        fSiPMResponseMode = SensorResponseMode::BuiltinCurve;
      }
      fMPT_SiPM_Active = nullptr;
  }

  void SetPMTCurveFile(const G4String& v) { fPMTCurveFile = v; fMPT_TopPMT_Active = nullptr; fMPT_BotPMT_Active = nullptr; }
  void SetSiPMCurveFile(const G4String& v) { fSiPMCurveFile = v; fMPT_SiPM_Active = nullptr; }

  static G4String ToLower(const G4String& s) {
      G4String out = s;
      std::transform(out.begin(), out.end(), out.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      return out;
  }

  void SetInnerRadiusCm(G4double v) { fInnerRadius = v * CLHEP::cm; }
  void SetInnerHeightCm(G4double v) { fInnerHeight = v * CLHEP::cm; }
  void SetOuterDiameterCm(G4double v) { fOuterDiameter = v * CLHEP::cm; }
  void SetOuterHeightCm(G4double v) { fOuterHeight = v * CLHEP::cm; }
  void SetOuterScintillation(G4bool enabled) { fOuterScintillation = enabled; }
  G4bool GetOuterScintillation() const { return fOuterScintillation; }

  void SetNeutronDetectors(const G4String& value) {
      std::vector<NeutronDetectorConfig> detectors;
      std::stringstream entries(value);
      std::string entry;
      G4int index = 0;

      while (std::getline(entries, entry, ';')) {
          entry.erase(std::remove_if(entry.begin(), entry.end(), ::isspace), entry.end());
          if (entry.empty()) {
              continue;
          }

          std::vector<std::string> parts;
          std::stringstream fields(entry);
          std::string field;
          while (std::getline(fields, field, ':')) {
              parts.push_back(field);
          }

          if (parts.size() == 2) {
              parts.insert(parts.begin(), "A" + std::to_string(index));
          }
          if (parts.size() != 3) {
              G4cerr << "Ignoring malformed neutron detector entry: " << entry << G4endl;
              continue;
          }

          detectors.push_back({
              parts[0],
              std::stod(parts[1]) * CLHEP::deg,
              std::stod(parts[2]) * CLHEP::cm,
          });
          ++index;
      }

      fNeutronDetectors = detectors;
  }

  G4bool IsNeutronDetectorName(const G4String& name) const {
      for (const auto& detector : fNeutronDetectors) {
          if (detector.name == name) {
              return true;
          }
      }
      return false;
  }

  void SetPMTSideCm(G4double v) { fPMTSide = v * CLHEP::cm; }
  void SetSiPMTileCm(G4double v) { fSiPMTile = v * CLHEP::cm; }

  void SetTOFWindowNs(G4double minNs, G4double maxNs) {
      if (minNs >= maxNs) {
          G4cerr << "Ignoring invalid TOF window: " << minNs << " " << maxNs << G4endl;
          return;
      }
      fTOFMinNs = minNs;
      fTOFMaxNs = maxNs;
      G4cout << "Configured TOF window: [" << fTOFMinNs << ", " << fTOFMaxNs << "] ns" << G4endl;
  }

  G4double GetTOFMinNs() const { return fTOFMinNs; }
  G4double GetTOFMaxNs() const { return fTOFMaxNs; }

  void SetPrimaryEnergyMeV(G4double energyMeV) {
      if (energyMeV <= 0.0) {
          G4cerr << "Ignoring invalid primary energy: " << energyMeV << " MeV" << G4endl;
          return;
      }
      fPrimaryEnergy = energyMeV * CLHEP::MeV;
      G4cout << "Configured primary energy: " << energyMeV << " MeV" << G4endl;
  }

  G4double GetPrimaryEnergy() const { return fPrimaryEnergy; }

  void SetAbsLengthCm(G4double v) {
      fLArAbsLength = v * CLHEP::cm;
      for (int i = 0; i <= 2 && i < static_cast<int>(absorptionLength.size()); i++) {
          absorptionLength[i] = fLArAbsLength;
      }
      UpdateConfigurableOpticalProperties();
  }

  void SetBirksMmPerMeV(G4double v) {
      birks = v * CLHEP::mm / CLHEP::MeV;
      UpdateConfigurableOpticalProperties();
  }

  void SetScintYieldPerMeV(G4double v) { fLArScintYieldPerMeV = v; UpdateConfigurableOpticalProperties(); }
  void SetScintYieldScale(G4double v) { fLArScintYieldScale = v; UpdateConfigurableOpticalProperties(); }
  void SetIonScintYieldScale(G4double v) { fLArIonScintYieldScale = v; }
  G4double GetIonScintYieldScale() const { return fLArIonScintYieldScale; }
  G4bool IsInsideScintillatingLAr(const G4ThreeVector& position) const {
      // The LAr cylinders' local z axes are rotated onto global y.
      const G4double radial2 = position.x() * position.x() +
                               position.z() * position.z();
      const G4bool insideInner =
          radial2 <= fInnerRadius * fInnerRadius &&
          std::abs(position.y()) <= 0.5 * fInnerHeight;
      const G4bool insideOuter = fOuterScintillation &&
          radial2 <= 0.25 * fOuterDiameter * fOuterDiameter &&
          std::abs(position.y()) <= 0.5 * fOuterHeight;
      return insideInner || insideOuter;
  }
  void SetFastTimeNs(G4double v) { fLArFastTime = v * CLHEP::ns; UpdateConfigurableOpticalProperties(); }
  void SetSlowTimeNs(G4double v) { fLArSlowTime = v * CLHEP::ns; UpdateConfigurableOpticalProperties(); }
  void SetFastFraction(G4double v) { fLArFastFraction = v; UpdateConfigurableOpticalProperties(); }

  void SetTPBEfficiency(G4double v) { fTPBEfficiency = v; UpdateConfigurableOpticalProperties(); }
  void SetReflectorReflectivity(G4double v) { fReflectorReflectivity = v; fMPT_ESR = nullptr; }

public:

  G4double tpb;
  std::vector<G4VisAttributes*> fVisAttributes;

  explicit DetectorConstruction();
  virtual ~DetectorConstruction();

  void SetAbsorberMaterial (const G4String&);
  void SetAbsorberThickness(G4double);
  void SetAbsorberSizeYZ   (G4double);

  void UpdateOuterCellSize(G4double cellDiameter);
  void SetAbsorberXpos(G4double);

  void SetBirksConstant(G4double value) {
      birks = value;
      UpdateConfigurableOpticalProperties();
  }

  void SetAbsorptionLength(G4double value) {
      fLArAbsLength = value;
      for (int i = 0; i <= 2 && i < static_cast<int>(absorptionLength.size()); i++) {
          absorptionLength[i] = value;
      }
      UpdateConfigurableOpticalProperties();
  }

  void SetOuterCellSize(G4double value) {
      outerDiameter = value;
      fOuterDiameter = value;
  }

  const void placeLiquidScintillator(G4ThreeVector pos, G4double rotatAng, G4double fiberDiam, G4double fiberLen, G4Material* fibMat, G4LogicalVolume* motherVol, G4String name) {

    G4Tubs* solidFiber = new G4Tubs(name,                      // name
                           0 * CLHEP::mm, 0.5 * fiberDiam,     // r1, r2
                           0.5 * fiberLen,                     // half-length
                           0., CLHEP::twopi);                  // theta1, theta2

    G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber,  // solid
                                     fibMat,                      // material
                                     name);                       // name

    // rotation around y axis
    // u, v, w are the daughter axes, projected on the mother frame
    G4double rotatAng2 = CLHEP::pi / 2 - rotatAng;
    G4ThreeVector u = G4ThreeVector(std::cos(rotatAng2), 0, std::sin(rotatAng2));
    G4ThreeVector v = G4ThreeVector(0, 1, 0);
    G4ThreeVector w = G4ThreeVector(std::sin(-rotatAng2), 0, std::cos(rotatAng2));

    // rotation around x axis
    // G4ThreeVector u1 = G4ThreeVector(1, 0, 0);
    // G4ThreeVector v1 = G4ThreeVector(0, std::cos(rotatAng), std::sin(-rotatAng));
    // G4ThreeVector w1 = G4ThreeVector(0, std::sin(rotatAng), std::cos(rotatAng));

    G4RotationMatrix rotMat = G4RotationMatrix(u, v, w);
    std::cout << "position before rotation = " << pos << std::endl;
    pos = rotateAroundY(pos, rotatAng);
    std::cout << "position after rotation = " << pos << std::endl;
    const G4Transform3D& trans = G4Transform3D(rotMat, pos);

    G4VPhysicalVolume* physiFiber = new G4PVPlacement(trans,
                                  logicFiber,                  // its logical volume
                                  name,                        // its name
                                  motherVol,                   // its mother
                                  false,                       // no boolean operat
                                  0);

    auto visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));   // make liquid scintilator turquoise
    logicFiber->SetVisAttributes(visAttributes);
    visAttributes->SetForceSolid(true);
    fVisAttributes.push_back(visAttributes);

    // After: G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fibMat, name);

    // Create (once) and reuse a region for all LS volumes
    //static G4Region* gLSRegion = nullptr;
    //if (!gLSRegion) {
      //gLSRegion = new G4Region("LSRegion");
      //auto cuts = new G4ProductionCuts();
      //cuts->SetProductionCut(2.0*CLHEP::mm, G4ProductionCuts::GetIndex("gamma"));
      //cuts->SetProductionCut(2.0*CLHEP::mm, G4ProductionCuts::GetIndex("e-"));
      //cuts->SetProductionCut(2.0*CLHEP::mm, G4ProductionCuts::GetIndex("e+"));
      //cuts->SetProductionCut(2.0*CLHEP::mm, G4ProductionCuts::GetIndex("proton"));
      //gLSRegion->SetProductionCuts(cuts);
    //}
    //gLSRegion->AddRootLogicalVolume(logicFiber);

}

  void SetWorldMaterial(const G4String&);
  void SetWorldSizeX   (G4double);
  void SetWorldSizeYZ  (G4double);

  void SetMagField(G4double);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void PrintGeomParameters();

  const G4Material* GetAbsorberMaterial() const {return fAbsorberMaterial;};
  G4double GetAbsorberThickness() const         {return fAbsorberThickness;};
  G4double GetAbsorberSizeYZ() const            {return fAbsorberSizeYZ;};

  G4double GetAbsorberXpos() const              {return fXposAbs;};
  G4double GetxstartAbs() const                 {return fXstartAbs;};
  G4double GetxendAbs() const                   {return fXendAbs;};
  
  std::vector<G4double> absorptionLength;
  G4double birks;
  G4double outerDiameter;
  G4double LArX;
  G4double LArY;
  G4double LArZ;
  G4double innerDiameter;
  G4double innerHeight;
  G4double outerHeight;
  G4double           fiberDiameter;
  G4double           fiberLength;
  G4Material *       fiberMat;
  G4Material*        tpbMat;


  G4Material*       GetTPBMaterial()   const { return tpbMat; }
  G4LogicalVolume*  GetInnerCellLV()   const { return innerCellLogic; }

  const G4Material* GetWorldMaterial() const    {return fWorldMaterial;};
  G4double GetWorldSizeX() const                {return fWorldSizeX;};

  const G4VPhysicalVolume* GetAbsorber() const  {return innerCellPhysi;};

  const G4ThreeVector & MatVectMul(G4RotationMatrix mat,const G4ThreeVector & vect){

      G4double x = mat.xx()*vect.x() + mat.xy()*vect.y()+mat.xz()*vect.z();
      G4double y = mat.yx()*vect.x() + mat.yy()*vect.y()+mat.yz()*vect.z();
      G4double z = mat.zx()*vect.x() + mat.zy()*vect.y()+mat.zz()*vect.z();

      const G4ThreeVector & ret = G4ThreeVector(x,y,z);
      std::cout << "ret = " << ret << std::endl;
      return ret;
  }

  const G4ThreeVector & rotateAroundY(const G4ThreeVector & v,G4double angle){

      G4ThreeVector x = G4ThreeVector(std::cos(angle), 0, std::sin(-angle));
      G4ThreeVector y = G4ThreeVector(0, 1, 0);
      G4ThreeVector z = G4ThreeVector(std::sin(angle),0,std::cos(angle));
      G4RotationMatrix rotm1  = G4RotationMatrix(x, y, z);

      const G4ThreeVector & ret = MatVectMul(rotm1,v);
      return ret;
  }
 

private:

  void DefineMaterials();
  void ComputeGeomParameters();
  void ChangeGeometry();

  G4Material*        Al;
  G4Material*        noScintMaterial;
  G4Material*        fAbsorberMaterial;
  G4bool             fOuterScintillation = false;
  G4double           fAbsorberThickness;
  G4double           fAbsorberSizeYZ;

  G4double           fXposAbs;
  G4double           fXstartAbs, fXendAbs;

  G4Material*        fWorldMaterial;
  G4double           fWorldSizeX;
  G4double           fWorldSizeYZ;

  G4Box*             fSolidWorld;
  G4LogicalVolume*   fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  G4Box*             airSolid;
  G4LogicalVolume*   airLogic;
  G4VPhysicalVolume* airPhysi;

  G4Box*             hydroSolid;
  G4LogicalVolume*   hydroLogic;
  G4VPhysicalVolume* hydroPhysi;

  G4Tubs*            outerCellSolid;
  G4LogicalVolume*   outerCellLogic;
  G4VPhysicalVolume* outerCellPhysi;

  G4Tubs*            innerLayerSolid;
  G4LogicalVolume*   innerLayerLogic;
  G4VPhysicalVolume*  innerLayerPhysi;

  G4Tubs*            outerLayerOneSolid;
  G4LogicalVolume*   outerLayerOneLogic;
  G4VPhysicalVolume* outerLayerOnePhysi;

  G4Tubs*            outerLayerTwoSolid;
  G4LogicalVolume*   outerLayerTwoLogic;
  G4VPhysicalVolume* outerLayerTwoPhysi;

  G4Box*            fSolidTeflon;
  G4LogicalVolume*   fLogicTeflon;
  G4VPhysicalVolume* fPhysiTeflon;

 // G4Box*             fSolidAbsorber;
  G4Tubs*            innerCellSolid;
  G4LogicalVolume*   innerCellLogic;
  G4VPhysicalVolume* innerCellPhysi;
 
  G4Tubs*          svol_fiber;
  G4LogicalVolume* lvol_fiber;
  G4VPhysicalVolume* pvol_fiber;
 
  G4Tubs*          new_svol_fiber;
  G4LogicalVolume* new_lvol_fiber;
  G4VPhysicalVolume* new_pvol_fiber;

   
  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
