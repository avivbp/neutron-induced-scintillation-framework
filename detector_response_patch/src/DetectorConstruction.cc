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
/// \file electromagnetic/TestEm5/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Run.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4UniformMagField.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Exception.hh"
#include <map>
#include <string>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cctype>
//#include "G4GeometryTolerance.hh"

namespace {
  constexpr int kNESR = 46;

  G4double kESR_E[kNESR] = {
  2.101*eV, 2.142*eV, 2.180*eV, 2.214*eV, 2.246*eV,
  2.269*eV, 2.297*eV, 2.335*eV, 2.360*eV, 2.385*eV,
  2.406*eV, 2.438*eV, 2.464*eV, 2.503*eV, 2.526*eV,
  2.549*eV, 2.585*eV, 2.615*eV, 2.646*eV, 2.684*eV,
  2.717*eV, 2.765*eV, 2.806*eV, 2.849*eV, 2.886*eV,
  2.955*eV, 3.044*eV, 3.121*eV, 3.157*eV, 3.175*eV,
  3.193*eV, 3.221*eV, 3.239*eV, 3.249*eV, 3.268*eV,
  3.278*eV, 3.287*eV, 3.297*eV, 3.307*eV, 3.327*eV,
  3.347*eV, 3.377*eV, 3.408*eV, 3.450*eV, 3.494*eV,
  3.538*eV
  };

  G4double kESR_R[kNESR] = {
  0.979, 0.982, 0.982, 0.982, 0.982,
  0.982, 0.982, 0.982, 0.982, 0.982,
  0.982, 0.982, 0.982, 0.982, 0.982,
  0.982, 0.979, 0.979, 0.979, 0.977,
  0.977, 0.977, 0.977, 0.977, 0.979,
  0.977, 0.972, 0.958, 0.944, 0.925,
  0.896, 0.816, 0.721, 0.633, 0.562,
  0.493, 0.422, 0.355, 0.296, 0.277,
  0.261, 0.246, 0.234, 0.223, 0.220,
  0.204
};


 
  namespace {
    G4double SafeDenominatorForReflectivity(G4double refl)
    {
      refl = std::clamp(refl, 0.0, 0.999999);
      return std::max(1.0 - refl, 1.0e-9);
    }

    G4double CorrectForGeant4DetectionBranching(G4double physicalPDE, G4double refl)
    {
      physicalPDE = std::clamp(physicalPDE, 0.0, 1.0);
      const G4double corrected = physicalPDE / SafeDenominatorForReflectivity(refl);
      return std::clamp(corrected, 0.0, 1.0);
    }

    G4String LowerCopy(const G4String& s)
    {
      G4String out = s;
      std::transform(out.begin(), out.end(), out.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      return out;
    }

    bool TryParseDouble(const std::string& token, G4double& value)
    {
      try {
        size_t pos = 0;
        value = std::stod(token, &pos);
        return pos > 0;
      }
      catch (...) {
        return false;
      }
    }

    std::vector<std::string> SplitCSVLine(const std::string& line)
    {
      std::vector<std::string> out;
      std::stringstream ss(line);
      std::string item;
      while (std::getline(ss, item, ',')) {
        // trim leading/trailing spaces
        const auto first = item.find_first_not_of(" \t\r\n");
        const auto last  = item.find_last_not_of(" \t\r\n");
        out.push_back(first == std::string::npos ? std::string() : item.substr(first, last - first + 1));
      }
      return out;
    }

    G4MaterialPropertiesTable* BuildCurveMPT(
        const std::vector<G4double>& energies,
        const std::vector<G4double>& physicalEfficiencies,
        G4double refl)
    {
      auto mpt = new G4MaterialPropertiesTable();

      if (energies.empty() || physicalEfficiencies.empty() || energies.size() != physicalEfficiencies.size()) {
        G4Exception("BuildCurveMPT", "DetectorResponse001", JustWarning,
                    "Empty or malformed sensor response curve. Falling back to zero efficiency.");
        G4double E2[2] = {2.0*eV, 3.5*eV};
        G4double Eff[2] = {0.0, 0.0};
        G4double Refl[2] = {std::clamp(refl, 0.0, 1.0), std::clamp(refl, 0.0, 1.0)};
        mpt->AddProperty("EFFICIENCY", E2, Eff, 2);
        mpt->AddProperty("REFLECTIVITY", E2, Refl, 2);
        return mpt;
      }

      std::vector<std::pair<G4double, G4double>> points;
      for (size_t i = 0; i < energies.size(); ++i) {
        points.emplace_back(energies[i], std::clamp(physicalEfficiencies[i], 0.0, 1.0));
      }
      std::sort(points.begin(), points.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });

      std::vector<G4double> e(points.size());
      std::vector<G4double> effCorrected(points.size());
      for (size_t i = 0; i < points.size(); ++i) {
        e[i] = points[i].first;
        effCorrected[i] = CorrectForGeant4DetectionBranching(points[i].second, refl);
      }

      G4double E2[2] = {e.front(), e.back()};
      G4double Refl[2] = {std::clamp(refl, 0.0, 1.0), std::clamp(refl, 0.0, 1.0)};
      mpt->AddProperty("EFFICIENCY", e.data(), effCorrected.data(), static_cast<G4int>(e.size()));
      mpt->AddProperty("REFLECTIVITY", E2, Refl, 2);
      return mpt;
    }

    G4double PeakOrOne(const std::vector<G4double>& v)
    {
      if (v.empty()) return 1.0;
      const auto it = std::max_element(v.begin(), v.end());
      return (*it > 0.0) ? *it : 1.0;
    }
  }

  G4MaterialPropertiesTable* MakeFlatMPT(G4double eff, G4double refl)
  {
    G4double E2[2] = {2.0*eV, 3.5*eV};

    auto mpt = new G4MaterialPropertiesTable();
    G4double Eff[2]  = {CorrectForGeant4DetectionBranching(eff, refl),
                        CorrectForGeant4DetectionBranching(eff, refl)};
    G4double Refl[2] = {std::clamp(refl, 0.0, 1.0), std::clamp(refl, 0.0, 1.0)};
    mpt->AddProperty("EFFICIENCY",   E2, Eff,  2);
    mpt->AddProperty("REFLECTIVITY", E2, Refl, 2);
    return mpt;
  }

  G4MaterialPropertiesTable* MakePMTMPT(G4double refl, G4double requestedPeakPDE)
  {
    G4int nPoints = 28;

    G4double energy[nPoints] = {
      1.98*eV, 2.00*eV, 2.01*eV, 2.02*eV, 2.05*eV, 2.09*eV, 2.12*eV,
      2.17*eV, 2.20*eV, 2.24*eV, 2.30*eV, 2.33*eV, 2.37*eV, 2.38*eV,
      2.41*eV, 2.45*eV, 2.51*eV, 2.63*eV, 2.69*eV, 2.73*eV, 2.80*eV,
      2.87*eV, 2.96*eV, 3.06*eV, 3.19*eV, 3.30*eV, 3.44*eV, 3.54*eV
    };

    std::vector<G4double> efficiency = {
      0.0158, 0.0184, 0.0212, 0.0241, 0.0295, 0.0373, 0.0450,
      0.0562, 0.0655, 0.0754, 0.0932, 0.109,  0.126,  0.136,
      0.153,  0.172,  0.184,  0.215,  0.228,  0.233,  0.244,
      0.256,  0.268,  0.271,  0.275,  0.281,  0.281,  0.281
    };

    const G4double scale = (requestedPeakPDE > 0.0) ? requestedPeakPDE / PeakOrOne(efficiency) : 1.0;
    for (auto& x : efficiency) x = std::clamp(scale * x, 0.0, 1.0);

    std::vector<G4double> energies(energy, energy + nPoints);
    return BuildCurveMPT(energies, efficiency, refl);
  }

  G4MaterialPropertiesTable* MakeSiPMMPT(G4double refl, G4double requestedPeakPDE)
  {
    G4int nPoints = 32;

    G4double energy[nPoints] = {
      2.00*eV, 2.03*eV, 2.07*eV, 2.11*eV, 2.14*eV, 2.18*eV, 2.21*eV, 2.26*eV,
      2.30*eV, 2.32*eV, 2.36*eV, 2.40*eV, 2.43*eV, 2.46*eV, 2.49*eV, 2.56*eV,
      2.63*eV, 2.71*eV, 2.78*eV, 2.87*eV, 2.96*eV, 3.05*eV, 3.12*eV, 3.17*eV,
      3.22*eV, 3.22*eV, 3.26*eV, 3.32*eV, 3.38*eV, 3.42*eV, 3.46*eV, 3.50*eV
    };

    std::vector<G4double> efficiency = {
      0.241, 0.253, 0.264, 0.279, 0.291, 0.305, 0.315, 0.329,
      0.344, 0.352, 0.360, 0.370, 0.378, 0.383, 0.388, 0.395,
      0.401, 0.402, 0.398, 0.390, 0.378, 0.361, 0.351, 0.340,
      0.323, 0.315, 0.305, 0.286, 0.271, 0.256, 0.242, 0.229
    };

    const G4double scale = (requestedPeakPDE > 0.0) ? requestedPeakPDE / PeakOrOne(efficiency) : 1.0;
    for (auto& x : efficiency) x = std::clamp(scale * x, 0.0, 1.0);

    std::vector<G4double> energies(energy, energy + nPoints);
    return BuildCurveMPT(energies, efficiency, refl);
  }

  G4MaterialPropertiesTable* MakeCurveMPTFromFile(const G4String& fileName, G4double refl, G4double requestedPeakPDE)
  {
    std::ifstream in(fileName);
    if (!in) {
      G4ExceptionDescription msg;
      msg << "Could not open sensor response curve file: " << fileName
          << ". Falling back to a flat response using the requested PDE.";
      G4Exception("MakeCurveMPTFromFile", "DetectorResponse002", JustWarning, msg);
      return MakeFlatMPT(requestedPeakPDE, refl);
    }

    std::vector<G4double> energies;
    std::vector<G4double> qes;
    bool headerSaysWavelength = false;
    bool headerSaysEnergy = false;

    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      auto parts = SplitCSVLine(line);
      if (parts.size() < 2) continue;

      G4double x = 0.0, qe = 0.0;
      const bool numeric = TryParseDouble(parts[0], x) && TryParseDouble(parts[1], qe);
      if (!numeric) {
        const std::string h0 = parts[0];
        const std::string h1 = parts[1];
        std::string header = h0 + "," + h1;
        std::transform(header.begin(), header.end(), header.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        headerSaysWavelength = header.find("wave") != std::string::npos || header.find("nm") != std::string::npos;
        headerSaysEnergy = header.find("energy") != std::string::npos || header.find("ev") != std::string::npos;
        continue;
      }

      // Accept qe either as fraction or percent.
      if (qe > 1.0) qe /= 100.0;

      G4double energy = 0.0;
      if (headerSaysWavelength || (!headerSaysEnergy && x > 20.0)) {
        // First column is wavelength in nm.
        energy = (1239.841984 / x) * eV;
      } else {
        // First column is photon energy in eV.
        energy = x * eV;
      }

      energies.push_back(energy);
      qes.push_back(std::clamp(qe, 0.0, 1.0));
    }

    if (energies.empty()) {
      G4ExceptionDescription msg;
      msg << "Sensor response curve file contained no valid points: " << fileName
          << ". Falling back to a flat response using the requested PDE.";
      G4Exception("MakeCurveMPTFromFile", "DetectorResponse003", JustWarning, msg);
      return MakeFlatMPT(requestedPeakPDE, refl);
    }

    if (requestedPeakPDE > 0.0) {
      const G4double scale = requestedPeakPDE / PeakOrOne(qes);
      for (auto& qe : qes) qe = std::clamp(scale * qe, 0.0, 1.0);
    }

    return BuildCurveMPT(energies, qes, refl);
  }

  G4MaterialPropertiesTable* MakeESRMPT()
  {
    auto mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", kESR_E, kESR_R, kNESR);

    // Explicitly set efficiency to 0 across the same grid
    std::vector<G4double> Eff(kNESR, 0.0);
    mpt->AddProperty("EFFICIENCY", kESR_E, Eff.data(), kNESR);

    // If you use unified model and want specular behavior:
    mpt->AddConstProperty("SPECULARSPIKECONSTANT", 1.0, true);
    mpt->AddConstProperty("SPECULARLOBECONSTANT",  0.0, true);
    mpt->AddConstProperty("BACKSCATTERCONSTANT",   0.0, true);

    return mpt;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberMaterial(nullptr),fWorldMaterial(nullptr),
   fSolidWorld(nullptr),fLogicWorld(nullptr),fPhysiWorld(nullptr),
   innerCellSolid(nullptr),innerCellLogic(nullptr),innerCellPhysi(nullptr),
   outerCellSolid(nullptr),outerCellLogic(nullptr),outerCellPhysi(nullptr),
   innerLayerSolid(nullptr),innerLayerLogic(nullptr),innerLayerPhysi(nullptr),
   outerLayerOneSolid(nullptr),outerLayerOneLogic(nullptr),outerLayerOnePhysi(nullptr),
   outerLayerTwoSolid(nullptr),outerLayerTwoLogic(nullptr),outerLayerTwoPhysi(nullptr),
   fDetectorMessenger(nullptr),lvol_fiber(nullptr),pvol_fiber(nullptr)
{
  // default parameter values of the calorimeter
  fAbsorberThickness = 75.*mm;
  fAbsorberSizeYZ    = 47.*mm;
  fXposAbs           = 0.*cm;
  fiberDiameter = 12.7*cm;
  fiberLength = fiberDiameter;
  ComputeGeomParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   ("G4_Galactic");
  SetAbsorberMaterial("G4_lAr");
  //SetAbsorberMaterial("G4_C");

  auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  //auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  std::vector<G4double> energy     = {9.5*eV,9.69*eV,9.89*eV,10.0*eV};
  // refractive index from Collins [https://indico.cern.ch/event/609440/contributions/2548164/attachments/1441824/2220090/Light_Propagation_in_Liquid_Argon.pdf]
  std::vector<G4double> rindex     = {1.4,1.42,1.43,1.44};
  // taking the value of 1ppm C impurities in lAr with cross section 1Mbarn per molecule
  // absorption length falls linearly with impurity concentration , cross section per molecule
  G4double absLen = fLArAbsLength;
  absorptionLength = {absLen, absLen, absLen, absLen};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // property independent of energy

  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  tpb = 0.266;
  // Doke et al. , NIMA 291,3(1990) [https://www.sciencedirect.com/science/article/pii/016890029090011T?via%3Dihub] 
  MPT->AddConstProperty("SCINTILLATIONYIELD", fLArScintYieldScale * fLArScintYieldPerMeV / MeV);
  MPT->AddProperty("RINDEX", energy, rindex);
  MPT->AddProperty("ABSLENGTH", energy, absorptionLength);

  // fano factor of LAr = 0.11 , Doke et al, NIM 134 (1976)353, taken from [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf]

  // The actual number of emitted photons during a step
  //fluctuates around the mean number of photons with
  //a width given by ResolutionScale x sqrt(MeanNumberOfPhotons)
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);

  MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", fLArFastTime);
  // we take an average of the following for the slow time components:
  // [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf] - 1400ns
  // [https://indico.cern.ch/event/609440/contributions/2548164/attachments/1441824/2220090/Light_Propagation_in_Liquid_Argon.pdf] - 1500ns
  // [Creus et al. , JINST 10 (2015) ] - 1600ns
  MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2", fLArSlowTime);
  // values taken from [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf] - yield fraction = 0.75 in old geant4 version
  MPT->AddConstProperty("SCINTILLATIONYIELD1", fLArFastFraction);
  MPT->AddConstProperty("SCINTILLATIONYIELD2", 1.0 - fLArFastFraction);
  std::vector<G4double> eVUV   = {9.50*eV, 9.65*eV, 9.70*eV, 9.80*eV, 9.90*eV};
  std::vector<G4double> relVUV = {0.00,    0.50,    1.00,    0.50,    0.00   };
  MPT->AddProperty("SCINTILLATIONCOMPONENT1", eVUV, relVUV);
  MPT->AddProperty("SCINTILLATIONCOMPONENT2", eVUV, relVUV);

  lAr->SetMaterialPropertiesTable(MPT);

  // Set the Birks Constant for the LAr scintillator, using C. William thesis [http://unizh.web.cern.ch/Publications/Articles/thesis_William.pdf]
  //lAr->GetIonisation()->SetBirksConstant(0.053* mm / MeV);
  // trying using Doke et al. 1990 paper
 
 //lAr->GetIonisation()->SetBirksConstant(0.89* mm / MeV);
 lAr->GetIonisation()->SetBirksConstant(birks);
  std::cout << "lAr Ion Pair Mean Energy = " << lAr->GetIonisation()->GetMeanEnergyPerIonPair()/CLHEP::eV << " eV" << std::endl;
  if(innerCellLogic) { innerCellLogic->SetMaterial(lAr); }

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  //auto nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  //nav->SetPushThreshold(1e-4*CLHEP::mm); // 0.1 µm-ish nudge

  //SetAbsorberMaterial("G4_lAr");
 //   fLogicAbsorber->SetMaterial(lAr);
    //G4RunManager::GetRunManager()->PhysicsHasBeenModified();
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this,run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

void DetectorConstruction::UpdateOuterCellSize(G4double d) {
  outerDiameter = d;
  auto rm = G4RunManager::GetRunManager();
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()       ->Clean();
  rm->DefineWorldVolume( Construct() );
  rm->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;
  
  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

  
  std::vector<G4double> energy3     = {1.5*eV,9.69*eV,9.89*eV,10.0*eV};
  std::vector<G4double> rindex3     = {1.4,1.42,1.43,1.44};
  G4double absLen3 = 1.5*m;
  std::vector<G4double> absorptionLength3 = {absLen3, absLen3, absLen3,absLen3};


  // Material properties (transport only)
  auto* lAr_noScint = new G4Material("LAr_NoScint", 1.390*g/cm3, 1);
  lAr_noScint->AddElement(Ar, 1);

  auto* mptNoScint = new G4MaterialPropertiesTable();
  mptNoScint->AddProperty("RINDEX", energy3, rindex3);
  mptNoScint->AddProperty("ABSLENGTH", energy3, absorptionLength3);
  lAr_noScint->SetMaterialPropertiesTable(mptNoScint);
  noScintMaterial = lAr_noScint;

  Al = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  auto* mptAL = new G4MaterialPropertiesTable();
  std::vector<G4double> energies = {1.5*eV,10.0*eV};
  std::vector<G4double> rindices = {1.0,1.0};
  std::vector<G4double> absLens = {7.*nm,7.*nm};
  mptAL->AddProperty("RINDEX",energies,rindices);
  mptAL->AddProperty("ABSLENGTH",energies,absLens);
  Al->SetMaterialPropertiesTable(mptAL);

  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = 
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  auto airMPT = new G4MaterialPropertiesTable();

  // Cover your visible band (adjust if you need wider/narrower)
  const G4double eMin = 1.5*CLHEP::eV;   // ~ 830 nm
  const G4double eMax = 10.2*CLHEP::eV;   // ~ 100 nm
  G4double ephoton[2] = { eMin, eMax };

  // Dry air ~1.00027 in visible; using a flat value is fine
  G4double rindexAir[2] = { 1.00027, 1.00027 };

  airMPT->AddProperty("RINDEX", ephoton, rindexAir, 2);
  Air->SetMaterialPropertiesTable(airMPT);

  // Scintillator

  G4Material* Sci =
       new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
     
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  std::vector<G4double> energy     = {9.5*eV,9.6*eV,10.0*eV};
  std::vector<G4double> rindex     = {1.42,1.43,1.44};
  std::vector<G4double> absorptionLength2 = {1500.*cm, 1500.*cm, 1500.*cm};

  G4MaterialPropertiesTable* MPT2 = new G4MaterialPropertiesTable();


  MPT2->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
  MPT2->AddProperty("RINDEX", energy, rindex);
  MPT2->AddProperty("ABSLENGTH", energy, absorptionLength2);


  MPT2->AddConstProperty("RESOLUTIONSCALE",1.0);

  MPT2->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 7. * ns);
  MPT2->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  MPT2->AddProperty("SCINTILLATIONCOMPONENT1", energy, { 1.0, 1.0, 1.0 });

  Sci->SetMaterialPropertiesTable(MPT2);
  fiberMat = Sci;

  auto TPB  = new G4Material("TPB", 1.20*g/cm3, 2);   // ~1.2 g/cc
  TPB->AddElement(C, 28);
  TPB->AddElement(H, 22);
  
  // --- TPB Material Properties Table (ascending energies; up to 10 eV) ---
  auto mptTPB = new G4MaterialPropertiesTable();
  
  // 1) Refractive index (kept flat here; supply dispersion if you have it)
  {
    static const G4int N = 12;
    static G4double E[N] = {
      1.50*eV, 2.00*eV, 2.50*eV, 3.00*eV, 3.50*eV, 4.00*eV,
      5.00*eV, 6.00*eV, 7.00*eV, 8.00*eV, 9.00*eV, 10.00*eV
    }; // strictly ascending (≈827 nm → 124 nm)
  
    static G4double n[N] = {
      1.65, 1.65, 1.65, 1.65, 1.65, 1.65,
      1.65, 1.65, 1.65, 1.65, 1.65, 1.65
    };
    mptTPB->AddProperty("RINDEX", E, n, N);
  }
  
  // 2) WLS absorption length (short in VUV; long in visible)
  //    Controls P_absorb together with your TPB film thickness.
  {
    static const G4int N = 12;
    static G4double E[N] = {
      1.50*eV, 2.00*eV, 2.50*eV, 3.00*eV, 3.50*eV, 4.00*eV,
      5.00*eV, 6.00*eV, 7.00*eV, 8.00*eV, 9.70*eV, 10.00*eV
    }; // ascending up to ≥ 128 nm (9.7 eV) and 124 nm (10 eV)
  
    static G4double L[N] = {
      10.*m, 10.*m, 10.*m, 10.*m, 10.*m, 5.*m,   // visible/near-UV ≫ film thickness
      1.*m,  1.*mm, 0.20*um, 0.18*um, 0.12*um, 0.10*um  // strong VUV absorption
    };
    mptTPB->AddProperty("WLSABSLENGTH", E, L, N);
  }
  
  // 3) WLS emission spectrum (visible band; zeros elsewhere)
  //    Shape peaked ~420–430 nm; values are relative (Geant4 normalizes internally).
  {
    // TPB emission (relative)
    G4int nTPB = 38;
    G4double tpbE[nTPB] = {
      2.304*eV, 2.323*eV, 2.345*eV, 2.364*eV, 2.388*eV, 2.410*eV, 2.429*eV, 2.453*eV, 2.476*eV, 2.498*eV,
      2.522*eV, 2.545*eV, 2.572*eV, 2.594*eV, 2.620*eV, 2.647*eV, 2.672*eV, 2.700*eV, 2.726*eV, 2.755*eV,
      2.785*eV, 2.815*eV, 2.844*eV, 2.873*eV, 2.905*eV, 2.952*eV, 2.999*eV, 3.033*eV, 3.069*eV, 3.106*eV,
      3.140*eV, 3.179*eV, 3.212*eV, 3.256*eV, 3.297*eV, 3.380*eV, 3.417*eV, 3.467*eV
    };
    G4double tpbWLS[nTPB] = {
      0.020, 0.063, 0.072, 0.081, 0.120, 0.133, 0.173, 0.203, 0.183, 0.171,
      0.216, 0.259, 0.290, 0.323, 0.373, 0.427, 0.556, 0.558, 0.632, 0.709,
      0.786, 0.832, 0.869, 0.909, 0.952, 1.000, 0.933, 0.816, 0.685, 0.486,
      0.278, 0.102, 0.053, 0.017, 0.050, 0.026, 0.011, 0.014
    };
    mptTPB->AddProperty("WLSCOMPONENT", tpbE, tpbWLS, nTPB);
  }
  
  // 4) WLS timing + intrinsic yield (“QE”)
  //    Set the mean photons emitted per photon absorbed (≤ 1 for PLQY).
  mptTPB->AddConstProperty("WLSTIMECONSTANT", 1.7*ns);
  mptTPB->AddConstProperty("WLSMEANNUMBERPHOTONS", fTPBEfficiency);
  
  TPB->SetMaterialPropertiesTable(mptTPB);

  tpbMat = TPB;
  
    //
  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeGeomParameters()
{
  // Compute derived parameters of the calorimeter
  fXstartAbs = fXposAbs-0.5*fAbsorberThickness; 
  fXendAbs   = fXposAbs+0.5*fAbsorberThickness;

  G4double xmax = std::max(std::abs(fXstartAbs), std::abs(fXendAbs)) +0.73*m+fiberLength;
  std::cout << "xmaxm = " << xmax << std::endl;
 // fWorldSizeX = 2.4*xmax; 
  //fWorldSizeYZ= fWorldSizeX;
  fWorldSizeX = fWorldSizeYZ =  20.*m;

  if(nullptr != fPhysiWorld) { ChangeGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{ 


  if(nullptr != fPhysiWorld) { return fPhysiWorld; }
  // World
  //
  std::cout << "xxxxxxxxxxxxxxxxxxx" << " world size = " << G4ThreeVector(fWorldSizeX,fWorldSizeYZ,fWorldSizeYZ)/CLHEP::cm << " cm" << std::endl; 

  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);  //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,          //its solid
                                   fWorldMaterial,        //its material
                                   "World");              //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(0.,0.,0.), //at (0,0,0)
                                 fLogicWorld,             //its logical volume
                                 "World",                 //its name
                                 0,                       //its mother  volume
                                 false,                   //no boolean operation
                                 0);                      //copy number
   
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;
  
  G4int ncomponents, natoms;
  G4double fractionmass;

  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4double airSizeX = fWorldSizeX/2-1.*m;
  G4double airSizeYZ = fWorldSizeYZ/2-1.*m;


  airSolid = new G4Box("air",
                      airSizeX,airSizeYZ,airSizeYZ);                              
  airLogic = new G4LogicalVolume(airSolid,    //its solid
                                Air,
                                "air");

  airPhysi = new G4PVPlacement(0,                      //no rotation
                              G4ThreeVector(0,0.,0.), 
                              airLogic,             //its logical volume
                             "air",                 //its name
                              fLogicWorld,                       //its mother  volume
                             false,                   //no boolean operation
                             0);


 
  auto* Teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
  std::vector<G4double> energy     = {9.69*eV,9.89*eV,10.0*eV};
  std::vector<G4double> rindex     = {1.35,1.37,1.38};
  G4MaterialPropertiesTable* MPT3 = new G4MaterialPropertiesTable();

  MPT3->AddProperty("RINDEX", energy, rindex);
  Teflon->SetMaterialPropertiesTable(MPT3);
  Teflon->GetMaterialPropertiesTable()->DumpTable();

  //fSolidTeflon = new G4Box("Teflon",0.5*LArX+1.*cm,0.5*LArY+1.*cm,0.5*LArZ+1.*cm);
    
  //fLogicTeflon = new G4LogicalVolume(fSolidTeflon,    //its solid
  //                             Teflon,                //its material
    //                             "Teflon");       //its name
    
  G4ThreeVector u = G4ThreeVector(1 , 0, 0);         
  G4ThreeVector v = G4ThreeVector(0 , 0, -1);        
  G4ThreeVector w = G4ThreeVector(0, 1, 0);
    
  G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
  G4ThreeVector position1 = G4ThreeVector(0,0,0);
  const G4Transform3D & transform1 = G4Transform3D(rotm1,position1);
    
  //fPhysiTeflon = new G4PVPlacement(transform1,
    //                 fLogicTeflon,     //its logical volume
      //                 "Teflon",         //its name
        //                 airLogic,        //its mother
          //                 false,              //no boulean operat
            //                 0);                 //copy number


  std::cout << "after rotm1" << std::endl;

  // outer cell
  outerDiameter = fOuterDiameter;
  outerHeight = fOuterHeight;

  outerCellSolid = new G4Tubs("outerCell",0., outerDiameter/2, outerHeight/2, 0., 360. * deg);

  outerCellLogic = new G4LogicalVolume(outerCellSolid,    //its solid
                                       noScintMaterial, //its material
                                       "outerCell");       //its name

  outerCellPhysi = new G4PVPlacement(transform1,
                                       outerCellLogic,     //its logical volume
                                       "outerCell",         //its name
                                       airLogic,        //its mother
                                       false,              //no boulean operat
                                       0);                 //copy number

  std::cout << "after outer cell" << std::endl;

  auto* steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  outerLayerOneSolid = new G4Tubs("outerLayerOne",outerDiameter/2, outerDiameter/2 + 1.5*mm, outerHeight/2, 0., 360. * deg);

  outerLayerOneLogic = new G4LogicalVolume(outerLayerOneSolid,    //its solid
                                         steel, //its material
                                         "outerLayerOne");       //its name

  outerLayerOnePhysi = new G4PVPlacement(transform1,
                                       outerLayerOneLogic,     //its logical volume
                                       "outerLayerOne",         //its name
                                       airLogic,        //its mother
                                       false,              //no boulean operat
                                       0); 
  
  std::cout << "after outer layer one" << std::endl;

  outerLayerTwoSolid = new G4Tubs("outerLayerTwo",outerDiameter/2 + 1.15*cm, outerDiameter/2 + 1.3*cm, outerHeight/2, 0., 360. * deg);

  outerLayerTwoLogic = new G4LogicalVolume(outerLayerTwoSolid,    //its solid
                                           steel, //its material
                                           "outerLayerTwo");       //its name

  outerLayerTwoPhysi = new G4PVPlacement(transform1,
                                           outerLayerTwoLogic,     //its logical volume
                                           "outerLayerTwo",         //its name
                                           airLogic,        //its mother
                                           false,              //no boulean operat
                                           0);

  std::cout << "after outer layer two" << std::endl;

  // inner cell
  innerDiameter = 2.0 * fInnerRadius;
  innerHeight = fInnerHeight;

  innerCellSolid = new G4Tubs("innerCell",0., innerDiameter/2, innerHeight/2, 0., 360. * deg);

  innerCellLogic = new G4LogicalVolume(innerCellSolid,    //its solid
                                       fAbsorberMaterial, //its material
                                       "innerCell");       //its name

  innerCellPhysi = new G4PVPlacement(0,                      //no rotation
                                     G4ThreeVector(0,0.,0.),
                                     innerCellLogic,     //its logical volume
                                     "innerCell",         //its name
                                     outerCellLogic,        //its mother
                                     false,              //no boulean operat
                                     0);                 //copy number


  std::cout << "after inner cell" << std::endl;

  innerLayerSolid = new G4Tubs("innerLayer",innerDiameter/2, innerDiameter/2 + 1.5*mm, innerHeight/2, 0., 360. * deg);

  innerLayerLogic = new G4LogicalVolume(innerLayerSolid,    //its solid
                                        Al, //its material
                                        "innerLayer");       //its name

  innerLayerPhysi = new G4PVPlacement(0,
                                      G4ThreeVector(),
                                      innerLayerLogic,     //its logical volume
                                      "innerLayer",         //its name
                                      outerCellLogic,        //its mother
                                      false,              //no boulean operat
                                      0);

  std::cout << "after inner layer" << std::endl;


  G4OpticalSurface* OpSurface = new G4OpticalSurface("name");
  //G4LogicalBorderSurface* Surface = new G4LogicalBorderSurface("name",innerCellPhysi,innerLayerPhysi,OpSurface);

  OpSurface->SetType(dielectric_metal);
  OpSurface->SetModel(unified);
  OpSurface->SetFinish(polished);
  std::vector<G4double> pp = {1.5*eV, 4.0*eV,7.0*eV,10.0*eV};
  //rindex = {1.35, 1.40};
  std::vector<G4double> reflectivity = {1.0, 1.0,1.0,1.0};
  //reflectivity = {0.01,0.01,0.01,0.01};

  G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

  //SMPT->AddProperty("RINDEX", pp, rindex);
  SMPT->AddProperty("REFLECTIVITY", pp, reflectivity);

  OpSurface->SetMaterialPropertiesTable(SMPT);
  //new G4LogicalSkinSurface("skinInnerLayer",innerLayerLogic,OpSurface);

  // Top/Bottom optical planes (same radius as inner cell), placed flush
  auto planeMat = Al; // any bulk; detection is on surface
  auto r = innerDiameter/2;
  auto t = 1.5*mm;

  auto topSolid = new G4Tubs("TopPlane", 0, r, t/2, 0, CLHEP::twopi);
  auto botSolid = new G4Tubs("BotPlane", 0, r, t/2, 0, CLHEP::twopi);
  auto topLV = new G4LogicalVolume(topSolid, planeMat, "TopPlaneLV");
  auto botLV = new G4LogicalVolume(botSolid, planeMat, "BotPlaneLV");

  auto topPV = new G4PVPlacement(nullptr, {0, 0, +innerHeight/2 + t/2}, topLV, "TopPlanePV", outerCellLogic, false, 0);
  auto botPV = new G4PVPlacement(nullptr, {0, 0, -innerHeight/2 - t/2}, botLV, "BotPlanePV", outerCellLogic, false, 0);

  new G4LogicalSkinSurface("skinTop",topLV,OpSurface);
  new G4LogicalSkinSurface("skinBot",botLV,OpSurface);


  auto mkReflector = [&](const char* name){
      auto s = new G4OpticalSurface(name);
      s->SetType(dielectric_metal);
      s->SetModel(unified);
      s->SetFinish(polished);
      auto mpt = new G4MaterialPropertiesTable();
      // Example: flat 95% reflectivity in visible; add your spectrum
      G4double e[4] = {1.5*eV, 4.0*eV,7.0*eV,10.0*eV};
      G4double R[4] = {1.0, 1.0,1.0,1.0};
      mpt->AddProperty("REFLECTIVITY", e, R, 2);
      s->SetMaterialPropertiesTable(mpt);
      return s;
  };
  auto reflTop = mkReflector("ReflTop");
  auto reflBot = mkReflector("ReflBot");

  // TPB Volumes
  const G4double tTPB    = 1.5*um;     // ~1–3 µm film is typical
  const G4double inset = 1.0*mm; // instead of 5 um

  // mother: LAr logical volume is LArLV (you already have it)
  
  // Side coating: a thin cylindrical shell hugging the wall
  auto tpbSideSolid = new G4Tubs("TPB_Side", r - tTPB - inset, r - inset, 0.5*innerHeight - tTPB - inset, 0.*deg, 360.*deg);
  auto tpbSideLV    = new G4LogicalVolume(tpbSideSolid, tpbMat, "TPB_Side_LV");
  auto tpbSidePV = new G4PVPlacement(nullptr, G4ThreeVector(), tpbSideLV, "TPB_Side_PV", innerCellLogic, true, 0);
  
  // Bottom disk
  auto tpbBotSolid = new G4Tubs("TPB_Bot", 0., r, 0.5*tTPB, 0.*deg, 360.*deg);
  auto tpbBotLV    = new G4LogicalVolume(tpbBotSolid, tpbMat, "TPB_Bot_LV");
  auto tpbBotPV = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-0.5*innerHeight + 0.5*tTPB + inset ), tpbBotLV,"TPB_Bot_PV", innerCellLogic, true, 0);
  
  // Top disk
  auto tpbTopSolid = new G4Tubs("TPB_Top", 0., r, 0.5*tTPB, 0.*deg, 360.*deg);
  auto tpbTopLV    = new G4LogicalVolume(tpbTopSolid, tpbMat, "TPB_Top_LV");
  auto tpbTopPV = new G4PVPlacement(nullptr, G4ThreeVector(0,0, +0.5*innerHeight - 0.5*tTPB - inset), tpbTopLV,"TPB_Top_PV", innerCellLogic, true, 0);
 
  // BUILDING PMTs + SiPMs configurations
  BuildPMTPatches();
  BuildSiPMPatches();

  ApplySensorConfig();
 
  //auto tpbUnion1 = new G4UnionSolid("TPB_union", tpbSideSolid, tpbTopSolid, nullptr, G4ThreeVector(0,0, +0.5*innerHeight - 0.5*tTPB - inset));
  //auto tpbUnion  = new G4UnionSolid("TPB_union2", tpbUnion1, tpbBotSolid, nullptr, G4ThreeVector(0,0,-0.5*innerHeight + 0.5*tTPB + inset ));

  //auto tpbLV = new G4LogicalVolume(tpbUnion, tpbMat, "TPB_LV");
  //auto tpbPV = new G4PVPlacement(nullptr, {}, tpbLV, "TPB_PV", innerCellLogic, true, 0);

  //new G4LogicalBorderSurface("LAr->Top",  tpbPV, topPV, reflTop);
  //new G4LogicalBorderSurface("LAr->Bot",  tpbPV, botPV, reflBot);

  // placing liquid scintillators at different angles to detect neutrons
  G4double angleArr[5] = {0.4363,0.6981,0.8727,1.0472,1.5708};
  G4double distances[5] = {1.0*m,1.0*m,1.0*m,1.0*m,1.0*m};
  G4ThreeVector vectArr[5];
  for(int i = 0;i<=4;i++){
      vectArr[i] = G4ThreeVector(distances[i],0.,0.);
  }   

  for(int i = 0;i<=4;i++){
      std::stringstream ss;
      ss << "A" << i;
      std::string fibName = ss.str();
      std::cout << "fiber diameter = " << fiberDiameter << " cm" << std::endl;
      placeLiquidScintillator(vectArr[i],angleArr[i],fiberDiameter,fiberLength,fiberMat,airLogic,fibName);
  } 


 auto visAttributes = new G4VisAttributes(G4Colour(1.,1.,1.)); // world is white
 //visAttributes->SetVisibility(false);
 //visAttributes->SetVisibility(true);
 //visAttributes->SetForceSolid(true);
 fLogicWorld->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);                                     

 visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,1.6));   // make inner LAr detector dark grey
 //visAttributes->SetVisibility(false);
 visAttributes->SetForceSolid(true);
 innerCellLogic->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,0.0));   // make Aluminium Layer yellow
 //visAttributes->SetVisibility(false);
 visAttributes->SetForceSolid(true);
 innerLayerLogic->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.4,0.4,0.4,0.9));   // make outer LAr detector dark grey and opaque
 visAttributes->SetForceSolid(true);
 outerCellLogic->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.9,0.1,0.9,0.6));  // make outer layer one some colour
 visAttributes->SetForceSolid(true);
 outerLayerOneLogic->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));  // make outer layer two some colour
 visAttributes->SetForceSolid(true);
 outerLayerTwoLogic->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));
 //visAttributes->SetVisibility(false);
 topLV->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));
 //visAttributes->SetVisibility(false);
 botLV->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));
 //visAttributes->SetVisibility(false);
 //tpbBotLV->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));
 //visAttributes->SetVisibility(false);
 //tpbTopLV->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.7,0.9,0.8,0.4));
 //visAttributes->SetVisibility(false);
 //tpbSideLV->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));   // make liquid scintilator turquoise
 //new_lvol_fiber->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));
 //testLogic->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 // PrintGeomParameters();         
  
  // ===================== Region-based production cuts =====================

// Inner LAr: keep dE/dx → scintillation; suppress soft EM secondaries.
{
 auto innerReg  = new G4Region("InnerLArRegion");
 innerReg->AddRootLogicalVolume(innerCellLogic);

// auto innerCuts = new G4ProductionCuts();
// innerCuts->SetProductionCut(5.0*CLHEP::mm, G4ProductionCuts::GetIndex("gamma"));
// innerCuts->SetProductionCut(5.0*CLHEP::mm, G4ProductionCuts::GetIndex("e-"));
// innerCuts->SetProductionCut(5.0*CLHEP::mm, G4ProductionCuts::GetIndex("e+"));
// innerCuts->SetProductionCut(5.0*CLHEP::mm, G4ProductionCuts::GetIndex("proton"));
// innerReg->SetProductionCuts(innerCuts);
}

// Outer LAr: very aggressive cuts (no EM analysis there).
{
 auto outerLArReg  = new G4Region("OuterLArRegion");
 outerLArReg->AddRootLogicalVolume(outerCellLogic);

// auto outerLArCuts = new G4ProductionCuts();
// outerLArCuts->SetProductionCut(5.0*CLHEP::cm, G4ProductionCuts::GetIndex("gamma"));
// outerLArCuts->SetProductionCut(5.0*CLHEP::cm, G4ProductionCuts::GetIndex("e-"));
// outerLArCuts->SetProductionCut(5.0*CLHEP::cm, G4ProductionCuts::GetIndex("e+"));
// outerLArCuts->SetProductionCut(5.0*CLHEP::cm, G4ProductionCuts::GetIndex("proton"));
// outerLArReg->SetProductionCuts(outerLArCuts);
}

// Air/world: huge cuts (don’t waste time on EM in air).
{
 auto airReg  = new G4Region("AirRegion");
 airReg->AddRootLogicalVolume(airLogic);

// auto airCuts = new G4ProductionCuts();
// airCuts->SetProductionCut(1.0*CLHEP::mm, G4ProductionCuts::GetIndex("gamma"));
// airCuts->SetProductionCut(1.0*CLHEP::mm, G4ProductionCuts::GetIndex("e-"));
// airCuts->SetProductionCut(1.0*CLHEP::mm, G4ProductionCuts::GetIndex("e+"));
// airCuts->SetProductionCut(0.01*CLHEP::mm, G4ProductionCuts::GetIndex("proton"));
// airReg->SetProductionCuts(airCuts);
}

//Steel shells: big cuts (no EM detail needed there).
{
 auto steelReg  = new G4Region("SteelRegion");
 steelReg->AddRootLogicalVolume(outerLayerOneLogic);
 steelReg->AddRootLogicalVolume(outerLayerTwoLogic);

// auto steelCuts = new G4ProductionCuts();
// steelCuts->SetProductionCut(1.0*CLHEP::cm, G4ProductionCuts::GetIndex("gamma"));
// steelCuts->SetProductionCut(1.0*CLHEP::cm, G4ProductionCuts::GetIndex("e-"));
// steelCuts->SetProductionCut(1.0*CLHEP::cm, G4ProductionCuts::GetIndex("e+"));
// steelCuts->SetProductionCut(1.0*CLHEP::cm, G4ProductionCuts::GetIndex("proton"));
// steelReg->SetProductionCuts(steelCuts);
}
// =======================================================================

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintGeomParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << fAbsorberMaterial << G4endl;
    
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeX,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is " 
         << G4BestUnit(fWorldSizeYZ,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(fAbsorberThickness,"Length")
         << " of " << fAbsorberMaterial->GetName();
  G4cout << ". The transverse size (YZ) is " 
         << G4BestUnit(fAbsorberSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the absorber "
         << G4BestUnit(fXposAbs,"Length");
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateConfigurableOpticalProperties()
{
  auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  if (lAr) {
    auto* mpt = lAr->GetMaterialPropertiesTable();
    if (mpt) {
      std::vector<G4double> energy = {9.5*eV, 9.69*eV, 9.89*eV, 10.0*eV};
      absorptionLength = {fLArAbsLength, fLArAbsLength, fLArAbsLength, fLArAbsLength};
      mpt->AddConstProperty("SCINTILLATIONYIELD", fLArScintYieldScale * fLArScintYieldPerMeV / MeV, false);
      mpt->AddProperty("ABSLENGTH", energy, absorptionLength, false);
      mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", fLArFastTime, false);
      mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", fLArSlowTime, false);
      mpt->AddConstProperty("SCINTILLATIONYIELD1", fLArFastFraction, false);
      mpt->AddConstProperty("SCINTILLATIONYIELD2", 1.0 - fLArFastFraction, false);
    }
    lAr->GetIonisation()->SetBirksConstant(birks);
  }

  if (tpbMat) {
    auto* mptTPB = tpbMat->GetMaterialPropertiesTable();
    if (mptTPB) {
      mptTPB->AddConstProperty("WLSMEANNUMBERPHOTONS", fTPBEfficiency, false);
    }
  }

  G4cout << "[Config] LAr yield=" << fLArScintYieldPerMeV
         << ", scale=" << fLArScintYieldScale
         << ", effective=" << fLArScintYieldPerMeV * fLArScintYieldScale
         << " photons/MeV, absLen=" << fLArAbsLength/cm
         << " cm, fastFrac=" << fLArFastFraction
         << ", TPB=" << fTPBEfficiency << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  std::cout << "argon40 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,40) << std::endl;
  std::cout << "argon36 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,36) << std::endl;
  std::cout << "argon38 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,38) << std::endl;
  std::cout << "argon39 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,39) << std::endl;
  std::cout << "argon41 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,41) << std::endl;

  if (pttoMaterial && fAbsorberMaterial != pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;

    if(innerCellLogic) { innerCellLogic->SetMaterial(fAbsorberMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial; 
    if(fLogicWorld) { fLogicWorld->SetMaterial(fWorldMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  fAbsorberThickness = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  fAbsorberSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeX(G4double val)
{
  fWorldSizeX = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  fWorldSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberXpos(G4double val)
{
  fXposAbs = val;
  ComputeGeomParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
  if ( fFieldMessenger.Get() == 0 ) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg =
      new G4GlobalMagFieldMessenger(fieldValue);
    //msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put( msg );        
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ChangeGeometry()
{
  fSolidWorld->SetXHalfLength(fWorldSizeX*0.5);
  fSolidWorld->SetYHalfLength(fWorldSizeYZ*0.5);
  fSolidWorld->SetZHalfLength(fWorldSizeYZ*0.5);

 // fSolidAbsorber->SetXHalfLength(fAbsorberThickness*0.5);
  //fSolidAbsorber->SetYHalfLength(fAbsorberSizeYZ*0.5);
  //fSolidAbsorber->SetZHalfLength(fAbsorberSizeYZ*0.5);

  svol_fiber->SetOuterRadius(fiberDiameter*0.5);
  svol_fiber->SetZHalfLength(fiberLength*0.5);
  svol_fiber->SetDeltaPhiAngle(2*M_PI);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildPMTPatches()
{
  // --- geometry numbers ---
  const G4double tile = fPMTSide;
  const G4double half = 0.5*tile;
  const G4double thick = 0.1*mm;  // thin patch
  const G4double eps   = 0.01*mm;  // inset into LAr so it's inside the cell

  const G4double R = 0.5*innerDiameter;
  const G4double H = innerHeight;

  // sanity: tiles must fit
  if (tile > 2*R) {
    G4cout << "[PMT] ERROR: PMT tile bigger than diameter.\n";
    return;
  }

  const G4double zTop = +0.5*H - 0.5*thick - eps;
  const G4double zBot = -0.5*H + 0.5*thick + eps;

  // gap between tiles (optional)
  const G4double gap = 0.0*mm;
  const G4double dx = half + 0.5*gap;
  const G4double dy = half + 0.5*gap;

  // 2x2 centers
  std::vector<G4ThreeVector> centers = {
    {+dx,+dy,0}, {+dx,-dy,0}, {-dx,+dy,0}, {-dx,-dy,0}
  };

  // dummy material; optical behavior is defined by the surface
  auto dummyMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  auto solid = new G4Box("PMTtileSolid", half, half, 0.5*thick);
  auto lv    = new G4LogicalVolume(solid, dummyMat, "PMTtileLV");

  fTopPMT_PV.clear(); fBotPMT_PV.clear();
  fTopPMT_Surf.clear(); fBotPMT_Surf.clear();

  for (int i=0;i<4;i++) {
    auto pvT = new G4PVPlacement(
      nullptr,
      centers[i] + G4ThreeVector(0,0,zTop),
      lv,
      "PMTtileTopPV",
      innerCellLogic,
      false,
      i,
      false
    );
    auto pvB = new G4PVPlacement(
      nullptr,
      centers[i] + G4ThreeVector(0,0,zBot),
      lv,
      "PMTtileBotPV",
      innerCellLogic,
      false,
      i,
      false
    );

    fTopPMT_PV.push_back(pvT);
    fBotPMT_PV.push_back(pvB);

    auto sTop = new G4OpticalSurface(("PMTsurfTop_"+std::to_string(i)).c_str());
    auto sBot = new G4OpticalSurface(("PMTsurfBot_"+std::to_string(i)).c_str());

    for (auto s : {sTop, sBot}) {
      s->SetType(dielectric_metal);
      s->SetModel(unified);
      s->SetFinish(polished);
 
      auto mptESR = MakeESRMPT();
      s->SetMaterialPropertiesTable(mptESR);
    }

    fTopPMT_Surf.push_back(sTop);
    fBotPMT_Surf.push_back(sBot);

    // Attach border surfaces (both directions)
    new G4LogicalBorderSurface(
      ("LAr_to_PMTtop_"+std::to_string(i)).c_str(),
      innerCellPhysi,
      pvT,
      sTop
    );
    new G4LogicalBorderSurface(
      ("PMTtop_to_LAr_"+std::to_string(i)).c_str(),
      pvT,
      innerCellPhysi,
      sTop
    );

    new G4LogicalBorderSurface(
      ("LAr_to_PMTbot_"+std::to_string(i)).c_str(),
      innerCellPhysi,
      pvB,
      sBot
    );
    new G4LogicalBorderSurface(
      ("PMTbot_to_LAr_"+std::to_string(i)).c_str(),
      pvB,
      innerCellPhysi,
      sBot
    );
  }

  G4cout << "[PMT] Built 4 top + 4 bottom PMT tiles.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildSiPMPatches()
{

  auto dummyMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  const G4double tile = fSiPMTile;
  const G4double thick = 0.1*mm; // radial thickness
  const G4double eps   = 0.01*mm; // inset into LAr

  const G4double R = 0.5*innerDiameter;
  const G4double H = innerHeight;

  const G4ThreeVector zHat(0,0,1);

  // Max rows that fit
  const int maxRows = static_cast<int>(std::floor(H / tile));
  // Tiles around circumference
  const G4double circumference = 2.0 * CLHEP::pi * R;
  const int nPhi = static_cast<int>(std::floor(circumference / tile));
  fSiPM_TilesPerRow = nPhi;

  G4cout << "[SiPM] maxRows=" << maxRows
         << " tilesPerRow=" << nPhi
         << " (H=" << H/cm << " cm, R=" << R/cm << " cm)\n";


  // Local box: X=tangential, Y=vertical (z), Z=radial thickness
  auto solid = new G4Box("SiPMtileSolid", 0.5*tile, 0.5*tile, 0.5*thick);
  auto lv    = new G4LogicalVolume(solid, dummyMat, "SiPMtileLV");
  
  auto visAttributes = new G4VisAttributes(G4Colour(1.,1.,1.)); // world is white
  visAttributes->SetVisibility(true);
  visAttributes->SetForceSolid(true);
  lv->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  fSiPM_PV.clear();
  fSiPM_Surf.clear();
  fSiPM_RowOfTile.clear();

  const G4double dphi = 2.0 * CLHEP::pi / nPhi;

  for (int iRow=0; iRow<maxRows; ++iRow) {

    const G4double z = -0.5*H + (iRow + 0.5)*tile;

    for (int j=0; j<nPhi; ++j) {

      const G4double phi = j * dphi;

      const G4double x = (R - eps - 0.5*thick) * std::cos(-phi);
      const G4double y = (R - eps - 0.5*thick) * std::sin(-phi);

      G4ThreeVector pos = G4ThreeVector(x,y,z);
      
    // Build rotation: columns = local axes expressed in global coordinates
    auto rot = new G4RotationMatrix();
    rot->rotateZ(phi);
    rot->rotateY(90.0*deg);

      const int copyNo = iRow*nPhi + j;

      auto pv = new G4PVPlacement(
        rot,
        pos,
        lv,
        "SiPMtilePV",
        innerCellLogic,
        false,
        copyNo,
        false
      );


      fSiPM_PV.push_back(pv);
      fSiPM_RowOfTile.push_back(iRow);

      auto s = new G4OpticalSurface(("SiPMsurf_"+std::to_string(copyNo)).c_str());
      s->SetType(dielectric_metal);
      s->SetModel(unified);
      s->SetFinish(polished);
      
      auto mptESR = MakeESRMPT();
      s->SetMaterialPropertiesTable(mptESR);
      fSiPM_Surf.push_back(s);

      new G4LogicalBorderSurface(
        ("LAr_to_SiPM_"+std::to_string(copyNo)).c_str(),
        innerCellPhysi,
        pv,
        s
      );
      new G4LogicalBorderSurface(
        ("SiPM_to_LAr_"+std::to_string(copyNo)).c_str(),
        pv,
        innerCellPhysi,
        s
      );
    }
  }

  G4cout << "[SiPM] Built tiles: " << fSiPM_PV.size()
         << " (rows=" << maxRows << ", perRow=" << nPhi << ")\n";
}

// ooOOOOooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooo

void DetectorConstruction::ApplySensorConfig()
{
  // ---------------------------
  // Clamp configuration knobs
  // ---------------------------
  if (fTopPMTs < 0) fTopPMTs = 0;
  if (fTopPMTs > 4) fTopPMTs = 4;

  if (fBotPMTs < 0) fBotPMTs = 0;
  if (fBotPMTs > 4) fBotPMTs = 4;

  if (fSiPMRows < 0) fSiPMRows = 0;
  if (fSiPMRows > 16) fSiPMRows = 16;

  // ---------------------------
  // Build / rebuild cached MPTs
  // ---------------------------

  // ESR MPT: build once (inactive behavior)
  if (!fMPT_ESR) {
    fMPT_ESR = MakeFlatMPT(0.0, fReflectorReflectivity); // inactive reflector: no detection, configurable reflectivity
  }

  // Top PMT active MPT: rebuild if PDE/top response mode/curve file or PMT reflectivity changed
  if (!fMPT_TopPMT_Active) 
  {
    if (fPMTResponseMode == SensorResponseMode::Flat) {
      fMPT_TopPMT_Active = MakeFlatMPT(fPDE_TopPMT, fRefl_PMT);
    } else if (fPMTResponseMode == SensorResponseMode::CsvCurve) {
      // Scale the CSV curve so its peak equals fPDE_TopPMT.
      fMPT_TopPMT_Active = MakeCurveMPTFromFile(fPMTCurveFile, fRefl_PMT, fPDE_TopPMT);
    } else if (fPMTResponseMode == SensorResponseMode::CsvCurveRaw) {
      // Use the CSV curve as an absolute QE/PDE curve; do not scale to fPDE_TopPMT.
      fMPT_TopPMT_Active = MakeCurveMPTFromFile(fPMTCurveFile, fRefl_PMT, 0.0);
    } else if (fPMTResponseMode == SensorResponseMode::BuiltinCurveRaw) {
      // Use the built-in PMT curve exactly as stored in the code; do not scale to fPDE_TopPMT.
      fMPT_TopPMT_Active = MakePMTMPT(fRefl_PMT, 0.0);
    } else {
      // Default: scale the built-in PMT curve so its peak equals fPDE_TopPMT.
      fMPT_TopPMT_Active = MakePMTMPT(fRefl_PMT, fPDE_TopPMT);
    }
  }

  // Bottom PMT active MPT: rebuild if PDE/bottom response mode/curve file or PMT reflectivity changed
  if (!fMPT_BotPMT_Active) 
  {
    if (fPMTResponseMode == SensorResponseMode::Flat) {
      fMPT_BotPMT_Active = MakeFlatMPT(fPDE_BotPMT, fRefl_PMT);
    } else if (fPMTResponseMode == SensorResponseMode::CsvCurve) {
      // Scale the CSV curve so its peak equals fPDE_BotPMT.
      fMPT_BotPMT_Active = MakeCurveMPTFromFile(fPMTCurveFile, fRefl_PMT, fPDE_BotPMT);
    } else if (fPMTResponseMode == SensorResponseMode::CsvCurveRaw) {
      // Use the CSV curve as an absolute QE/PDE curve; do not scale to fPDE_BotPMT.
      fMPT_BotPMT_Active = MakeCurveMPTFromFile(fPMTCurveFile, fRefl_PMT, 0.0);
    } else if (fPMTResponseMode == SensorResponseMode::BuiltinCurveRaw) {
      // Use the built-in PMT curve exactly as stored in the code; do not scale to fPDE_BotPMT.
      fMPT_BotPMT_Active = MakePMTMPT(fRefl_PMT, 0.0);
    } else {
      // Default: scale the built-in PMT curve so its peak equals fPDE_BotPMT.
      fMPT_BotPMT_Active = MakePMTMPT(fRefl_PMT, fPDE_BotPMT);
    }
  }

  // SiPM active MPT: rebuild if PDE/response mode/curve file or SiPM reflectivity changed
  if (!fMPT_SiPM_Active) 
  {
    if (fSiPMResponseMode == SensorResponseMode::Flat) {
      fMPT_SiPM_Active = MakeFlatMPT(fPDE_SiPM, fRefl_SiPM);
    } else if (fSiPMResponseMode == SensorResponseMode::CsvCurve) {
      // Scale the CSV curve so its peak equals fPDE_SiPM.
      fMPT_SiPM_Active = MakeCurveMPTFromFile(fSiPMCurveFile, fRefl_SiPM, fPDE_SiPM);
    } else if (fSiPMResponseMode == SensorResponseMode::CsvCurveRaw) {
      // Use the CSV curve as an absolute QE/PDE curve; do not scale to fPDE_SiPM.
      fMPT_SiPM_Active = MakeCurveMPTFromFile(fSiPMCurveFile, fRefl_SiPM, 0.0);
    } else if (fSiPMResponseMode == SensorResponseMode::BuiltinCurveRaw) {
      // Use the built-in SiPM curve exactly as stored in the code; do not scale to fPDE_SiPM.
      fMPT_SiPM_Active = MakeSiPMMPT(fRefl_SiPM, 0.0);
    } else {
      // Default: scale the built-in SiPM curve so its peak equals fPDE_SiPM.
      fMPT_SiPM_Active = MakeSiPMMPT(fRefl_SiPM, fPDE_SiPM);
    }
  }

  // ---------------------------
  // Apply to PMT surfaces
  // ---------------------------
  for (int i = 0; i < 4; ++i) {
    const bool topOn = (i < fTopPMTs);
    const bool botOn = (i < fBotPMTs);

    if (fTopPMT_Surf[i]) {
      fTopPMT_Surf[i]->SetMaterialPropertiesTable(topOn ? fMPT_TopPMT_Active : fMPT_ESR);
    }
    if (fBotPMT_Surf[i]) {
      fBotPMT_Surf[i]->SetMaterialPropertiesTable(botOn ? fMPT_BotPMT_Active : fMPT_ESR);
    }
  }

  // ---------------------------
  // Apply to SiPM surfaces
  // ---------------------------
  for (size_t k = 0; k < fSiPM_Surf.size(); ++k) {
    const int row = fSiPM_RowOfTile[k];
    const bool on = (row < fSiPMRows);

    if (fSiPM_Surf[k]) {
      fSiPM_Surf[k]->SetMaterialPropertiesTable(on ? fMPT_SiPM_Active : fMPT_ESR);
    }
  }

  std::cout << "number of active top PMTs: " << fTopPMTs << "\n number of active bottom PMTs: " << fBotPMTs << "\n number of active SiPM rows: " << fSiPMRows << std::endl;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

// oooOOOOOOOOOooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
