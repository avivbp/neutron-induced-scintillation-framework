#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include <sstream>

DetectorMessenger::DetectorMessenger(DetectorConstruction* det, Run* run)
    : G4UImessenger(), fDetector(det), fRun(run) {


  fDetDir = new G4UIdirectory("/det/");
  fDetDir->SetGuidance("Detector configuration");

  fDetectorModelCmd = new G4UIcmdWithAString("/det/setDetectorModel", this);
  fDetectorModelCmd->SetGuidance("Set detector model: nested_cell or box_cryostat.");
  fDetectorModelCmd->SetParameterName("DetectorModel", false);
  fDetectorModelCmd->AvailableForStates(G4State_PreInit);

  fBoxDimensionsCmd = new G4UIcmdWithAString("/det/setBoxDimensionsCm", this);
  fBoxDimensionsCmd->SetGuidance("Set full active-LAr box dimensions in cm: x y z.");
  fBoxDimensionsCmd->SetParameterName("BoxDimensionsCm", false);
  fBoxDimensionsCmd->AvailableForStates(G4State_PreInit);

  fBoxFiducialMarginCmd = new G4UIcmdWithAString("/det/setBoxFiducialMarginCm", this);
  fBoxFiducialMarginCmd->SetGuidance("Set inset fiducial margins in cm: x y z.");
  fBoxFiducialMarginCmd->SetParameterName("BoxFiducialMarginCm", false);
  fBoxFiducialMarginCmd->AvailableForStates(G4State_PreInit);

  fBoxCryostatThicknessCmd = new G4UIcmdWithADouble(
      "/det/setBoxCryostatThicknessCm", this);
  fBoxCryostatThicknessCmd->SetGuidance("Set box cryostat shell thickness in cm.");
  fBoxCryostatThicknessCmd->SetParameterName("BoxCryostatThicknessCm", false);
  fBoxCryostatThicknessCmd->SetRange("BoxCryostatThicknessCm>0.");
  fBoxCryostatThicknessCmd->AvailableForStates(G4State_PreInit);

  fBoxSensorLayoutsCmd = new G4UIcmdWithAString("/det/setBoxSensorLayouts", this);
  fBoxSensorLayoutsCmd->SetGuidance(
      "Set semicolon-separated box sensor layouts generated from YAML.");
  fBoxSensorLayoutsCmd->SetParameterName("BoxSensorLayouts", false);
  fBoxSensorLayoutsCmd->AvailableForStates(G4State_PreInit);

  // Geometry and optical-material commands are intended to be used before /run/initialize.
  fInnerRadiusCmd = new G4UIcmdWithADouble("/det/setInnerRadiusCm", this);
  fInnerRadiusCmd->SetGuidance("Set inner LAr radius in cm.");
  fInnerRadiusCmd->SetParameterName("InnerRadiusCm", false);
  fInnerRadiusCmd->SetRange("InnerRadiusCm>0.");
  fInnerRadiusCmd->AvailableForStates(G4State_PreInit);

  fInnerHeightCmd = new G4UIcmdWithADouble("/det/setInnerHeightCm", this);
  fInnerHeightCmd->SetGuidance("Set inner LAr height in cm.");
  fInnerHeightCmd->SetParameterName("InnerHeightCm", false);
  fInnerHeightCmd->SetRange("InnerHeightCm>0.");
  fInnerHeightCmd->AvailableForStates(G4State_PreInit);

  fOuterDiameterCmd = new G4UIcmdWithADouble("/det/setOuterDiameterCm", this);
  fOuterDiameterCmd->SetGuidance("Set outer LAr diameter in cm.");
  fOuterDiameterCmd->SetParameterName("OuterDiameterCm", false);
  fOuterDiameterCmd->SetRange("OuterDiameterCm>0.");
  fOuterDiameterCmd->AvailableForStates(G4State_PreInit);

  fOuterHeightCmd = new G4UIcmdWithADouble("/det/setOuterHeightCm", this);
  fOuterHeightCmd->SetGuidance("Set outer LAr height in cm.");
  fOuterHeightCmd->SetParameterName("OuterHeightCm", false);
  fOuterHeightCmd->SetRange("OuterHeightCm>0.");
  fOuterHeightCmd->AvailableForStates(G4State_PreInit);

  fOuterScintillationCmd = new G4UIcmdWithABool("/det/setOuterScintillation", this);
  fOuterScintillationCmd->SetGuidance(
      "Enable scintillation photon production in the outer LAr volume.");
  fOuterScintillationCmd->SetParameterName("OuterScintillation", false);
  fOuterScintillationCmd->AvailableForStates(G4State_PreInit);

  fNeutronDetectorsCmd = new G4UIcmdWithAString("/det/setNeutronDetectors", this);
  fNeutronDetectorsCmd->SetGuidance("Set neutron liquid scintillators as label:angle_deg:distance_cm;...");
  fNeutronDetectorsCmd->SetParameterName("NeutronDetectors", false);
  fNeutronDetectorsCmd->AvailableForStates(G4State_PreInit);

  fTOFWindowCmd = new G4UIcmdWithAString("/det/setTOFWindowNs", this);
  fTOFWindowCmd->SetGuidance("Set event TOF selection window in ns as: min max");
  fTOFWindowCmd->SetParameterName("TOFWindowNs", false);
  fTOFWindowCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPrimaryEnergyCmd = new G4UIcmdWithADouble("/det/setPrimaryEnergyMeV", this);
  fPrimaryEnergyCmd->SetGuidance("Set primary neutron kinetic energy in MeV.");
  fPrimaryEnergyCmd->SetParameterName("PrimaryEnergyMeV", false);
  fPrimaryEnergyCmd->SetRange("PrimaryEnergyMeV>0.");
  fPrimaryEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPMTSideCmd = new G4UIcmdWithADouble("/det/setPMTSideCm", this);
  fPMTSideCmd->SetGuidance("Set square PMT side length in cm.");
  fPMTSideCmd->SetParameterName("PMTSideCm", false);
  fPMTSideCmd->SetRange("PMTSideCm>0.");
  fPMTSideCmd->AvailableForStates(G4State_PreInit);

  fSiPMTileCmd = new G4UIcmdWithADouble("/det/setSiPMTileCm", this);
  fSiPMTileCmd->SetGuidance("Set square SiPM tile side length in cm.");
  fSiPMTileCmd->SetParameterName("SiPMTileCm", false);
  fSiPMTileCmd->SetRange("SiPMTileCm>0.");
  fSiPMTileCmd->AvailableForStates(G4State_PreInit);

  fAbsLengthCmCmd = new G4UIcmdWithADouble("/det/setAbsLengthCm", this);
  fAbsLengthCmCmd->SetGuidance("Set LAr absorption length in cm.");
  fAbsLengthCmCmd->SetParameterName("AbsLengthCm", false);
  fAbsLengthCmCmd->SetRange("AbsLengthCm>0.");
  fAbsLengthCmCmd->AvailableForStates(G4State_PreInit);

  fBirksMmPerMeVCmd = new G4UIcmdWithADouble("/det/setBirksMmPerMeV", this);
  fBirksMmPerMeVCmd->SetGuidance("Set Birks constant in mm/MeV.");
  fBirksMmPerMeVCmd->SetParameterName("BirksMmPerMeV", false);
  fBirksMmPerMeVCmd->SetRange("BirksMmPerMeV>=0.");
  fBirksMmPerMeVCmd->AvailableForStates(G4State_PreInit);

  fScintYieldCmd = new G4UIcmdWithADouble("/det/setScintYieldPerMeV", this);
  fScintYieldCmd->SetGuidance("Set LAr scintillation yield in photons/MeV before scaling.");
  fScintYieldCmd->SetParameterName("ScintYieldPerMeV", false);
  fScintYieldCmd->SetRange("ScintYieldPerMeV>=0.");
  fScintYieldCmd->AvailableForStates(G4State_PreInit);

  fScintScaleCmd = new G4UIcmdWithADouble("/det/setScintYieldScale", this);
  fScintScaleCmd->SetGuidance("Set multiplicative scintillation yield scale.");
  fScintScaleCmd->SetParameterName("ScintYieldScale", false);
  fScintScaleCmd->SetRange("ScintYieldScale>=0.");
  fScintScaleCmd->AvailableForStates(G4State_PreInit);

  fIonScintScaleCmd = new G4UIcmdWithADouble("/det/setIonScintYieldScale", this);
  fIonScintScaleCmd->SetGuidance(
      "Set the LAr ion/electron scintillation yield ratio after Birks quenching.");
  fIonScintScaleCmd->SetParameterName("IonScintYieldScale", false);
  fIonScintScaleCmd->SetRange("IonScintYieldScale>=0. && IonScintYieldScale<=1.");
  fIonScintScaleCmd->AvailableForStates(G4State_PreInit);

  fFastTimeCmd = new G4UIcmdWithADouble("/det/setFastTimeNs", this);
  fFastTimeCmd->SetGuidance("Set fast scintillation time constant in ns.");
  fFastTimeCmd->SetParameterName("FastTimeNs", false);
  fFastTimeCmd->SetRange("FastTimeNs>0.");
  fFastTimeCmd->AvailableForStates(G4State_PreInit);

  fSlowTimeCmd = new G4UIcmdWithADouble("/det/setSlowTimeNs", this);
  fSlowTimeCmd->SetGuidance("Set slow scintillation time constant in ns.");
  fSlowTimeCmd->SetParameterName("SlowTimeNs", false);
  fSlowTimeCmd->SetRange("SlowTimeNs>0.");
  fSlowTimeCmd->AvailableForStates(G4State_PreInit);

  fFastFractionCmd = new G4UIcmdWithADouble("/det/setFastFraction", this);
  fFastFractionCmd->SetGuidance("Set fast scintillation fraction.");
  fFastFractionCmd->SetParameterName("FastFraction", false);
  fFastFractionCmd->SetRange("FastFraction>=0. && FastFraction<=1.");
  fFastFractionCmd->AvailableForStates(G4State_PreInit);

  fTPBEfficiencyCmd = new G4UIcmdWithADouble("/det/setTPBEfficiency", this);
  fTPBEfficiencyCmd->SetGuidance("Set TPB WLS mean number of photons per absorbed photon.");
  fTPBEfficiencyCmd->SetParameterName("TPBEfficiency", false);
  fTPBEfficiencyCmd->SetRange("TPBEfficiency>=0. && TPBEfficiency<=1.");
  fTPBEfficiencyCmd->AvailableForStates(G4State_PreInit);

  fReflectorReflCmd = new G4UIcmdWithADouble("/det/setReflectorReflectivity", this);
  fReflectorReflCmd->SetGuidance("Set inactive reflector reflectivity.");
  fReflectorReflCmd->SetParameterName("ReflectorReflectivity", false);
  fReflectorReflCmd->SetRange("ReflectorReflectivity>=0. && ReflectorReflectivity<=1.");
  fReflectorReflCmd->AvailableForStates(G4State_PreInit);

  fTopPMTCmd = new G4UIcmdWithAnInteger("/det/setTopPMTs", this);
  fTopPMTCmd->SetGuidance("Active PMTs on top (0..4)");
  fTopPMTCmd->SetParameterName("nTopPMT",false);
  fTopPMTCmd->SetRange("nTopPMT>=0 && nTopPMT<=4");
  fTopPMTCmd->AvailableForStates(G4State_Idle);

  fBotPMTCmd = new G4UIcmdWithAnInteger("/det/setBotPMTs", this);
  fBotPMTCmd->SetGuidance("Active PMTs on bottom (0..4)");
  fBotPMTCmd->SetParameterName("nBotPMT",false);
  fBotPMTCmd->SetRange("nBotPMT>=0 && nBotPMT<=4");
  fBotPMTCmd->AvailableForStates(G4State_Idle);

  fSiPMRowsCmd = new G4UIcmdWithAnInteger("/det/setSiPMRows", this);
  fSiPMRowsCmd->SetGuidance("Active SiPM rows (0..16)");
  fSiPMRowsCmd->SetParameterName("nSiPM",false);
  fSiPMRowsCmd->SetRange("nSiPM>=0 && nSiPM<=16");
  fSiPMRowsCmd->AvailableForStates(G4State_Idle);

  fPDETopPMTCmd = new G4UIcmdWithADouble("/det/setPDETopPMT", this);
  fPDETopPMTCmd->SetGuidance("Top PMT PDE (0..1)");
  fPDETopPMTCmd->SetParameterName("TopPMTPDE",false);
  fPDETopPMTCmd->SetRange("TopPMTPDE>=0.0 && TopPMTPDE<=1.0");
  fPDETopPMTCmd->AvailableForStates(G4State_Idle);

  fPDEBotPMTCmd = new G4UIcmdWithADouble("/det/setPDEBotPMT", this);
  fPDEBotPMTCmd->SetGuidance("Bot PMT PDE (0..1)");
  fPDEBotPMTCmd->SetParameterName("BotPMTPDE",false);
  fPDEBotPMTCmd->SetRange("BotPMTPDE>=0.0 && BotPMTPDE<=1.0");
  fPDEBotPMTCmd->AvailableForStates(G4State_Idle);

  fPDESiPMCmd = new G4UIcmdWithADouble("/det/setPDESiPM", this);
  fPDESiPMCmd->SetGuidance("SiPM PDE (0..1)");
  fPDESiPMCmd->SetParameterName("SiPMPDE",false);
  fPDESiPMCmd->SetRange("SiPMPDE>=0.0 && SiPMPDE<=1.0");
  fPDESiPMCmd->AvailableForStates(G4State_Idle);

  fReflESRCmd = new G4UIcmdWithADouble("/det/setReflESR", this);
  fReflESRCmd->SetGuidance("ESR reflectivity used when a tile is inactive (0..1)");
  fReflESRCmd->SetParameterName("RefleESR", false);
  fReflESRCmd->SetRange("RefleESR>=0.0 && RefleESR<=1.0");
  fReflESRCmd->AvailableForStates(G4State_Idle);

  fReflPMTCmd = new G4UIcmdWithADouble("/det/setReflPMT", this);
  fReflPMTCmd->SetGuidance("PMT patch reflectivity when active (0..1)");
  fReflPMTCmd->SetParameterName("PMTRefle",false);
  fReflPMTCmd->SetRange("PMTRefle>=0.0 && PMTRefle<=1.0");
  fReflPMTCmd->AvailableForStates(G4State_Idle);

  fReflSiPMCmd = new G4UIcmdWithADouble("/det/setReflSiPM", this);
  fReflSiPMCmd->SetGuidance("SiPM patch reflectivity when active (0..1)");
  fReflSiPMCmd->SetParameterName("SiPMRefle",false);
  fReflSiPMCmd->SetRange("SiPMRefle>=0.0 && SiPMRefle<=1.0");
  fReflSiPMCmd->AvailableForStates(G4State_Idle);

  fPMTResponseModeCmd = new G4UIcmdWithAString("/det/setPMTResponseMode", this);
  fPMTResponseModeCmd->SetGuidance("PMT response mode: flat, builtin_curve, or csv_curve.");
  fPMTResponseModeCmd->SetParameterName("PMTResponseMode", false);
  fPMTResponseModeCmd->AvailableForStates(G4State_Idle);

  fSiPMResponseModeCmd = new G4UIcmdWithAString("/det/setSiPMResponseMode", this);
  fSiPMResponseModeCmd->SetGuidance("SiPM response mode: flat, builtin_curve, or csv_curve.");
  fSiPMResponseModeCmd->SetParameterName("SiPMResponseMode", false);
  fSiPMResponseModeCmd->AvailableForStates(G4State_Idle);

  fPMTCurveFileCmd = new G4UIcmdWithAString("/det/setPMTCurveFile", this);
  fPMTCurveFileCmd->SetGuidance("CSV file for PMT response curve. Columns: wavelength_nm,qe or energy_eV,qe.");
  fPMTCurveFileCmd->SetParameterName("PMTCurveFile", false);
  fPMTCurveFileCmd->AvailableForStates(G4State_Idle);

  fSiPMCurveFileCmd = new G4UIcmdWithAString("/det/setSiPMCurveFile", this);
  fSiPMCurveFileCmd->SetGuidance("CSV file for SiPM response curve. Columns: wavelength_nm,qe or energy_eV,qe.");
  fSiPMCurveFileCmd->SetParameterName("SiPMCurveFile", false);
  fSiPMCurveFileCmd->AvailableForStates(G4State_Idle);

    fSimDir = new G4UIdirectory("/mysim/");
    fSimDir->SetGuidance("Commands to control simulation parameters");

    fBirksCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setBirks", this);
    fBirksCmd->SetGuidance("Set Birks constant");
    fBirksCmd->SetParameterName("birks", false);
    fBirksCmd->SetRange("birks > 0.0");
    fBirksCmd->SetUnitCategory("Length");
    fBirksCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setAbsLength", this);
    fAbsLengthCmd->SetGuidance("Set absorption length");
    fAbsLengthCmd->SetParameterName("absLength", false);
    fAbsLengthCmd->SetRange("absLength > 0.0");
    fAbsLengthCmd->SetUnitCategory("Length");
    fAbsLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fTpbEficCmd = new G4UIcmdWithADouble("/mysim/setTpbEfic", this);
    fTpbEficCmd->SetGuidance("Set TPB efficiency");
    fTpbEficCmd->SetParameterName("tpbEfic", false);
    fTpbEficCmd->SetRange("tpbEfic >= 0.0 && tpbEfic <= 1.0");
    fTpbEficCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fDetSizeCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setDetSize", this);
    fDetSizeCmd->SetGuidance("Set outer LAr cell diameter");
    fDetSizeCmd->SetParameterName("outerDiameter", false);
    fDetSizeCmd->SetRange("outerDiameter > 0.0");
    fDetSizeCmd->SetUnitCategory("Length");
    fDetSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

DetectorMessenger::~DetectorMessenger() {
  delete fBirksCmd;
  delete fAbsLengthCmd;
  delete fTpbEficCmd;
  delete fDetSizeCmd;
  delete fDetectorModelCmd;
  delete fBoxDimensionsCmd;
  delete fBoxFiducialMarginCmd;
  delete fBoxCryostatThicknessCmd;
  delete fBoxSensorLayoutsCmd;
  delete fInnerRadiusCmd;
  delete fInnerHeightCmd;
  delete fOuterDiameterCmd;
  delete fOuterHeightCmd;
  delete fOuterScintillationCmd;
  delete fNeutronDetectorsCmd;
  delete fTOFWindowCmd;
  delete fPrimaryEnergyCmd;
  delete fPMTSideCmd;
  delete fSiPMTileCmd;
  delete fAbsLengthCmCmd;
  delete fBirksMmPerMeVCmd;
  delete fScintYieldCmd;
  delete fScintScaleCmd;
  delete fIonScintScaleCmd;
  delete fFastTimeCmd;
  delete fSlowTimeCmd;
  delete fFastFractionCmd;
  delete fTPBEfficiencyCmd;
  delete fReflectorReflCmd;
  delete fSimDir;
  delete fTopPMTCmd;
  delete fBotPMTCmd;
  delete fSiPMRowsCmd;
  delete fPDETopPMTCmd;
  delete fPDEBotPMTCmd;
  delete fPDESiPMCmd;
  delete fReflESRCmd;
  delete fReflPMTCmd;
  delete fReflSiPMCmd;
  delete fPMTResponseModeCmd;
  delete fSiPMResponseModeCmd;
  delete fPMTCurveFileCmd;
  delete fSiPMCurveFileCmd;
  delete fDetDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fDetectorModelCmd) {
        fDetector->SetDetectorModel(newValue);
    } else if (command == fBoxDimensionsCmd) {
        fDetector->SetBoxDimensionsCm(newValue);
    } else if (command == fBoxFiducialMarginCmd) {
        fDetector->SetBoxFiducialMarginCm(newValue);
    } else if (command == fBoxCryostatThicknessCmd) {
        fDetector->SetBoxCryostatThicknessCm(
            fBoxCryostatThicknessCmd->GetNewDoubleValue(newValue));
    } else if (command == fBoxSensorLayoutsCmd) {
        fDetector->SetBoxSensorLayouts(newValue);
    } else if (command == fInnerRadiusCmd) {
        fDetector->SetInnerRadiusCm(fInnerRadiusCmd->GetNewDoubleValue(newValue));
    } else if (command == fInnerHeightCmd) {
        fDetector->SetInnerHeightCm(fInnerHeightCmd->GetNewDoubleValue(newValue));
    } else if (command == fOuterDiameterCmd) {
        fDetector->SetOuterDiameterCm(fOuterDiameterCmd->GetNewDoubleValue(newValue));
    } else if (command == fOuterHeightCmd) {
        fDetector->SetOuterHeightCm(fOuterHeightCmd->GetNewDoubleValue(newValue));
    } else if (command == fOuterScintillationCmd) {
        fDetector->SetOuterScintillation(
            fOuterScintillationCmd->GetNewBoolValue(newValue));
    } else if (command == fNeutronDetectorsCmd) {
        fDetector->SetNeutronDetectors(newValue);
    } else if (command == fTOFWindowCmd) {
        std::stringstream ss(newValue);
        G4double minNs = 0.0;
        G4double maxNs = 0.0;
        if (ss >> minNs >> maxNs) {
            fDetector->SetTOFWindowNs(minNs, maxNs);
        } else {
            G4cerr << "Invalid TOF window: " << newValue << G4endl;
        }
    } else if (command == fPrimaryEnergyCmd) {
        fDetector->SetPrimaryEnergyMeV(fPrimaryEnergyCmd->GetNewDoubleValue(newValue));
    } else if (command == fPMTSideCmd) {
        fDetector->SetPMTSideCm(fPMTSideCmd->GetNewDoubleValue(newValue));
    } else if (command == fSiPMTileCmd) {
        fDetector->SetSiPMTileCm(fSiPMTileCmd->GetNewDoubleValue(newValue));
    } else if (command == fAbsLengthCmCmd) {
        fDetector->SetAbsLengthCm(fAbsLengthCmCmd->GetNewDoubleValue(newValue));
    } else if (command == fBirksMmPerMeVCmd) {
        fDetector->SetBirksMmPerMeV(fBirksMmPerMeVCmd->GetNewDoubleValue(newValue));
    } else if (command == fScintYieldCmd) {
        fDetector->SetScintYieldPerMeV(fScintYieldCmd->GetNewDoubleValue(newValue));
    } else if (command == fScintScaleCmd) {
        fDetector->SetScintYieldScale(fScintScaleCmd->GetNewDoubleValue(newValue));
    } else if (command == fIonScintScaleCmd) {
        fDetector->SetIonScintYieldScale(
            fIonScintScaleCmd->GetNewDoubleValue(newValue));
    } else if (command == fFastTimeCmd) {
        fDetector->SetFastTimeNs(fFastTimeCmd->GetNewDoubleValue(newValue));
    } else if (command == fSlowTimeCmd) {
        fDetector->SetSlowTimeNs(fSlowTimeCmd->GetNewDoubleValue(newValue));
    } else if (command == fFastFractionCmd) {
        fDetector->SetFastFraction(fFastFractionCmd->GetNewDoubleValue(newValue));
    } else if (command == fTPBEfficiencyCmd) {
        fDetector->SetTPBEfficiency(fTPBEfficiencyCmd->GetNewDoubleValue(newValue));
    } else if (command == fReflectorReflCmd) {
        fDetector->SetReflectorReflectivity(fReflectorReflCmd->GetNewDoubleValue(newValue));
    } else 
    if (command == fBirksCmd) {
        fDetector->SetBirksConstant(fBirksCmd->GetNewDoubleValue(newValue));
    } else if (command == fAbsLengthCmd) {
        fDetector->SetAbsorptionLength(fAbsLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == fTpbEficCmd) {
        fRun->SetTpbEfficiency(fTpbEficCmd->GetNewDoubleValue(newValue));
    }
      else if (command == fDetSizeCmd)  {
        G4double newLen = fDetSizeCmd->GetNewDoubleValue(newValue);
        //fDetector->SetOuterCellSize(newLen);
        //fDetector->UpdateOuterCellSize(newLen);
        fDetector->outerDiameter = newLen;
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
    if (command == fTopPMTCmd) {
    fDetector->SetTopPMTs(fTopPMTCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fBotPMTCmd) {
    fDetector->SetBotPMTs(fBotPMTCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fSiPMRowsCmd) {
    fDetector->SetSiPMRows(fSiPMRowsCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fPDETopPMTCmd) {
    fDetector->SetPDETopPMT(fPDETopPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fPDEBotPMTCmd) {
    fDetector->SetPDEBotPMT(fPDEBotPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }

  else if (command == fPDESiPMCmd) {
    fDetector->SetPDESiPM(fPDESiPMCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  //else if (command == fReflESRCmd) {
    //fDetector->SetReflESR(fReflESRCmd->GetNewDoubleValue(newValue));
    //fDetector->ApplySensorConfig();
 // }
  else if (command == fReflPMTCmd) {
    fDetector->SetReflPMT(fReflPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fReflSiPMCmd) {
    fDetector->SetReflSiPM(fReflSiPMCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fPMTResponseModeCmd) {
    fDetector->SetPMTResponseMode(newValue);
    fDetector->ApplySensorConfig();
  }
  else if (command == fSiPMResponseModeCmd) {
    fDetector->SetSiPMResponseMode(newValue);
    fDetector->ApplySensorConfig();
  }
  else if (command == fPMTCurveFileCmd) {
    fDetector->SetPMTCurveFile(newValue);
    fDetector->ApplySensorConfig();
  }
  else if (command == fSiPMCurveFileCmd) {
    fDetector->SetSiPMCurveFile(newValue);
    fDetector->ApplySensorConfig();
  }
}
