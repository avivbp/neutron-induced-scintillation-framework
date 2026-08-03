#ifndef DetectorMessenger_HH
#define DetectorMessenger_HH

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4RunManager.hh"

class DetectorConstruction;
class Run;

class DetectorMessenger : public G4UImessenger {
public:
    DetectorMessenger(DetectorConstruction* det, Run* run);
    ~DetectorMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    DetectorConstruction* fDetector;
    G4UIdirectory* fDetDir = nullptr;
    Run* fRun;

    // integer config
    G4UIcmdWithAnInteger* fTopPMTCmd  = nullptr;
    G4UIcmdWithAnInteger* fBotPMTCmd  = nullptr;
    G4UIcmdWithAnInteger* fSiPMRowsCmd = nullptr;

    // doubles
    G4UIcmdWithADouble* fPDETopPMTCmd    = nullptr;
    G4UIcmdWithADouble* fPDEBotPMTCmd    = nullptr;
    G4UIcmdWithADouble* fPDESiPMCmd   = nullptr;

    G4UIcmdWithADouble* fReflESRCmd   = nullptr;
    G4UIcmdWithADouble* fReflPMTCmd   = nullptr;
    G4UIcmdWithADouble* fReflSiPMCmd  = nullptr;

    G4UIcmdWithAString* fPMTResponseModeCmd  = nullptr;
    G4UIcmdWithAString* fSiPMResponseModeCmd = nullptr;
    G4UIcmdWithAString* fPMTCurveFileCmd     = nullptr;
    G4UIcmdWithAString* fSiPMCurveFileCmd    = nullptr;

    G4UIcmdWithADouble* fInnerRadiusCmd      = nullptr;
    G4UIcmdWithADouble* fInnerHeightCmd      = nullptr;
    G4UIcmdWithADouble* fOuterDiameterCmd    = nullptr;
    G4UIcmdWithADouble* fOuterHeightCmd      = nullptr;
    G4UIcmdWithABool* fOuterScintillationCmd = nullptr;
    G4UIcmdWithAString* fNeutronDetectorsCmd = nullptr;
    G4UIcmdWithADouble* fPMTSideCmd          = nullptr;
    G4UIcmdWithAString* fTOFWindowCmd        = nullptr;
    G4UIcmdWithADouble* fPrimaryEnergyCmd    = nullptr;
    G4UIcmdWithADouble* fSiPMTileCmd         = nullptr;
    G4UIcmdWithADouble* fAbsLengthCmCmd      = nullptr;
    G4UIcmdWithADouble* fBirksMmPerMeVCmd    = nullptr;
    G4UIcmdWithADouble* fScintYieldCmd       = nullptr;
    G4UIcmdWithADouble* fScintScaleCmd       = nullptr;
    G4UIcmdWithADouble* fIonScintScaleCmd    = nullptr;
    G4UIcmdWithADouble* fFastTimeCmd         = nullptr;
    G4UIcmdWithADouble* fSlowTimeCmd         = nullptr;
    G4UIcmdWithADouble* fFastFractionCmd     = nullptr;
    G4UIcmdWithADouble* fTPBEfficiencyCmd    = nullptr;
    G4UIcmdWithADouble* fReflectorReflCmd    = nullptr;

    G4UIdirectory* fSimDir;
    G4UIcmdWithADoubleAndUnit* fBirksCmd;
    G4UIcmdWithADoubleAndUnit* fAbsLengthCmd;
    G4UIcmdWithADouble* fTpbEficCmd;
    G4UIcmdWithADoubleAndUnit* fDetSizeCmd;
};

#endif
