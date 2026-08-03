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
/// \file electromagnetic/TestEm5/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4EmSecondaryParticleType.hh"
#include "G4Neutron.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction* detector, EventAction* EA)
 : G4UserStackingAction(), fDetector(detector), fEventAction(EA)
{
  fStackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  const auto* definition = aTrack->GetDefinition();
  const G4bool isNuclear =
      definition->GetParticleType() == "nucleus" ||
      definition == G4Neutron::NeutronDefinition();
  fNuclearTrack[aTrack->GetTrackID()] = isNuclear;

  // Geant4 11 makes its built-in particle-dependent scintillation and Birks
  // saturation mutually exclusive. Preserve the fitted Birks response, then
  // Bernoulli-thin only LAr scintillation photons whose parent was an ion.
  // For a Poisson photon population this produces the correct Poisson mean and
  // variance at Leff*mean before any optical transport or detector response.
  if (definition == G4OpticalPhoton::OpticalPhotonDefinition() &&
      aTrack->GetCreatorProcess() &&
      aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation") {
    const auto parent = fNuclearTrack.find(aTrack->GetParentID());
    const G4bool isLAr =
        fDetector->IsInsideScintillatingLAr(aTrack->GetPosition());
    if (isLAr && parent != fNuclearTrack.end() && parent->second &&
        G4UniformRand() > fDetector->GetIonScintYieldScale()) {
      return fKill;
    }
  }

  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //keep primary particle
  if (aTrack->GetParentID() == 0) { return fUrgent; }

  //if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
      //if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
        //  auto lv = aTrack->GetLogicalVolumeAtVertex();            // volume at creation
         // auto reg = lv ? lv->GetRegion() : nullptr;
          //if (reg) {
            //const auto& rname = reg->GetName();
            //const auto t = aTrack->GetGlobalTime();
            //if (rname=="InnerLArRegion" && t > fEventAction->tZeroEx)
             // return fKill;
            //if (rname=="LSRegion" && t > fEventAction->tOneEx)
             // return fKill;
          //}
      //}
  return fUrgent;
  //}

  //G4int procID = aTrack->GetCreatorProcess()->GetProcessSubType();
  //G4int modelID = aTrack->GetCreatorModelID();

  //count secondary particles
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
  run->CountParticles(aTrack->GetDefinition());
  /*
  G4cout << "###StackingAction: new " 
         << aTrack->GetDefinition()->GetParticleName()
         << " E(MeV)= " << aTrack->GetKineticEnergy()
         << "  " << aTrack->GetMomentumDirection() << G4endl;
  */
  //
  //energy spectrum of secondaries
  //
  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();

 // if (charge != 0.) {
 //   analysisManager->FillH1(2,energy);
   // analysisManager->FillH1(4,energy);
   // if(procID >= 51 && procID <= 65) { 
     // analysisManager->FillH1(58,energy);
     // analysisManager->FillH1(60,energy);
   // } else if(_AugerElectron == modelID) {
     // analysisManager->FillH1(50,energy);
 //     analysisManager->FillH1(52,energy);
   // } else if(_ePIXE == modelID) {
   //   analysisManager->FillH1(54,energy);
     // analysisManager->FillH1(56,energy);
   // }
 // }

  //if (aTrack->GetDefinition() == G4Gamma::Gamma()) {
   // analysisManager->FillH1(3,energy);
    //analysisManager->FillH1(5,energy);
    //if(procID >= 51 && procID <= 65) { 
     // analysisManager->FillH1(59,energy);
      //analysisManager->FillH1(61,energy);
    //} else if(_Fluorescence == modelID) {
     // analysisManager->FillH1(51,energy);
     // analysisManager->FillH1(53,energy);
    //} else if(_GammaPIXE == modelID) {
     // analysisManager->FillH1(55,energy);
      //analysisManager->FillH1(57,energy);
   // }
  //}  

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (0 < fKillSecondary) {
    if (fKillSecondary == 1) {
      fEventAction->AddEnergy(energy);
    }  
    status = fKill;
  }
    
  return status;
}

void StackingAction::PrepareNewEvent()
{
  fNuclearTrack.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
