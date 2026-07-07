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
/// \file electromagnetic/TestEm5/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "EventAction.hh"
#include <cmath>
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Scintillation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessVector.hh"
#include "G4Cerenkov.hh"
#include "G4OpBoundaryProcess.hh"
class DetectorConstruction;
class RunAction;
class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{

  private:
   G4Material*      fTPBMat        = nullptr;
   G4LogicalVolume* fInnerCellLV   = nullptr;

  public:
   std::ofstream csvfile;
   G4double tpbEfic;
   SteppingAction(DetectorConstruction*,EventAction*);
  ~SteppingAction();

   virtual void UserSteppingAction(const G4Step*);
   G4bool withinArapuca(const G4Step* aStep,G4double height);
   G4bool withinPMTs(const G4Step* aStep,G4double height);

   G4double calcStepLength(const G4Step* aStep){
       const G4ThreeVector & pre = aStep->GetPreStepPoint()->GetPosition()/CLHEP::cm;
       const G4ThreeVector & post = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
       G4double xDiff = post.getX() - pre.getX();
       G4double yDiff = post.getY() - pre.getY();
       G4double zDiff = post.getZ() - pre.getZ();
       G4double stepLength = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);
       return stepLength;
   };

   G4double calcAngle(const G4ThreeVector & v1, const G4ThreeVector & v2){
       G4double dotProduct = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
       G4double angle = acos(dotProduct);
       return angle;
   }

   G4double calcRecoilEnergy(G4double E0,G4int A,G4double scatAngle){
       return (2*E0*A * (1 - cos(scatAngle)) / pow(A+1,2))/CLHEP::keV;
   }

   void printEventStats(const G4Step* aStep,G4int numMax){

       G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
       std::cout << "current eventID = " << eventID << std::endl;
       std::cout << "current step number = " << aStep->GetTrack()->GetCurrentStepNumber() << std::endl;
       //std::cout << "time before step = " << aStep->GetPreStepPoint()->GetGlobalTime()/CLHEP::ns << " ns" << std::endl;
       const G4String & initialParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
       G4double preEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;
       G4double postEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV;

       std::cout << "initial particle : " << initialParticleName << " with energy before step: " << preEnergy << " MeV " << ", and energy after step: " << postEnergy << " MeV" <<  std::endl;
       const G4String & volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
       std::cout << "particle interaction in " << volumeName << std::endl;
       //std::cout << "initial particle trackID = " << aStep->GetTrack()->GetTrackID() << std::endl;
       const auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
       if ((process->GetProcessName() == "neutronInelastic") || (process->GetProcessName() == "nCapture")){
           auto hadronicProcess = dynamic_cast<const G4HadronicProcess*>(process);
           G4Nucleus nucleus = *(hadronicProcess->GetTargetNucleus());
           auto name = nucleus.GetIsotope()->GetName();
           std::cout << "target nucleus name: " << name << std::endl;
       }
       else if (process->GetProcessName() == "Transportation"){
           std::cout << "post step position: " << aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm << std::endl;
       }

       std::cout << "process name: " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
       std::cout << "secondaries:" << std::endl;
  
       const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
       G4int size  = (int) (secondaries->size());

       if (numMax){
           if (size > numMax){
               size = numMax;
           }
       }


       for(int i = 0; i < size; i++){ 
           const G4Track* secondary = (*secondaries)[i];
           const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
           std::cout << i << " particle name = " << particleName << std::endl;
           std::cout << "parentID = " << secondary->GetParentID() << ", trackID = " << secondary->GetTrackID() << std::endl;
           std::cout << " momentum direction: " << secondary->GetMomentumDirection() << std::endl;
           std::cout << "kinetic energy = " << secondary->GetKineticEnergy()/CLHEP::MeV << " MeV" <<  std::endl;
           G4int charge = secondary->GetParticleDefinition()->GetAtomicNumber();
           if (charge > 4){
               int Z = charge;
               int A = secondary->GetParticleDefinition()->GetAtomicMass();
               std::cout << "Z = " << Z << ", and A = " << A << std::endl;
               G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
               std::cout << " recoil nucleus angle = " << Beta*180/CLHEP::pi << " degrees" << std::endl;
               G4double alpha = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());
               std::cout << "neutron scatter angle = " << alpha*180/CLHEP::pi << " degrees" << std::endl;
           }
       }
       std::cout<<std::endl;
  };
 
 std::string GetBoundaryProcessStatusString(G4OpBoundaryProcessStatus status) {
    switch (status) {
        case Undefined: return "Undefined";
        case Transmission: return "Transmission";
        case FresnelRefraction: return "FresnelRefraction";
        case FresnelReflection: return "FresnelReflection";
        case TotalInternalReflection: return "TotalInternalReflection";
        case LambertianReflection: return "LambertianReflection";
        case LobeReflection: return "LobeReflection";
        case SpikeReflection: return "SpikeReflection";
        case Absorption: return "Absorption";
        case Detection: return "Detection";
        case NotAtBoundary: return "NotAtBoundary";
        case SameMaterial: return "SameMaterial";
        case StepTooSmall: return "StepTooSmall";
        case NoRINDEX: return "NoRINDEX";
        case PolishedLumirrorAirReflection: return "PolishedLumirrorAirReflection";
        case PolishedLumirrorGlueReflection: return "PolishedLumirrorGlueReflection";
        case PolishedAirReflection: return "PolishedAirReflection";
        case PolishedTeflonAirReflection: return "PolishedTeflonAirReflection";
        case PolishedTiOAirReflection: return "PolishedTiOAirReflection";
        case PolishedTyvekAirReflection: return "PolishedTyvekAirReflection";
        case PolishedVM2000AirReflection: return "PolishedVM2000AirReflection";
        case PolishedVM2000GlueReflection: return "PolishedVM2000GlueReflection";
        case EtchedLumirrorAirReflection: return "EtchedLumirrorAirReflection";
        case EtchedLumirrorGlueReflection: return "EtchedLumirrorGlueReflection";
        case EtchedAirReflection: return "EtchedAirReflection";
        case EtchedTeflonAirReflection: return "EtchedTeflonAirReflection";
        case EtchedTiOAirReflection: return "EtchedTiOAirReflection";
        case EtchedTyvekAirReflection: return "EtchedTyvekAirReflection";
        case EtchedVM2000AirReflection: return "EtchedVM2000AirReflection";
        case EtchedVM2000GlueReflection: return "EtchedVM2000GlueReflection";
        case GroundLumirrorAirReflection: return "GroundLumirrorAirReflection";
        case GroundLumirrorGlueReflection: return "GroundLumirrorGlueReflection";
        case GroundAirReflection: return "GroundAirReflection";
        case GroundTeflonAirReflection: return "GroundTeflonAirReflection";
        case GroundTiOAirReflection: return "GroundTiOAirReflection";
        case GroundTyvekAirReflection: return "GroundTyvekAirReflection";
        case GroundVM2000AirReflection: return "GroundVM2000AirReflection";
        case GroundVM2000GlueReflection: return "GroundVM2000GlueReflection";
        default: return "Unknown";
    }
}
 // Generic function to find if an element of any type exists in list
  template <typename T>
  bool contains(std::list<T> & listOfElements, const T & element)
  {
      // Find the iterator if element in list
      auto it = std::find(listOfElements.begin(), listOfElements.end(), element);
      //return if iterator points to end or not. It points to end then it means element
      // does not exists in list
      return it != listOfElements.end();
  }


  private:
    DetectorConstruction* fDetector;
    EventAction*          fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
