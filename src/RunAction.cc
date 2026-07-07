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
/// \file electromagnetic/TestEm5/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include "G4OpticalParameters.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4Threading.hh"   // only if you use IsMasterThread()
#include "G4ios.hh"         // for G4cout / G4endl

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fRun(0), fHistoManager(0)
{
  numPassed = 0;
  std::cout << "Absorber starting x position is : " << fDetector->GetxstartAbs()/CLHEP::cm << " cm" << std::endl;
  std::cout << "Absorber end x position is : " << fDetector->GetxendAbs()/CLHEP::cm << " cm" << std::endl;
 // csvfile.open("results.csv"); 
  //std::cout << "opened csvfile in run action" << std::endl;

 // csvfile.open("freePathDist.csv");
  //std::cout << "opened 2nd csvfile in run action" << std::endl;
  // Book predefined histograms

 // auto man = G4AnalysisManager::Instance();
  //man->CreateNtuple("TestEm5Aviv", "Edep and TrackL");
  //man->CreateNtupleDColumn("Ek-Electron");
  //man->CreateNtupleDColumn("DeltaE");
  //man->FinishNtuple();
  //fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)

{
  csvfile << "process,distance" << "," << std::endl;
  // show Rndm status
  if (IsMaster()) {
      G4Random::showEngineStatus();
      auto* pm = G4Electron::Electron()->GetProcessManager();
      auto* p  = pm ? pm->GetProcess("Cerenkov") : nullptr;

      if (p && pm) {
          G4bool active = pm->GetProcessActivation(p);
          G4cout << "[Run] Cerenkov process: "
                 << (active ? "ACTIVE" : "INACTIVE")
                 << G4endl;
      }

  }
 
  // keep run condition
  if ( fPrimary ) { 
    G4ParticleDefinition* particle 
      = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }
  
  //histograms
  //        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* grun)
{  
  csvfile.close();
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun();    
      
  // save histograms

  // show Rndm status
  if (IsMaster()){

      const auto* run = static_cast<const Run*>(grun);
      // --- TPB validation (minimal) ---
//      const auto made_in_inner   = run->GetChkCreatedInnerScint();
//      const auto wls_abs_from_in = run->GetChkWLSAbsFromInner();
//      const auto wls_emit_from_in= run->GetChkWLSEmitFromInner();
//
//      G4cout << "\n=== TPB VALIDATION (Scint in innerCell → WLS in TPB) ===\n"
//       << "created (Scint, innerCell)        : " << made_in_inner << G4endl
//       << "WLS absorbed from those tracks    : " << wls_abs_from_in << G4endl
//       << "WLS emitted photons from those    : " << wls_emit_from_in << G4endl
//       << std::fixed << std::setprecision(4)
//       << "ratio (WLS_abs_from_inner / created_inner) = "
//       << (made_in_inner ? double(wls_abs_from_in)/made_in_inner : 0.0) << G4endl
//       << "ratio (WLS_emit_from_inner / WLS_abs_from_inner) = "
//       << (wls_abs_from_in ? double(wls_emit_from_in)/wls_abs_from_in : 0.0) << G4endl;
//

//    {
//      std::ofstream out("photHits.csv", std::ios::out | std::ios::trunc);
//      out << "eventID,trackID,surface,u,v,tAbs\n";       // header
//      //out << run->GetPhotHitsCSV();                      // payload
//    }
//    {
//      std::ofstream out("evtStats.csv", std::ios::out | std::ios::trunc);
//      out << "eventID,detector,tOne,numElasticSensitive,numInelasticSensitive,ExtScatter,numInelastic,good\n";
//      //out << run->GetEvtStatsCSV();
//    }
      G4Random::showEngineStatus();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
