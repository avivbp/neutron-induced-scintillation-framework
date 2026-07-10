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
/// \file electromagnetic/TestEm5/TestEm5.cc
/// \brief Main program of the electromagnetic/TestEm5 example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"
#include <csignal>
#include "G4RunManagerFactory.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonPhysics.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"
#include "G4Scintillation.hh"
#include "G4OpticalParameters.hh"
#include "G4Cerenkov.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4GeometryTolerance.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void SignalHandler(int) {
 //   G4cout << "\n>>> SIGINT received, aborting current run cleanly...\n" << G4endl;
  //  G4RunManager::GetRunManager()->AbortRun(true); // true = soft abort
//}

int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);
  
  //Use SteppingVerbose with Unit
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);
 
  //std::signal(SIGINT, SignalHandler);
 
  //Creating run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  // Argument modes:
//   ./LArLightSim                         -> interactive, default threads
//   ./LArLightSim 8                       -> interactive, 8 threads
//   ./LArLightSim macro.mac               -> batch macro, default threads
//   ./LArLightSim macro.mac 8             -> batch macro, 8 threads

G4String macroFile = "";

  if (argc == 2) {
    G4String arg1 = argv[1];

    // If argument ends with .mac, treat it as a macro file.
    // Otherwise treat it as number of threads for interactive mode.
    if (arg1.contains(".mac")) {
        macroFile = arg1;
    } else {
        G4int nThreads = G4UIcommand::ConvertToInt(argv[1]);
        runManager->SetNumberOfThreads(nThreads);

        G4cout << "===== LArLightSim is started with "
               << runManager->GetNumberOfThreads()
               << " threads =====" << G4endl;

        ui = new G4UIExecutive(argc, argv);
        }
  }

  if (argc == 3) {
    macroFile = argv[1];

    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);

    G4cout << "===== LArLightSim is started with "
           << runManager->GetNumberOfThreads()
           << " threads =====" << G4endl;
  }
  //set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  //G4PhysListFactory physListFactory;
  //const G4String plName = "QGSP_BERT_HP";
  //G4VModularPhysicsList* pList = physListFactory.GetReferencePhysList(plName);
  //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //pList->RegisterPhysics(opticalPhysics);
  //runManager->SetUserInitialization(pList);

  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  //G4RadioactiveDecayPhysics* radPhysics = new G4RadioactiveDecayPhysics();
  //G4DecayPhysics* decayPhysics = new G4DecayPhysics();
  //G4IonPhysics* ionPhysics = new G4IonPhysics();
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

  auto opticalParams = G4OpticalParameters::Instance();
  
  opticalParams->SetWLSTimeProfile("exponential");
  
  opticalParams->SetScintTrackSecondariesFirst(true);
  
  opticalParams->SetCerenkovMaxPhotonsPerStep(100);
  //opticalParams->SetCerenkovTrackSecondariesFirst(true);
 
  physicsList->RegisterPhysics(opticalPhysics);
  physicsList->RegisterPhysics(new G4EmLowEPPhysics());
  //physicsList->RegisterPhysics(radPhysics);
  //physicsList->RegisterPhysics(decayPhysics);
  //physicsList->RegisterPhysics(ionPhysics);
  runManager->SetUserInitialization(physicsList);

 // runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(detector);

  //set user action classes
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // Initialize visualization before processing macros so /vis commands are available
  // in both interactive sessions and visualization macro runs.
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (macroFile != "") {
  // batch mode
  G4String command = "/control/execute ";
  UImanager->ApplyCommand(command + macroFile);
}
else if (ui) {
  // interactive mode
  ui->SessionStart();
  delete ui;
}
  //job termination
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
