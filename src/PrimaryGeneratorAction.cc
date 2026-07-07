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
/// \file electromagnetic/TestEm5/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
  :G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),fDetector(DC),fRndmBeam(1),fGunMessenger(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;  
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematic()
{    
  // default particle kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4ParticleDefinition *particle = particleTable->FindParticle("geantino");
  //fParticleGun->SetParticleMomentum(0.*MeV);
  //fParticleGun->SetParticleDefinition(particle);

  G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
  //static G4ParticleDefinition* particle = G4OpticalPhoton::OpticalPhotonDefinition();
  //G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  //G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
  //G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(18,40,0);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(2.5*MeV);
  //fParticleGun->SetParticleEnergy(117.8*keV);
  //fParticleGun->SetParticleEnergy(511.0*keV);
  //fParticleGun->SetParticlePolarization(G4ThreeVector(0,1,0));
  //fPolarized = true;
  //G4double x0 = -0.5*(fDetector->GetWorldSizeX());
  G4double x0 = -0.1*m;
  x0 = -3.*m;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, 0.0, 0.0));
  //G4double xDir = 2*G4UniformRand()-1.;
  //G4double yDir = 2*G4UniformRand()-1.;
  //G4double zDir = 2*G4UniformRand()-1.;
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
 // std::cout << "particle starting positioN:" << x0 << std::endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();

  if(particle == G4Geantino::Geantino())
  {
      G4int Z = 3;
      G4int A = 7;

      G4double charge   = 0.*eplus;
      G4double energy = 0.*keV;

      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,energy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(charge);
      //fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  G4double cylinderHeight = 7.87*cm;
  G4double cylinderRadius = 3.8*cm;

  //G4double xDir = 2*G4UniformRand()-1.;
  G4double yDir = (2*G4UniformRand()-1.)*0.045;
  G4double zDir = (2*G4UniformRand()-1.)*0.045;

  if (fParticleGun->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    if (fPolarized)
      SetOptPhotonPolar(fPolarization);
    else
      SetOptPhotonPolar();
  }


  if(0 == anEvent->GetEventID()) {
    //G4double x0 = -0.5*(fDetector->GetWorldSizeX());
    G4double x0 = -3.*m;
    //x0 = 0.*m;
    //G4double y0 = -2.*m;
    //x0 = -10.*m;
    fParticleGun->SetParticlePosition(G4ThreeVector(x0, 0.0, 0.0));  
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  }
  //this function is called at the begining of event
  //
  //randomize the beam, if requested.
  fRndmBeam = 0;
  if (fRndmBeam > 0.) 
    {
      G4ThreeVector oldPosition = fParticleGun->GetParticlePosition();    
      G4double rbeam = 7.*cm;
      G4double x0 = -1.*m;
      //G4double y0 = (2*G4UniformRand()-1.)*rbeam;
      //G4double z0 = (2*G4UniformRand()-1.)*rbeam;
     // std::cout << "primary starting position:" << G4ThreeVector(x0,y0,z0)/CLHEP::cm << std::endl;

      // pencil beam
      //fParticleGun->SetParticlePosition(G4ThreeVector(x0,0.0,0.0));
      G4double rndEnergy = (1.45 + (2*G4UniformRand()-1.)*0.085) * MeV;
      fParticleGun->SetParticleEnergy(rndEnergy);
      fParticleGun->SetParticlePosition(G4ThreeVector(x0,0,0));      
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }
  else  fParticleGun->GeneratePrimaryVertex(anEvent); 
}

void PrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  if (fParticleGun->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    G4ExceptionDescription ed;
    ed << "The particleGun is not an opticalphoton.";
    G4Exception("PrimaryGeneratorAction::SetOptPhotonPolar", "OpNovice2_004", JustWarning, ed);
    return;
  }

  fPolarized = true;
  fPolarization = angle;

  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2 = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
