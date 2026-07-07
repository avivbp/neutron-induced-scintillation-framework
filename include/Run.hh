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
/// \file electromagnetic/TestEm5/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
//#include "CsvWriterThreadSafe.hh"
//#include <memory>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <unordered_set>
class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{

 // --- TPB validation (minimal) ---
private:
  G4int fChk_CreatedInnerScint = 0;  // # optical photons born by Scintillation in innerCell
  G4int fChk_WLSAbs_fromInner  = 0;  // # OpWLS absorptions of those same tracks
  G4int fChk_WLSEmit_fromInner = 0;  // # OpWLS-emitted secondaries from those absorptions

public:
  inline void AddChkCreatedInnerScint(int n=1){ fChk_CreatedInnerScint += n; }
  inline void AddChkWLSAbsFromInner(int n=1){  fChk_WLSAbs_fromInner  += n; }
  inline void AddChkWLSEmitFromInner(int n=1){ fChk_WLSEmit_fromInner += n; }

  inline G4int GetChkCreatedInnerScint() const { return fChk_CreatedInnerScint; }
  inline G4int GetChkWLSAbsFromInner()  const { return fChk_WLSAbs_fromInner; }
  inline G4int GetChkWLSEmitFromInner() const { return fChk_WLSEmit_fromInner; }

  std::unordered_set<int> completedEvents;
  inline bool completedEvent(int eventID) {
    return completedEvents.insert(eventID).second; // true only the first time
  }

 public:
   G4int numInteractions;
   G4int numElastic;
   G4int numInelastic;
   G4int numNNPrime;
   G4int numNGamma;
   G4int numNP;
   G4int numSuc;
   G4int numBackground;
   G4int single_elastic; 
   G4int multiple_elastic; 
   G4int elastic_NNPrime;
   G4int multiple_NNPrime;
   G4int no_interaction;
   G4int numPassed;
   G4double dEdx;
   G4double tpb;
   std::ofstream csvf;
   std::ofstream csvfi;
   std::ofstream csvfil;
   std::ofstream csfile;
   Run(DetectorConstruction*);
  ~Run();

 public:
   
   void SetTpbEfficiency(G4double value) {
       tpb = value;
   }
 
   // Appenders used by EventAction at EndOfEvent
   void AppendPhotHitsCSV(const std::string& s) { fPhotHitsCSV += s; }
   void AppendEvtStatsCSV(const std::string& s) { fEvtStatsCSV += s; }

   const std::string& GetPhotHitsCSV() const { return fPhotHitsCSV; }
   const std::string& GetEvtStatsCSV() const { return fEvtStatsCSV; }

   void addNumInteractions(){numInteractions+=1;}; 
   void addNumElastic(){numElastic+=1;};
   void addNumInelastic(){numInelastic+=1;};
   void addNumNNPrime(){numNNPrime+=1;};
   void addNumNGamma(){numNGamma+=1;};
   void addNumNP(){numNP+=1;};
   void SetPrimary(G4ParticleDefinition* particle, G4double energy);

   void AddEnergy (G4double edep)
              {fEnergyDeposit += edep; fEnergyDeposit2 += edep*edep;};

   void AddTrakLenCharg (G4double length)
              {fTrakLenCharged += length; fTrakLenCharged2 += length*length;};

   void AddTrakLenNeutr (G4double length)
              {fTrakLenNeutral += length; fTrakLenNeutral2 += length*length;};

   void AddMscProjTheta (G4double theta)
              {if (std::abs(theta) <= fMscThetaCentral) { fMscEntryCentral++;
                 fMscProjecTheta += theta;  fMscProjecTheta2 += theta*theta;}
              };

   void CountStepsCharg (G4int nSteps)
              {fNbStepsCharged += nSteps; fNbStepsCharged2 += nSteps*nSteps;};

   void CountStepsNeutr (G4int nSteps)
              {fNbStepsNeutral += nSteps; fNbStepsNeutral2 += nSteps*nSteps;};

   void CountParticles (G4ParticleDefinition* part)
              {     if (part == G4Gamma::Gamma())       fNbGamma++ ;
               else if (part == G4Electron::Electron()) fNbElect++ ;
               else if (part == G4Positron::Positron()) fNbPosit++ ; };

   void CountTransmit (G4int flag)
              {     if (flag == 1)  fTransmit[0]++;
               else if (flag == 2) {fTransmit[0]++; fTransmit[1]++; }};

   void CountReflect (G4int flag)
              {     if (flag == 1)  fReflect[0]++;
               else if (flag == 2) {fReflect[0]++; fReflect[1]++; }};
    
   void AddEnergyLeak (G4double eleak, G4int index)
            {fEnergyLeak[index] += eleak; fEnergyLeak2[index] += eleak*eleak;};
            
   G4double ComputeMscHighland();
               
   virtual void Merge(const G4Run*);
   
   void EndOfRun();   

 private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin;

    std::string fPhotHitsCSV;
    std::string fEvtStatsCSV;
                           
    G4double fEnergyDeposit,  fEnergyDeposit2;
    G4double fTrakLenCharged, fTrakLenCharged2;
    G4double fTrakLenNeutral, fTrakLenNeutral2;
    G4double fNbStepsCharged, fNbStepsCharged2;
    G4double fNbStepsNeutral, fNbStepsNeutral2;
    G4double fMscProjecTheta, fMscProjecTheta2;
    G4double fMscThetaCentral;
    
    G4int    fNbGamma, fNbElect, fNbPosit;
    G4int    fTransmit[2],   fReflect[2];
    G4int    fMscEntryCentral;
    
    G4double fEnergyLeak[2],  fEnergyLeak2[2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

