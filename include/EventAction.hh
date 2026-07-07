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
/// \file electromagnetic/TestEm5/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <list>
#include <sstream>      // NEW
#include <string>       // NEW
#include <vector>       // NEW
#include <unordered_set>
#include <unordered_map>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{

struct WLSHitSig { char s; int u,v,t; };
inline bool nearly_same(const WLSHitSig& a, const WLSHitSig& b) {
  if (a.s != b.s) return false;
  auto d = [](int x,int y){ return std::abs(x-y); };
  return d(a.u,b.u) <= 1 && d(a.v,b.v) <= 1 && d(a.t,b.t) <= 1; // or <=1 for time too
}

private:
  std::unordered_map<int, WLSHitSig> fLastWLSSigByTrack; // trackID -> last signature
public:
  inline void ClearWLSDeDupe() { fLastWLSSigByTrack.clear(); }
  inline bool ShouldLogWLS(int tid, char s, int u, int v, int t) {
    WLSHitSig cur{ s,u,v,t };
    auto it = fLastWLSSigByTrack.find(tid);
    if (it != fLastWLSSigByTrack.end() && nearly_same(it->second, cur)) return false;
    fLastWLSSigByTrack[tid] = cur;
    return true;
  }

// Tracks (IDs) of optical photons born by Scintillation in innerCell (per event)
private:
  std::unordered_set<int> fBornInnerScint;

public:
  inline void ClearBornInnerScint()              { fBornInnerScint.clear(); }
  inline void MarkBornInnerScint(int tid)        { fBornInnerScint.insert(tid); }
  inline bool IsBornInnerScint(int tid) const    { return fBornInnerScint.count(tid) > 0; }

  public:
    G4int numInteractions;
    G4bool wroteToFile;
    G4double sumStepLengths;
    G4int totNeutronNum;
    G4int totProtonNum;
    G4double gammaEnergyDeposit;
    G4double secondarygammasDeposit;
    G4bool passed;
    G4bool performedNPrime;
    G4int numElastic;
    G4int numInelasticSensitive;
    G4bool detected;
    G4bool scatteredElastically;
    G4int numElasticSensitive;
    G4bool scatteredInelastically;
    G4double nucleusRecoilEnergy;
    G4double nScatterAngle;
    G4double tZero;
    G4double tOne;
    G4double outerCellScatterAngle;
    G4double outerCellEDep;
    G4double CryoScatterAngle;
    G4double CryoEDep;
    G4double innerLayerEDep;
    G4double innerCellEDep;
    G4double incomingAng;
    G4bool scatteredNotSensitive;
    G4bool extInelastic;
    G4bool Coincedence;
    G4bool fail;
    G4bool aborted;
    G4bool nCapture;
    G4bool secondaryNeutron;
    G4bool NNPrime;
    G4bool NP;
    G4bool NNPrimeP;
    G4bool N2N;
    G4bool N3N;
    std::ostringstream buf_hits;
    G4int ExtScatter;
    G4int CryoScatter;
    G4int innerLayerScatter;
    G4int numPhotLS;
    CLHEP::Hep3Vector pos0;
    CLHEP::Hep3Vector pos1;
    G4double fiducialIncomingEn;
    G4double experimentalRecoilEnergy;
    G4double tZeroEx;
    G4double tOneEx;
    G4double tHelp;
    G4double escapeEn;
    G4int nScint;
    G4int nCher;
    G4int numPhotons;
    G4int numPE;
    G4int largestTrackID;
    G4int lastTracked;
    G4int numReachedUp;
    G4int numReachedDown;
    G4int num;
    G4int numInelastic;
    G4int nnPrime;
    G4int numSteps;
    G4int numSurface;
    G4int numBases;
    G4int extScatterBefore;
    G4int extScatterAfter;
    G4double totEnergy;
    G4double tUp;
    G4double tDown;
    G4double eDep;
    G4double enDep;
    G4double sumStepLength;
    G4double coincidenceTime;
    G4String detector;
    std::ofstream csvfile;
    std::list<int> parentIDS;
    EventAction();
   ~EventAction();

  public:

    void ClearSeenUV() { fSeenUV.clear(); }
    bool MarkSeenUV(G4int trackID) { return fSeenUV.insert(trackID).second; }

    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
 
    inline void BufferHitRow(const std::vector<std::string>& fields) {
    for (size_t i=0;i<fields.size();++i) { if(i) buf_hits<<','; buf_hits<<fields[i]; }
    buf_hits << '\n';
    }
    void numInteractionsPP(){numInteractions += 1;};
    void setNPrime(G4bool b){performedNPrime = b;};
    void resetNumInteractions(){numInteractions = 0;};
    void setPassed(G4bool b){passed = b;};
    void addStepLength(G4double stepL){sumStepLengths+=stepL;};
    void addNeutronNum(G4int num){totNeutronNum += num;};
    void addProtonNum(G4int num){totProtonNum += num;};
    void resetStepLengths(){sumStepLengths = 0.;};
    void resetNeutronNum(){totNeutronNum = 0;};
    void resetProtonNum(){totProtonNum = 0;};
    void setWroteToFile(G4bool boo ){wroteToFile = boo;};
    void addEDep(G4double gammaEDep){gammaEnergyDeposit+=gammaEDep;};
    void addDeposit(G4double gammaEDep){secondarygammasDeposit+=gammaEDep;};
    void elasticPP(){numElastic += 1;};

    void AddEnergy      (G4double edep)   {fEnergyDeposit  += edep;};
    void AddTrakLenCharg(G4double length) {fTrakLenCharged += length;};
    void AddTrakLenNeutr(G4double length) {fTrakLenNeutral += length;};
    
    void CountStepsCharg ()               {fNbStepsCharged++ ;};
    void CountStepsNeutr ()               {fNbStepsNeutral++ ;};
    
    void SetTransmitFlag (G4int flag) 
                           {if (flag > fTransmitFlag) fTransmitFlag = flag;};
    void SetReflectFlag  (G4int flag) 
                           {if (flag > fReflectFlag)   fReflectFlag = flag;};
                                             
        
  private:
    std::unordered_set<G4int> fSeenUV; // TrackIDs of UV photons already counted
    G4double fEnergyDeposit;
    G4double fTrakLenCharged, fTrakLenNeutral;
    G4int    fNbStepsCharged, fNbStepsNeutral;
    G4int    fTransmitFlag,   fReflectFlag;        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
