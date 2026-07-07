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
/// \file electromagnetic/TestEm5/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include <sstream>
#include <unordered_map>
#include "G4RunManager.hh"
#include "G4Event.hh"
// --- Per-thread CSV output (drop-in) ---
#include <fstream>
#include <memory>
#include "G4Threading.hh"
#include "G4UImanager.hh"

namespace {
  // One set of files *per worker thread* (no mutex required).
  thread_local std::unique_ptr<std::ofstream> tls_hits;
  thread_local std::unique_ptr<std::ofstream> tls_stats;
  thread_local bool tls_opened = false;

  inline int ThreadIdForName() {
    // Geant4 helper: 0..N-1 for workers; for sequential builds this returns 0.
    // If unavailable in your version, you can fall back to 0 for ST and
    // std::hash<std::thread::id>{}(std::this_thread::get_id()) for MT.
    return G4Threading::G4GetThreadId();
  }

  inline void OpenPerThreadFilesOnce() {
    if (tls_opened) return;
    tls_opened = true;

    const int tid = ThreadIdForName();
    // Per-thread filenames
    std::ostringstream nh, ns;
    nh << "photHits_t" << tid << ".csv";
    ns << "evtStats_t" << tid << ".csv";

    // Open in trunc mode so each run starts fresh for this worker
    tls_hits = std::make_unique<std::ofstream>(nh.str(), std::ios::out | std::ios::trunc);
    tls_stats= std::make_unique<std::ofstream>(ns.str(), std::ios::out | std::ios::trunc);

    // Headers (once per thread-file)
    (*tls_hits)  << "event,track,surfC,u,v,t_10ps";
    (*tls_stats) << "event,detector,tOne,numElasticSensitive,numInelasticSensitive,ExtScatter,numInelastic,good\n";

    tls_hits->flush();
    tls_stats->flush();
  }

  inline void AppendHitsPerThread(const std::string& payload) {
    // payload may contain several '\n'-terminated rows
    (*tls_hits) << payload;
    tls_hits->flush();      // flush each event for Ctrl-C safety
  }

  inline void AppendStatsPerThread(const std::string& oneLine) {
    (*tls_stats) << oneLine << '\n'; 
    tls_stats->flush();
  }
}

//G4bool wroteDistance = false;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
 fEnergyDeposit(0.),
 fTrakLenCharged(0.), fTrakLenNeutral(0.),
 fNbStepsCharged(0), fNbStepsNeutral(0),
 fTransmitFlag(0), fReflectFlag(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt )
{
 //std::cout << "in begin of event :" << evt->GetEventID() << std::endl;
 // initialisation per event
 //if (evt->GetEventID() >=17994) { // or any condition
    //G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 2");
 //}
 setWroteToFile(false);
 resetStepLengths();
 resetNumInteractions();
 resetNeutronNum();
 resetProtonNum();
 setNPrime(false);
 gammaEnergyDeposit = secondarygammasDeposit = 0.;
// std::cout << passed << std::endl;
// setPassed(false);
// std::cout << passed << std::endl;
 fEnergyDeposit  = experimentalRecoilEnergy = sumStepLength = 0.;
 pos0 = pos1 = G4ThreeVector(0.,0.,0.);
 fTrakLenCharged = fTrakLenNeutral = tHelp = eDep = enDep = tUp = tDown = coincidenceTime = 0.; 
 fNbStepsCharged = fNbStepsNeutral = numElastic = numElasticSensitive = detected = scatteredElastically = scatteredInelastically = nucleusRecoilEnergy = nScatterAngle = numInelasticSensitive = 0;
 nScint = nCher = numPhotons = largestTrackID = lastTracked = numReachedUp = numReachedDown = numPE = num = numSteps = totEnergy = ExtScatter = CryoScatter = innerLayerScatter = 0;
 fTransmitFlag   = fReflectFlag    = scatteredNotSensitive = extInelastic = fail = numInelastic = nnPrime = numSurface = numBases = aborted = nCapture = secondaryNeutron = 0;
 outerCellScatterAngle = CryoScatterAngle = outerCellEDep = CryoEDep = innerLayerEDep = innerCellEDep = fiducialIncomingEn = escapeEn = incomingAng = 0.;
 tOne = tZero = tZeroEx = tOneEx = 10000.*CLHEP::ns;
 extScatterBefore = extScatterAfter = NNPrime = NP = NNPrimeP = N2N = N3N = numPhotLS = 0;
 detector = ""; 
 Coincedence = true;
 std::ostringstream().swap(buf_hits);
 ClearWLSDeDupe();
 //ClearBornInnerScint();
 //ClearSeenUV();

 //std::cout << "after zeroing in begin eventAction " << std::endl;
 //csvfile.open("scintCheck.csv");
 //csvfile << "energy" << "," << "time" << "," << "theta" << "," << "phi" << "," << "xpos" << "," << "ypos" << "," << "zpos" << "," << std::endl;
 //csvfile.close();

 //csvfile.open("primaryPos.csv");
 //csvfile << "eventID" << "," << "xpos" << "," << "ypos" << "," << "zpos" << "," << std::endl;
 //csvfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int eventID = evt->GetEventID();

  //std::cout << "end of eventAction number: " << eventID << std::endl;

  run->completedEvent(eventID);
  
  scatteredNotSensitive = extInelastic || ExtScatter || CryoScatter || innerLayerScatter;
  //std::cout << "after completedEvent" << std::endl;
  //std::cout << "num photons reached PMT = " << num << std::endl;
  //std::cout << "num photons reached top PMT array = " << numReachedUp << std::endl;
  //std::cout << "num photons reached bottom PMT = " << numReachedDown << std::endl;
  //std::cout << "tot num = " << numReachedDown + numReachedUp << std::endl;
  //if (tZeroEx != 10000){
      //std::cout << std::endl;
      //std::cout << "end of event, tZeroEx = " << tZeroEx << " ns" << std::endl;
      //std::cout << std::endl;
  //}

  //run->fCsv.FlushBuffered();

  const PrimaryGeneratorAction* prim = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const DetectorConstruction* det = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double energy = prim->GetParticleGun()->GetParticleEnergy()/CLHEP::MeV;

  G4int nPMTsTop = det->GetTopPMTs();
  G4int nPMTsBot = det->GetBotPMTs();
  G4int nSiPMRows = det->GetSiPMRows();
  //std::cout << "event #" << eventID << ", nPMTsTop = " << nPMTsTop << ", nPMTsBot = " << nPMTsBot << ", nSiPMRows = " << nSiPMRows << std::endl;
 
 
  G4bool good = 0;
  // if only elastic, (n,n') happened then good
  if (!N3N && !N2N && !NNPrimeP && !NP){
      good = 1;
  }

  //std::stringstream ss2;
  //ss2 << "scatterStats_" << energy << "_MeV_" << det->outerDiameter/CLHEP::cm << "_cm.csv";
  //std::string filename2 = ss2.str();

  //run->csvfi.open(filename2,std::ios_base::app);
  //run->csvfi << eventID << "," << ExtScatter << "," << CryoScatter << "," << innerLayerScatter << "," << outerCellScatterAngle << "," << CryoScatterAngle << "," << outerCellEDep << "," << CryoEDep << "," << innerLayerEDep << "," << numElasticSensitive << "," << numInelasticSensitive << "," << innerCellEDep << "," << fiducialIncomingEn << "," << extScatterBefore << "," << extScatterAfter << "," << escapeEn << "," << incomingAng << "," << good << "," << std::endl;
  //run->csvfi.close();

  // to be changed based on single elastic scatter tof
  G4int earlyCutoff = 40.*CLHEP::ns;
  G4int lateCutoff = 50.*CLHEP::ns;
  G4double tof = tOne - tZeroEx; // using tOne to save runtime on optical photons, tOne - tOneEx ~ O(0.1 [ns])

  //if (tOne < 1000 && tZeroEx < 1000) std::cout << "tOne = " << tOne << ", tZeroEx = " << tZeroEx << ", tof = " << tOne - tZeroEx << std::endl;
  //std::cout << "after tof calc" << std::endl;

  G4int earlyCutoffNoScint = 40.*CLHEP::ns;
  G4int lateCutoffNoScint = 50.*CLHEP::ns;
  G4double tofNoScint = tOne - tZero;
  // single elastic scatter, neutron detected by liquid scintillator
  if ((tofNoScint > earlyCutoffNoScint) && (tofNoScint < lateCutoffNoScint) && detected){
      run->csvfi.open("tof.csv",std::ios_base::app);
      run->csvfi << eventID << "," << tofNoScint << "," << detector << "," << numElasticSensitive << "," << ExtScatter << "," << CryoScatter << "," << innerLayerScatter << "," << std::endl;
      run->csvfi.close();
  }

  //std::stringstream ss;
  //ss << "numPE_" << det->absorptionLength[0]/CLHEP::cm << "_cm_absLen_" << energy << "_MeV.csv";
  //std::string filename = ss.str();

  std::stringstream ss;
  ss << "numPE_" << nPMTsTop << "_topPMTs_" << nPMTsBot << "_botPMTs_" << nSiPMRows << "_SiPMRows.csv";
  std::string filename = ss.str();

  if (numPE != 0 && detected && (tof > earlyCutoff) && (tof < lateCutoff)){
      run->csvfi.open(filename,std::ios_base::app);
      run->csvfi << eventID << "," << numPE << "," << numPhotons << "," << eDep << "," << numElasticSensitive << "," << numInelastic << "," << nCapture << "," << tof << "," << detector << "," << nucleusRecoilEnergy << "," << scatteredNotSensitive << "," << numInelasticSensitive << "," << std::endl;
      run->csvfi.close();
  }

  //auto payload = buf_hits.str();
  //if (!payload.empty()){
      //G4cout << "bufhits size : " << payload.size() << G4endl;
      //if (detected) std::cout << "also detected" << std::endl;
      //if ((tof > earlyCutoff) && (tof < lateCutoff)) std::cout << "also tof" << std::endl;
  //}
  //if (detected && (tof > earlyCutoff) && (tof < lateCutoff) && !payload.empty()) {
   // OpenPerThreadFilesOnce();
   // AppendHitsPerThread(payload);
    //G4cout << "bufhits size : " << payload.size() << G4endl;


    //std::cout << "after buf hits" << std::endl;
    //std::ofstream f("photHits.csv", std::ios::app);
    //f << buf_hits.str();

    //run->AppendPhotHitsCSV(buf_hits.str());
    //buf_hits.str("");
    //buf_hits.clear();

    //std::ostringstream evtline;
    //evtline << eventID << ',' << detector << ',' << tOne << ',' << numElasticSensitive << ',' << numInelasticSensitive << ',' << ExtScatter << "," << numInelastic << "," << good; 
   // AppendStatsPerThread(evtline.str());
   
    //std::cout << "after appendstatsperthread" << std::endl;

    //run->AppendEvtStatsCSV(evtline.str());
   
    //run->csvfi.open("evtStats.csv",std::ios_base::app);
    //run->csvfi << eventID << "," << detector << "," << tOne << "," << numElasticSensitive << "," << numInelasticSensitive << "," << ExtScatter << "," << numInelastic << "," << good << "," <<std::endl;
    //run->csvfi.close();
  //}
 // std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout<<std::endl;

  //G4int numCompletedEvents = run->completedEvents.size();
  if (eventID % 1000 == 0){
      std::cout << "event " << eventID << std::endl;
  }


 run->AddEnergy(fEnergyDeposit);
 run->AddTrakLenCharg(fTrakLenCharged);
 run->AddTrakLenNeutr(fTrakLenNeutral);

 run->CountStepsCharg(fNbStepsCharged);
 run->CountStepsNeutr(fNbStepsNeutral);

 run->CountTransmit (fTransmitFlag);
 run->CountReflect  (fReflectFlag);
 
// if (fEnergyDeposit > 0.){
 //   G4AnalysisManager::Instance()->FillH1(1,fEnergyDeposit);
  //  G4AnalysisManager::Instance()->FillNtupleDColumn(1,fEnergyDeposit);
   // G4AnalysisManager::Instance()->AddNtupleRow();
 //   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

