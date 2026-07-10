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
/// \file electromagnetic/TestEm5/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4EventManager.hh"
#include "G4AtomicShells.hh"
#include "G4AutoLock.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Threading.hh"
#include "G4ProcessManager.hh"
#include "G4LogicalBorderSurface.hh"

namespace {
  G4ThreadLocal G4OpBoundaryProcess* tlBoundaryProc = nullptr;
}

static G4OpBoundaryProcess* GetBoundaryProc()
{
  if (tlBoundaryProc) return tlBoundaryProc;

  auto* pm = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  if (!pm) return nullptr;

  const auto* pv = pm->GetProcessList();
  for (size_t i = 0; i < pv->size(); ++i) {
    auto* p = (*pv)[i];
    if (p->GetProcessName() == "OpBoundary") {
      tlBoundaryProc = static_cast<G4OpBoundaryProcess*>(p);
      break;
    }
  }
  return tlBoundaryProc;
}
//namespace { G4Mutex wlsPrintMtx = G4MUTEX_INITIALIZER; }
//namespace { G4Mutex bndMtx = G4MUTEX_INITIALIZER; }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DET,
                               EventAction* EA)
:G4UserSteppingAction(),fDetector(DET), fEventAction(EA)
{ fTPBMat = fDetector->GetTPBMaterial();
  fInnerCellLV = fDetector->GetInnerCellLV();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int parentID = aStep->GetTrack()->GetParentID();
  G4double gammaEDep = 0;
  G4bool NP = 0;
  G4bool NNprime = 0;
  G4bool NNPrimeP = 0;
  G4bool N2N = 0;
  G4bool N3N = 0;
  G4bool nPrimeThisStep = false;
  G4double emittedNeutronAngle = 0.0;
 // std::string positron("e+");
  std::string proton("proton");
  std::string neutron("neutron");
 // std::string gamma("gamma");
 // std::string phot("phot");
  std::string trans("Transportation");
  std::string elastic("hadElastic");
  std::string inelastic("neutronInelastic");
  std::string capture("nCapture");
  std::string scintillator("liquidScintillator");
  std::string sensitiveLAr("Absorber");
  std::string Absorber("biggerAbsorber");


  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4bool verbose2 = false;
  G4bool verbose1 = true;
  if (verbose2) {
      printEventStats(aStep,10);
  }
  const G4String & processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  const G4String & incomingParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  auto* def = aStep->GetTrack()->GetParticleDefinition();
  bool isNeutron = (def == G4Neutron::Definition());
  bool isProton  = (def == G4Proton::Definition());
  const G4String & volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
  G4int size  = (int) (secondaries->size());

  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
            aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  G4double pre_Ekin = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
  G4double post_Ekin = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
  G4double eDiff = pre_Ekin - post_Ekin;
  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
  G4double postStepTime = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
  G4double preStepTime = aStep->GetPreStepPoint()->GetGlobalTime()/CLHEP::ns;

  auto pre = aStep->GetPreStepPoint();
  auto post = aStep->GetPostStepPoint();
  auto prePV  = pre->GetPhysicalVolume();
  auto postPV = post->GetPhysicalVolume();
  if (!prePV || !postPV) return;
  auto preLV  = prePV->GetLogicalVolume();
  auto postLV = postPV->GetLogicalVolume();
  auto preName = prePV->GetName();
  auto postName = postPV->GetName();
  const G4bool isNeutronDetector = fDetector->IsNeutronDetectorName(volumeName);
  G4Track* track = aStep->GetTrack();
  G4int stepNo = track->GetCurrentStepNumber();

  if (verbose2) {
      std::cout << "after stepNo" << std::endl;
  }

  const G4ThreeVector & initMomentum = aStep->GetPreStepPoint()->GetMomentumDirection();
  const G4ThreeVector & finMomentum = aStep->GetPostStepPoint()->GetMomentumDirection();

  // Parameters (tweak as you like)
  constexpr int   kPosScale_cm_to_mm = 10;     // 1 digit after decimal in cm
  constexpr int   kTimeScale_10ps    = 100;    // ns → 10 ps ticks
  constexpr double kTwoPi = 2.0*CLHEP::pi;
  const G4VProcess* creator = track->GetCreatorProcess();

  // --- Stepping “fuses” (add after you define pre/post/prePV/postPV/etc.) ---
static thread_local G4int   lastEvt   = -1;
static thread_local G4int   lastTrack = -1;
static thread_local int     zeroLenStreak = 0;
static thread_local int     pingPongStreak = 0;
static thread_local const G4LogicalVolume* prevA = nullptr;
static thread_local const G4LogicalVolume* prevB = nullptr;
static thread_local unsigned long eventOptStepBudget = 0;

// Reset per-event/per-track state
const G4int evtId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
if (evtId != lastEvt) {
  lastEvt = evtId;
  eventOptStepBudget = 0;
}
if (trackID != lastTrack) {
  lastTrack = trackID;
  zeroLenStreak = 0;
  pingPongStreak = 0;
  prevA = prevB = nullptr;
}

const auto stepLen = track->GetStepLength();
const auto stepStat = post->GetStepStatus();

// A) Zero-length boundary streak (tunable thresholds)
const G4double epsLen = 1e-6*CLHEP::mm;                // ~1 nm
const int maxZeroStreak = 100;                         // e.g. 100 consecutive
if (stepStat == fGeomBoundary) {
  if (stepLen < epsLen) ++zeroLenStreak; else zeroLenStreak = 0;
  if (zeroLenStreak >= maxZeroStreak) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }
}

  if (verbose2) {
      std::cout << "after A)" << std::endl;
  }
// B) Ping-pong A<->B detection at boundary (microscopic steps)
const int maxPingPong = 200;                           // e.g. 400 flips
if (stepStat == fGeomBoundary) {
  const auto* A = prePV->GetLogicalVolume();
  const auto* B = postPV->GetLogicalVolume();
  if (prevA == B && prevB == A) ++pingPongStreak;
  else                          pingPongStreak = 0;
  prevA = A; prevB = B;
  if (pingPongStreak >= maxPingPong) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }
}
if (verbose2){
    std::cout << "after B)" << std::endl;
}

  // ---- BEGIN: universal guards (paste near top of UserSteppingAction) ----
const bool isOpt = (particleDef == opticalphoton);

if (track->GetGlobalTime() > 10000.0 * CLHEP::ns) {

    track->SetTrackStatus(fStopAndKill);
    return;
}
// 0) Hard cap for non-opticals too (you already have this block; keep it here early)
if (!isOpt && stepNo > 2500) {
  //G4cout << "particle : " << particleDef->GetParticleName()
    //     << " exceeded 2500 steps, pre->post " << prePV->GetName() << "->" << postPV->GetName()
      //   << ", stepLen " << track->GetStepLength()/CLHEP::um << " um, status " << post->GetStepStatus()
        // << G4endl;
  track->SetTrackStatus(fStopAndKill);
  return;
}

if (verbose2){
std::cout << "after 0" << std::endl;
}
// 1) Kill very old tracks (optical + non-optical)
if (track->GetGlobalTime() > 5.0*CLHEP::microsecond) {
  // If not a WLS optical, tell me who it is
  if (!isOpt || (isOpt && track->GetCreatorProcess() &&
                 track->GetCreatorProcess()->GetProcessName() != "OpWLS")) {
      const auto* postLV = post->GetPhysicalVolume()->GetLogicalVolume();
      const G4String postName = post->GetPhysicalVolume()->GetName();

      const bool inAirOrWorld =
          (postName.contains("World") || postName.contains("air") || postName.contains("Air"));

      if (inAirOrWorld) {
          // gentle but effective: old OR very slow neutrons in air don’t matter for LAr scint
          if (track->GetGlobalTime() > 3000.*CLHEP::ns) { track->SetTrackStatus(fStopAndKill); return; }

          const auto ekin = post->GetKineticEnergy();
          if (ekin < 1.*CLHEP::keV) { track->SetTrackStatus(fStopAndKill); return; }
      }
   // const G4VProcess* cr = track->GetCreatorProcess();
  //  G4cout << "[>5us] name=" << particleDef->GetParticleName()
    //       << " creator=" << (cr ? cr->GetProcessName() : "(primary/no-creator)")
      //     << " pre->post " << prePV->GetName() << "->" << postPV->GetName()
        //   << G4endl;
  }
}
if (verbose2){
std::cout << "after 1)" << std::endl;
}
// 2) Extra optical fuses
if (isOpt) {
  // 2a) Per-track zero-length churn guard
  static thread_local G4int lastTrkId = -1;
  static thread_local int   zeroStreak = 0;
  if (track->GetTrackID() != lastTrkId) { zeroStreak = 0; lastTrkId = track->GetTrackID(); }

  const bool atBoundary = (post->GetStepStatus() == fGeomBoundary);
  const G4double sl = aStep->GetStepLength();
  if (atBoundary && sl < 1e-9*CLHEP::mm) {
    if (++zeroStreak > 200) {
        G4cout << "track reached zero length streak of 200" << G4endl;
        track->SetTrackStatus(fStopAndKill); return; }
  } else {
    zeroStreak = 0;
  }
if (verbose2){
std::cout << "after 2a)" << std::endl;
}
  // 2b) Per-track boundary count + cumulative path length
  struct TrkState { int bounces=0; double accLen=0.0; };
  static thread_local int      ts_lastId = -1;
  static thread_local TrkState ts;
  if (track->GetTrackID() != ts_lastId) { ts = TrkState{}; ts_lastId = track->GetTrackID(); }

  if (atBoundary) { ++ts.bounces; }
  ts.accLen += sl;

  if (ts.bounces > 250)                     { track->SetTrackStatus(fStopAndKill); return; }
  if (ts.accLen  > 5.0*CLHEP::m)             { track->SetTrackStatus(fStopAndKill); return; }

  // 2c) Your existing optical step/time caps (keep them, but they can now be lower)
  if (stepNo > 2500)                         { track->SetTrackStatus(fStopAndKill); return; }
}
if (verbose2){
std::cout << "after 2B)" << std::endl;
}
// ---- END: universal guards ----

 // if (eventID == 0) {
  //const int tid = G4Threading::G4GetThreadId();
  //G4cout << "eventID= " << eventID << " , [hb] tid=" << tid <<
    //" trk=" << track->GetTrackID()
    // << " pre=" << prePV->GetName()
     // << " post="<< postPV->GetName()
     // << " len=" << aStep->GetStepLength()/CLHEP::um << " um"
      //<< " status=" << post->GetStepStatus()
      // << G4endl;
//}

  // tracking each scintillation photon individually
  G4int checkByHand = 1;
  if (checkByHand){
  if(particleDef == opticalphoton){

            if(creator && creator->GetProcessName() == "OpWLS" && post->GetStepStatus() == fGeomBoundary){

            auto boundary = GetBoundaryProc();
            if (boundary) {
                G4OpBoundaryProcessStatus status = boundary->GetStatus();
                G4bool DETECTION = (status == Detection);
                G4bool ABSORPTION = (status == Absorption);
                G4bool REFLECTION = ( status == FresnelReflection || status == TotalInternalReflection || status == LambertianReflection || status == SpikeReflection );
                const auto* border = G4LogicalBorderSurface::GetSurface(prePV, postPV);
                const G4SurfaceProperty* surfProp = border ? border->GetSurfaceProperty() : nullptr;
                const auto* optSurf = dynamic_cast<const G4OpticalSurface*>(surfProp);

                if (optSurf && ((optSurf->GetName().rfind("SiPM",0) == 0) || (optSurf->GetName().rfind("PMT",0) == 0))) {
                    if (DETECTION){
                        fEventAction->numPE += 1;
                        //std::cout << "optical surface : " << optSurf->GetName() << std::endl;
                        //std::cout << "num PE = " << fEventAction->numPE << std::endl;
                        track->SetTrackStatus(fStopAndKill);
                        // photon detected by surface
                    }
               
                    else if (ABSORPTION) {
                        //std::cout << "absorbed" << std::endl;
                    }
                    else if (REFLECTION){
                        //std::cout << "reflected at " << preName << " - " << postName << " boundary " << std::endl;
                    }
                }
            }              
            G4double t = post->GetGlobalTime()/CLHEP::ns;
            if (t < fEventAction->tZeroEx) {fEventAction->tZeroEx = t;}
        }

        //if(creator && creator->GetProcessName() == "OpWLS" ){

          //if (verbose1) std::cout << "created by OpWLS" << std::endl;

          //double t = post->GetGlobalTime()/CLHEP::ns; // Time of flight at boundary

	  //if (post->GetStepStatus() == fGeomBoundary) {
              //if (verbose1) std::cout << "at fGeomBoundary" << std::endl;

            // Decide based on QE and reflectivity (depending on your setup)
            //  if (G4UniformRand() < 0.25 && t < fEventAction->tOneEx) {  // Apply some QE, e.g. 25%
            //    fEventAction->tOneEx = t;
            //  }

            //G4OpBoundaryProcess* boundary = nullptr;

            //auto procMan = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
            //if (procMan) {
             // const G4ProcessVector* pv = procMan->GetProcessList();
              //for (size_t i = 0; i < pv->size(); ++i) {
              //  auto p = (*pv)[i];
              //  if (p->GetProcessName() == "OpBoundary") {
              //    boundary = static_cast<G4OpBoundaryProcess*>(p);
              //    break;
               // }
             // }
           // }
           // if (boundary) {
           //     G4OpBoundaryProcessStatus status = boundary->GetStatus();

            //    if (status == Detection) {
             //     printEventStats(aStep,10);
              //    fEventAction->numPE += 1;
               //   track->SetTrackStatus(fStopAndKill);
                  // photon detected by surface
               // }
               // else if (status == Absorption) {
                //  std::cout << "absorbed" << std::endl;
                  // absorbed
               // }
               // else if (status == FresnelReflection ||
                //         status == TotalInternalReflection ||
                //         status == LambertianReflection ||
                 //        status == SpikeReflection) {
                  //       std::cout << "reflected" << std::endl;
                  // reflected
               // }
           // }              


            //const bool fromTPB = (preName.rfind("TPB", 0) == 0);   // any TPB_* PV
            //if (fromTPB) G4cout << "from tpb material name = " << preLV->GetMaterial()->GetName() << G4endl;
            //const bool toInner = (postName == "innerCell");
            //if (toInner) G4cout << "to inner post logical volume name = " << postLV->GetName() << G4endl;
            // accept only crossings from tpb to the inner cell
           // if (fromTPB && toInner){ 
                //if (verbose1) std::cout << "from TPB to inner" << std::endl;
                //G4cout << "from tpb to inner" << G4endl;

                //fEventAction->numPhotons += 1;
                // classify surface based on which TPB piece we just left (prePV)
             //   char surfC = 'O';
               // const G4String& preName = prePV->GetName();
               // if (preName.find("Side") != std::string::npos)     surfC = 'S';
               // else if (preName.find("Top") != std::string::npos)      surfC = 'T';
               // else if (preName.find("Bot") != std::string::npos) surfC = 'B';
                // else 'O' (other) by default

                //const auto p  = post->GetPosition();
                //if (verbose1) std::cout << "preName = " << preName << ", t = " << t << " ns" << std::endl;
	        //if (surfC != 'O') {
                 // if (t < fEventAction->tZeroEx){
                   //   fEventAction->tZeroEx = t;
                      //std::cout << "tZeroEx = " << t <<  " ns" << std::endl;
                  //}
               // }

                // Helpers
               // auto wrap = [](double a, double w) {
               //   double r = std::fmod(a, w);
               //   if (r < 0) r += w;
               //   return r;
               // };
                
                //int u_int=0, v_int=0;   // quantized positions
                
                // Bases: (u,v)=(x_cm,z_cm) → mm ticks
              //  if (surfC=='T' || surfC=='B') {
              //    u_int = static_cast<int>(std::llround( (p.x()/CLHEP::cm) * kPosScale_cm_to_mm ));
               //   v_int = static_cast<int>(std::llround( (p.z()/CLHEP::cm) * kPosScale_cm_to_mm ));
               // }
                // Side: (u,v)=(phi_rad,y_cm) → tile index & mm ticks
               // else if (surfC=='S') {
                  //tile count for R=5 cm, tile=0.6 cm
                 // constexpr double R_cm = 5.0;
                 // constexpr double tile_cm = 0.6;
                 // const int nphi = static_cast<int>(std::floor((kTwoPi*R_cm)/tile_cm)); // 52
                 // const double dphi = kTwoPi / nphi;
                
                  //double phi = std::atan2(p.z(), p.x());      // [-pi, pi]
                 // phi = wrap(phi, kTwoPi);                    // [0, 2pi)
                 // int k = static_cast<int>(std::floor(phi / dphi));
                 // if (k >= nphi) k = nphi - 1;                // guard
                
                 // u_int = k;  // store tile index instead of φ (massive savings)
                 // v_int = static_cast<int>(std::llround( (p.y()/CLHEP::cm) * kPosScale_cm_to_mm ));

                  //u_int = static_cast<int>(std::llround( (p.x()/CLHEP::cm) * kPosScale_cm_to_mm ));
                  //v_int = static_cast<int>(std::llround( (p.z()/CLHEP::cm) * kPosScale_cm_to_mm ));

               // }
                // Time in 10 ps ticks
                //int t_10ps = static_cast<int>(std::llround(t * kTimeScale_10ps ));
                
                //if (verbose1) std::cout << "before log WLS" << std::endl;
                //if (fEventAction->ShouldLogWLS(trackID, surfC, u_int, v_int, t_10ps)) {
                    //if (verbose1) std::cout << "should log WLS" << std::endl;
                    //G4cout << "should log WLS" << G4endl;
                    //if (track->GetCurrentStepNumber() > 2500) G4cout << "when logging wls, step number is : " << track->GetCurrentStepNumber() << G4endl;
                  //  fEventAction->BufferHitRow({std::to_string(eventID),std::to_string(trackID),std::string(1,surfC),std::to_string(u_int),std::to_string(v_int),std::to_string(t_10ps)});
                //}
 
	  //}    
     // }
     // }
  }
  }
if (verbose2){
std::cout << "after optical)" << std::endl;
}
  G4int hmm2 = 1;
  if (hmm2 && isNeutron){
      //printEventStats(aStep,10);
      if (inelastic.compare(processName) == 0){ 
          fEventAction->numInelastic += 1;
          //fEventAction->aborted = true;
          G4EventManager::GetEventManager()->AbortCurrentEvent();
          //run->addNumInelastic();
     
          G4int ns = 0;
          G4int ps = 0;
          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              if (neutron.compare(particleName) == 0){
                  ns += 1;
                  fEventAction->secondaryNeutron = true;
              } 
              else if (particleName == "proton"){
                  ps += 1;
              }
          }
          if (ns == 0 && ps == 1){
              NP = 1;
              fEventAction->NP = 1;
          }
          else if (ns == 1 && ps == 1){
              NNPrimeP = 1;
              fEventAction->NNPrimeP = 1;
          }
          else if (ns == 1 && ps == 0){
              NNprime = 1;
              fEventAction->NNPrime = 1;
          }
          else if (ns == 2 && ps == 0){
              N2N = 1;
              fEventAction->N2N = 1;
          }
          else if (ns == 3 && ps == 0){
              N3N = 1;
              fEventAction->N3N = 1;
          }
          else {
              //printEventStats(aStep,20);
              //G4EventManager::GetEventManager()->AbortCurrentEvent();
          }
      }
  }
  if (verbose2){
std::cout << "after hmm2)" << std::endl;
}
  G4int hmm = 1;
  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (hmm && isNeutron){

      G4VPhysicalVolume* postVolume = aStep->GetPostStepPoint()->GetPhysicalVolume();

      if (postVolume && postVolume->GetName() == "World" && volumeName == "air" && processName == "Transportation") {

          fEventAction->escapeEn = post_Ekin;
      }

      if (volumeName == "innerLayer" && trans.compare(processName) != 0){
          //std::cout << "scattered off of aluminum layer with scattering angle = " << calcAngle(initMomentum,finMomentum) << std::endl;
          fEventAction->innerLayerScatter += 1;

          if (NNprime){
              fEventAction->innerLayerEDep += pre_Ekin - (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
          }

          else{
              fEventAction->innerLayerEDep += eDiff;
          }
      }

      if (volumeName == "outerCell" && trans.compare(processName) != 0){
          fEventAction->ExtScatter += 1;
 
          // before interaction in fiducial volume
          if (fEventAction->tZero == 10000){
              fEventAction->extScatterBefore += 1;
          }
          else{
              fEventAction->extScatterAfter += 1;
          }
          
          if (NNprime){
              fEventAction->outerCellScatterAngle = calcAngle(initMomentum,(*secondaries)[0]->GetMomentumDirection());
              fEventAction->outerCellEDep += pre_Ekin - (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
          }

          else{
              fEventAction->outerCellScatterAngle = calcAngle(initMomentum,finMomentum);
              fEventAction->outerCellEDep += eDiff;
          }
      }

      if ((volumeName == "outerLayerOne" || volumeName == "outerLayerTwo") && trans.compare(processName) != 0){
          fEventAction->CryoScatter += 1;

          if (NNprime){
              fEventAction->CryoScatterAngle = calcAngle(initMomentum,(*secondaries)[0]->GetMomentumDirection());
              fEventAction->CryoEDep += pre_Ekin - (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
          }

          else{
              fEventAction->CryoScatterAngle = calcAngle(initMomentum,finMomentum);
              fEventAction->CryoEDep += eDiff;
          }
      }

      // calculate neutron angle relative to beam line as it goes into fiducial volume
      if (volumeName == "innerLayer" && postVolume->GetName() == "innerCell" && processName == "Transportation" && fEventAction->tZero == 10000){
          G4double incomingAngle = calcAngle(G4ThreeVector(1,0,0),initMomentum);
          fEventAction->incomingAng = incomingAngle;
      }

      if (volumeName == "innerCell" && trans.compare(processName) != 0 && fEventAction->tZero > postStepTime){
          //std::cout << "interaction in fiducial volume" << std::endl;

          fEventAction->fiducialIncomingEn = pre_Ekin;
          fEventAction->tZero = postStepTime;
          //std::cout << "tZero = " << postStepTime << " ns" << std::endl;
      }

      if (volumeName == "innerCell" && elastic.compare(processName) == 0){
          fEventAction->nucleusRecoilEnergy += pre_Ekin - post_Ekin;
          fEventAction->numElasticSensitive += 1;
          fEventAction->innerCellEDep += eDiff;
      }

      if (volumeName == "innerCell" && inelastic.compare(processName) == 0){
          //ONLY FOR TOF RUN, REMOVE LATER!!!!
          //G4EventManager::GetEventManager()->AbortCurrentEvent();
          fEventAction->numInelasticSensitive += 1;
          if (NNprime){
              fEventAction->innerCellEDep += pre_Ekin - (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
          }
      }

       if (processName == "nCapture"){
           fEventAction->nCapture = true;
       }


      // if a neutron has an interaction with the liquid scintillator
      if (!fEventAction->detected && trans.compare(processName) != 0 && isNeutronDetector){ 
          //std::cout << "interaction with Ax" << std::endl;
          //printEventStats(aStep,0);

          if (secondaries && !secondaries->empty()) {
              const G4Track* secondary = (*secondaries)[0];
              if (secondary->GetParticleDefinition()->GetParticleName().compare(proton) == 0) {
                  // printEventStats(aStep);
                  G4double en_trans = pre_Ekin - post_Ekin;
                  fEventAction->enDep += en_trans;
                  fEventAction->detector = volumeName;
                  fEventAction->detected = true;

                  // stop timer when neutron deposits energy in the organic scintillator
                  fEventAction->tOne = aStep->GetPostStepPoint()->GetGlobalTime() / CLHEP::ns;
                  //std::cout << "in event " << eventID << ", tOne = " << fEventAction->tOne << " ns" << std::endl;
              }
          }
      }

  }
if (verbose2){
std::cout << "after hmm)" << std::endl;
}

if (verbose2){
auto* proc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
G4String pname = proc ? proc->GetProcessName() : "None";
G4cout << "Step end: proc=" << pname
       << " vol=" << aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()
       << " step#=" << aStep->GetTrack()->GetCurrentStepNumber()
       << G4endl;

}
}
   G4bool SteppingAction::withinArapuca(const G4Step* aStep,G4double height){
       G4double xLen = 5; // ARAPUCA half side length in cm
       const G4ThreeVector & post = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
       G4double epsilon = 0.01;
       G4bool within = ((post.y() - height) < epsilon) && fabs(post.x()) < xLen && fabs(post.z()) < xLen;
    
       return within;
   }

   G4bool SteppingAction::withinPMTs(const G4Step* aStep,G4double height){
       G4double PMTSize = 2.05; // photocathode size in cm is 2.05x2.05 cm squared
       G4double totArea = pow(PMTSize,2) * 4; // 4 PMTs 
       G4double xLen = sqrt(totArea); // total effective area's side length

       const G4ThreeVector & post = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
       G4double epsilon = 0.01;
       G4bool within = ((fabs(post.y() - height) < epsilon) || (fabs(post.y() + height) < epsilon)) && fabs(post.x()) < xLen && fabs(post.z()) < xLen;

       return within;
   }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
