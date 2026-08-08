#include "ParticleBookkeeping.hh"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void Require(bool condition, const std::string& message)
{
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace

int main()
{
  Require(ClassifyParticle("neutron", "baryon", 0, 1, 0) ==
              ParticleClass::PrimaryNeutron,
          "primary neutron was misclassified");
  Require(ClassifyParticle("neutron", "baryon", 0, 1, 1) ==
              ParticleClass::SecondaryNeutron,
          "secondary neutron was misclassified");
  Require(ClassifyParticle("gamma", "gamma", 0, 0, 1) ==
              ParticleClass::Gamma,
          "gamma was misclassified");
  Require(ClassifyParticle("e+", "lepton", 0, 0, 1) ==
              ParticleClass::ElectronPositron,
          "positron was misclassified");
  Require(ClassifyParticle("proton", "baryon", 1, 1, 1) ==
              ParticleClass::Proton,
          "proton was misclassified");
  Require(ClassifyParticle("alpha", "nucleus", 2, 4, 1) ==
              ParticleClass::Alpha,
          "alpha was misclassified");
  Require(ClassifyParticle("Ar40", "nucleus", 18, 40, 1) ==
              ParticleClass::NuclearRecoil,
          "argon recoil was misclassified");
  Require(ClassifyParticle("mu-", "lepton", 0, 0, 1) ==
              ParticleClass::Other,
          "fallback particle was misclassified");

  Require(UsesIonScintillationYieldScale(ParticleClass::NuclearRecoil),
          "nuclear recoil did not use ion scintillation scale");
  Require(!UsesIonScintillationYieldScale(ParticleClass::Proton),
          "proton behavior changed from the existing thinning policy");
  Require(IsElectronicRecoilOrigin(ParticleClass::ElectronPositron),
          "electron was not grouped as electronic recoil");
  Require(IsNuclearRecoilOrigin(ParticleClass::NuclearRecoil),
          "nucleus was not grouped as nuclear recoil");
  Require(!IsNuclearRecoilOrigin(ParticleClass::Proton),
          "proton unexpectedly changed recoil grouping policy");

  EventParticleBookkeeping bookkeeping;
  bookkeeping.eventId = 7;
  bookkeeping.energyDepositMeV[static_cast<std::size_t>(
      ParticleClass::NuclearRecoil)] = 1.25;
  bookkeeping.energyDepositMeV[static_cast<std::size_t>(
      ParticleClass::Gamma)] = 0.75;
  bookkeeping.scintillationPhotons[static_cast<std::size_t>(
      ParticleClass::NuclearRecoil)] = 10;
  bookkeeping.scintillationPhotons[static_cast<std::size_t>(
      ParticleClass::Gamma)] = 20;
  bookkeeping.detectedPhotoelectrons[static_cast<std::size_t>(
      ParticleClass::NuclearRecoil)] = 3;
  bookkeeping.detectedPhotoelectrons[static_cast<std::size_t>(
      ParticleClass::ElectronPositron)] = 5;
  bookkeeping.detectedPhotoelectrons[static_cast<std::size_t>(
      ParticleClass::Proton)] = 2;

  Require(bookkeeping.TotalEnergyDepositMeV() == 2.0,
          "energy total does not match particle classes");
  Require(bookkeeping.TotalScintillationPhotons() == 30,
          "photon total does not match particle classes");
  Require(bookkeeping.TotalDetectedPhotoelectrons() == 10,
          "photoelectron total does not match particle classes");
  const auto row = ParticleBookkeepingCsvRow(bookkeeping);
  Require(row.find("7,2,30,") == 0, "CSV totals are incorrect");
  Require(ParticleBookkeepingCsvHeader().find(
              "nuclear_recoil_energy_deposit_MeV") != std::string::npos,
          "CSV header is missing nuclear-recoil energy");
  Require(ParticleBookkeepingCsvHeader().find(
              "nuclear_recoil_detected_photoelectrons") != std::string::npos,
          "CSV header is missing nuclear-recoil PE");
  Require(ParticleBookkeepingCsvFilename("C 01/test") ==
              "particle_class_summary_C_01_test.csv",
          "run label was not made filename-safe");

  const auto response = DuneEventResponseCsvRow(bookkeeping, 1, 2, 1, true);
  Require(response.find("7,10,5,3,2,30,20,10,0,2,0.75,1.25,0,1,2,1,1") ==
              0,
          "DUNE response did not preserve ER, NR, and other components");
  Require(DuneEventResponseCsvHeader().find("nuclear_recoil_num_pe") !=
              std::string::npos,
          "DUNE response header is missing nuclear-recoil PE");
  Require(DuneEventResponseCsvFilename("DUNE 00/test") ==
              "dune_event_response_DUNE_00_test.csv",
          "DUNE response run label was not made filename-safe");

  std::cout << "Particle bookkeeping tests passed\n";
  return EXIT_SUCCESS;
}
