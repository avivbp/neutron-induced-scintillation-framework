#ifndef ParticleBookkeeping_h
#define ParticleBookkeeping_h 1

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

enum class ParticleClass : std::size_t {
  PrimaryNeutron = 0,
  SecondaryNeutron,
  Gamma,
  ElectronPositron,
  Proton,
  Alpha,
  NuclearRecoil,
  Other,
  Count,
};

constexpr std::size_t kParticleClassCount =
    static_cast<std::size_t>(ParticleClass::Count);

ParticleClass ClassifyParticle(
    const std::string& particleName, const std::string& particleType,
    int atomicNumber, int atomicMass, int parentId);
const char* ParticleClassName(ParticleClass particleClass);
bool UsesIonScintillationYieldScale(ParticleClass particleClass);
bool IsElectronicRecoilOrigin(ParticleClass particleClass);
bool IsNuclearRecoilOrigin(ParticleClass particleClass);

struct EventParticleBookkeeping {
  int eventId = 0;
  std::array<double, kParticleClassCount> energyDepositMeV{};
  std::array<std::uint64_t, kParticleClassCount> scintillationPhotons{};
  std::array<std::uint64_t, kParticleClassCount> detectedPhotoelectrons{};

  double TotalEnergyDepositMeV() const;
  std::uint64_t TotalScintillationPhotons() const;
  std::uint64_t TotalDetectedPhotoelectrons() const;
};

std::string ParticleBookkeepingCsvHeader();
std::string ParticleBookkeepingCsvRow(
    const EventParticleBookkeeping& bookkeeping);
std::string ParticleBookkeepingCsvFilename(const std::string& runLabel);

bool InitializeParticleBookkeepingCsv(
    const std::string& filename = "particle_class_summary.csv");
bool AppendParticleBookkeepingCsv(
    const EventParticleBookkeeping& bookkeeping,
    const std::string& filename = "particle_class_summary.csv");

std::string DuneEventResponseCsvHeader();
std::string DuneEventResponseCsvRow(
    const EventParticleBookkeeping& bookkeeping, int fiducialElasticCount,
    int totalInelasticCount, int fiducialInelasticCount, bool neutronCapture);
std::string DuneEventResponseCsvFilename(const std::string& runLabel);

bool InitializeDuneEventResponseCsv(
    const std::string& filename = "dune_event_response.csv");
bool AppendDuneEventResponseCsv(
    const EventParticleBookkeeping& bookkeeping, int fiducialElasticCount,
    int totalInelasticCount, int fiducialInelasticCount, bool neutronCapture,
    const std::string& filename = "dune_event_response.csv");

#endif
