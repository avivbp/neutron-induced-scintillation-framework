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

struct EventParticleBookkeeping {
  int eventId = 0;
  std::array<double, kParticleClassCount> energyDepositMeV{};
  std::array<std::uint64_t, kParticleClassCount> scintillationPhotons{};

  double TotalEnergyDepositMeV() const;
  std::uint64_t TotalScintillationPhotons() const;
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

#endif
