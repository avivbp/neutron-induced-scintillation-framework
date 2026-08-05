#include "ParticleBookkeeping.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <sstream>

namespace {

std::mutex particleBookkeepingCsvMutex;

std::string LowerCopy(const std::string& value)
{
  std::string result = value;
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}

}  // namespace

ParticleClass ClassifyParticle(
    const std::string& particleName, const std::string& particleType,
    int atomicNumber, int atomicMass, int parentId)
{
  const auto name = LowerCopy(particleName);
  const auto type = LowerCopy(particleType);
  if (name == "neutron") {
    return parentId == 0 ? ParticleClass::PrimaryNeutron
                         : ParticleClass::SecondaryNeutron;
  }
  if (name == "gamma") {
    return ParticleClass::Gamma;
  }
  if (name == "e-" || name == "e+") {
    return ParticleClass::ElectronPositron;
  }
  if (name == "proton") {
    return ParticleClass::Proton;
  }
  if (name == "alpha" || (atomicNumber == 2 && atomicMass == 4)) {
    return ParticleClass::Alpha;
  }
  if (type == "nucleus" || atomicMass > 1) {
    return ParticleClass::NuclearRecoil;
  }
  return ParticleClass::Other;
}

const char* ParticleClassName(ParticleClass particleClass)
{
  switch (particleClass) {
    case ParticleClass::PrimaryNeutron:
      return "primary_neutron";
    case ParticleClass::SecondaryNeutron:
      return "secondary_neutron";
    case ParticleClass::Gamma:
      return "gamma";
    case ParticleClass::ElectronPositron:
      return "electron_positron";
    case ParticleClass::Proton:
      return "proton";
    case ParticleClass::Alpha:
      return "alpha";
    case ParticleClass::NuclearRecoil:
      return "nuclear_recoil";
    case ParticleClass::Other:
      return "other";
    case ParticleClass::Count:
      break;
  }
  return "other";
}

bool UsesIonScintillationYieldScale(ParticleClass particleClass)
{
  return particleClass == ParticleClass::PrimaryNeutron ||
         particleClass == ParticleClass::SecondaryNeutron ||
         particleClass == ParticleClass::Alpha ||
         particleClass == ParticleClass::NuclearRecoil;
}

double EventParticleBookkeeping::TotalEnergyDepositMeV() const
{
  return std::accumulate(energyDepositMeV.begin(), energyDepositMeV.end(),
                         0.0);
}

std::uint64_t EventParticleBookkeeping::TotalScintillationPhotons() const
{
  return std::accumulate(scintillationPhotons.begin(),
                         scintillationPhotons.end(), std::uint64_t{0});
}

std::string ParticleBookkeepingCsvHeader()
{
  std::ostringstream header;
  header << "event_id,total_active_lar_energy_deposit_MeV,"
         << "total_scintillation_photons";
  for (std::size_t index = 0; index < kParticleClassCount; ++index) {
    const auto particleClass = static_cast<ParticleClass>(index);
    header << ',' << ParticleClassName(particleClass)
           << "_energy_deposit_MeV";
    header << ',' << ParticleClassName(particleClass)
           << "_scintillation_photons";
  }
  header << '\n';
  return header.str();
}

std::string ParticleBookkeepingCsvRow(
    const EventParticleBookkeeping& bookkeeping)
{
  std::ostringstream row;
  row << std::setprecision(12)
      << bookkeeping.eventId << ','
      << bookkeeping.TotalEnergyDepositMeV() << ','
      << bookkeeping.TotalScintillationPhotons();
  for (std::size_t index = 0; index < kParticleClassCount; ++index) {
    row << ',' << bookkeeping.energyDepositMeV[index]
        << ',' << bookkeeping.scintillationPhotons[index];
  }
  row << '\n';
  return row.str();
}

std::string ParticleBookkeepingCsvFilename(const std::string& runLabel)
{
  if (runLabel.empty()) {
    return "particle_class_summary.csv";
  }

  std::string safeLabel = runLabel;
  std::transform(safeLabel.begin(), safeLabel.end(), safeLabel.begin(),
                 [](unsigned char c) {
                   return std::isalnum(c) || c == '_' || c == '-' ? c : '_';
                 });
  return "particle_class_summary_" + safeLabel + ".csv";
}

bool InitializeParticleBookkeepingCsv(const std::string& filename)
{
  std::lock_guard<std::mutex> lock(particleBookkeepingCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::trunc);
  if (!output) {
    return false;
  }
  output << ParticleBookkeepingCsvHeader();
  return output.good();
}

bool AppendParticleBookkeepingCsv(
    const EventParticleBookkeeping& bookkeeping, const std::string& filename)
{
  const auto row = ParticleBookkeepingCsvRow(bookkeeping);
  std::lock_guard<std::mutex> lock(particleBookkeepingCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::app);
  if (!output) {
    return false;
  }
  output << row;
  return output.good();
}
