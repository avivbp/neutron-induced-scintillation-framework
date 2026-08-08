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
std::mutex duneEventResponseCsvMutex;

std::string LowerCopy(const std::string& value)
{
  std::string result = value;
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}

template <typename Array>
typename Array::value_type SumByOrigin(
    const Array& values, bool (*predicate)(ParticleClass))
{
  typename Array::value_type total{};
  for (std::size_t index = 0; index < kParticleClassCount; ++index) {
    const auto particleClass = static_cast<ParticleClass>(index);
    if (predicate(particleClass)) {
      total += values[index];
    }
  }
  return total;
}

std::string SafeRunLabel(const std::string& runLabel)
{
  std::string safeLabel = runLabel;
  std::transform(safeLabel.begin(), safeLabel.end(), safeLabel.begin(),
                 [](unsigned char c) {
                   return std::isalnum(c) || c == '_' || c == '-' ? c : '_';
                 });
  return safeLabel;
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

bool IsElectronicRecoilOrigin(ParticleClass particleClass)
{
  return particleClass == ParticleClass::Gamma ||
         particleClass == ParticleClass::ElectronPositron;
}

bool IsNuclearRecoilOrigin(ParticleClass particleClass)
{
  // Match the particle classes whose scintillation is treated with the
  // configured ion-yield scale. Proton light retains the historical policy
  // and is therefore reported in the explicit "other" component.
  return UsesIonScintillationYieldScale(particleClass);
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

std::uint64_t EventParticleBookkeeping::TotalDetectedPhotoelectrons() const
{
  return std::accumulate(detectedPhotoelectrons.begin(),
                         detectedPhotoelectrons.end(), std::uint64_t{0});
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
    header << ',' << ParticleClassName(particleClass)
           << "_detected_photoelectrons";
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
        << ',' << bookkeeping.scintillationPhotons[index]
        << ',' << bookkeeping.detectedPhotoelectrons[index];
  }
  row << '\n';
  return row.str();
}

std::string ParticleBookkeepingCsvFilename(const std::string& runLabel)
{
  if (runLabel.empty()) {
    return "particle_class_summary.csv";
  }

  return "particle_class_summary_" + SafeRunLabel(runLabel) + ".csv";
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

std::string DuneEventResponseCsvHeader()
{
  return "event_id,total_num_pe,electronic_recoil_num_pe,"
         "nuclear_recoil_num_pe,other_num_pe,"
         "total_scintillation_photons,"
         "electronic_recoil_scintillation_photons,"
         "nuclear_recoil_scintillation_photons,"
         "other_scintillation_photons,"
         "total_active_lar_energy_deposit_MeV,"
         "electronic_recoil_energy_deposit_MeV,"
         "nuclear_recoil_energy_deposit_MeV,"
         "other_energy_deposit_MeV,fiducial_elastic_count,"
         "total_inelastic_count,fiducial_inelastic_count,n_capture\n";
}

std::string DuneEventResponseCsvRow(
    const EventParticleBookkeeping& bookkeeping, int fiducialElasticCount,
    int totalInelasticCount, int fiducialInelasticCount, bool neutronCapture)
{
  const auto erPE = SumByOrigin(bookkeeping.detectedPhotoelectrons,
                                IsElectronicRecoilOrigin);
  const auto nrPE = SumByOrigin(bookkeeping.detectedPhotoelectrons,
                                IsNuclearRecoilOrigin);
  const auto otherPE = SumByOrigin(
      bookkeeping.detectedPhotoelectrons,
      [](ParticleClass particleClass) {
        return !IsElectronicRecoilOrigin(particleClass) &&
               !IsNuclearRecoilOrigin(particleClass);
      });
  const auto erPhotons = SumByOrigin(bookkeeping.scintillationPhotons,
                                     IsElectronicRecoilOrigin);
  const auto nrPhotons = SumByOrigin(bookkeeping.scintillationPhotons,
                                     IsNuclearRecoilOrigin);
  const auto otherPhotons = SumByOrigin(
      bookkeeping.scintillationPhotons, [](ParticleClass particleClass) {
        return !IsElectronicRecoilOrigin(particleClass) &&
               !IsNuclearRecoilOrigin(particleClass);
      });
  const auto erEnergy = SumByOrigin(bookkeeping.energyDepositMeV,
                                    IsElectronicRecoilOrigin);
  const auto nrEnergy = SumByOrigin(bookkeeping.energyDepositMeV,
                                    IsNuclearRecoilOrigin);
  const auto otherEnergy = SumByOrigin(
      bookkeeping.energyDepositMeV, [](ParticleClass particleClass) {
        return !IsElectronicRecoilOrigin(particleClass) &&
               !IsNuclearRecoilOrigin(particleClass);
      });

  std::ostringstream row;
  row << std::setprecision(12) << bookkeeping.eventId << ','
      << bookkeeping.TotalDetectedPhotoelectrons() << ',' << erPE << ','
      << nrPE << ',' << otherPE << ','
      << bookkeeping.TotalScintillationPhotons() << ',' << erPhotons << ','
      << nrPhotons << ',' << otherPhotons
      << ',' << bookkeeping.TotalEnergyDepositMeV() << ',' << erEnergy << ','
      << nrEnergy << ',' << otherEnergy << ','
      << fiducialElasticCount << ',' << totalInelasticCount << ','
      << fiducialInelasticCount << ',' << (neutronCapture ? 1 : 0) << '\n';
  return row.str();
}

std::string DuneEventResponseCsvFilename(const std::string& runLabel)
{
  if (runLabel.empty()) {
    return "dune_event_response.csv";
  }
  return "dune_event_response_" + SafeRunLabel(runLabel) + ".csv";
}

bool InitializeDuneEventResponseCsv(const std::string& filename)
{
  std::lock_guard<std::mutex> lock(duneEventResponseCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::trunc);
  if (!output) {
    return false;
  }
  output << DuneEventResponseCsvHeader();
  return output.good();
}

bool AppendDuneEventResponseCsv(
    const EventParticleBookkeeping& bookkeeping, int fiducialElasticCount,
    int totalInelasticCount, int fiducialInelasticCount, bool neutronCapture,
    const std::string& filename)
{
  const auto row = DuneEventResponseCsvRow(
      bookkeeping, fiducialElasticCount, totalInelasticCount,
      fiducialInelasticCount, neutronCapture);
  std::lock_guard<std::mutex> lock(duneEventResponseCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::app);
  if (!output) {
    return false;
  }
  output << row;
  return output.good();
}
