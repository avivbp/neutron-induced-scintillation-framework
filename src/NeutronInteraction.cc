#include "NeutronInteraction.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>

namespace {

std::mutex interactionCsvMutex;

std::string LowerCopy(const std::string& value)
{
  std::string result = value;
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}

std::string EscapeCsv(const std::string& value)
{
  if (value.find_first_of(",\"\r\n") == std::string::npos) {
    return value;
  }

  std::string escaped = "\"";
  for (const char c : value) {
    if (c == '"') {
      escaped += "\"\"";
    } else {
      escaped += c;
    }
  }
  escaped += '"';
  return escaped;
}

std::string EncodeSecondaries(
    const std::vector<NeutronSecondaryTruth>& secondaries)
{
  std::ostringstream encoded;
  encoded << std::setprecision(12);
  for (std::size_t index = 0; index < secondaries.size(); ++index) {
    if (index != 0) {
      encoded << '|';
    }
    const auto& secondary = secondaries[index];
    encoded << secondary.particle << ':' << secondary.pdgCode << ':'
            << secondary.kineticEnergyKeV;
  }
  return encoded.str();
}

}  // namespace

NeutronInteractionChannel ClassifyNeutronInteraction(
    const std::string& processName, bool isHadronic, bool exitsWorld)
{
  const auto lowerName = LowerCopy(processName);
  if (exitsWorld &&
      (lowerName.empty() || lowerName == "none" ||
       lowerName == "transportation")) {
    return NeutronInteractionChannel::TransportExit;
  }
  if (lowerName == "hadelastic") {
    return NeutronInteractionChannel::Elastic;
  }
  if (lowerName == "neutroninelastic") {
    return NeutronInteractionChannel::Inelastic;
  }
  if (lowerName == "ncapture") {
    return NeutronInteractionChannel::Capture;
  }
  if (lowerName.find("fission") != std::string::npos) {
    return NeutronInteractionChannel::Fission;
  }
  if (isHadronic) {
    return NeutronInteractionChannel::OtherHadronic;
  }
  if (!lowerName.empty() && lowerName != "none" &&
      lowerName != "transportation") {
    return NeutronInteractionChannel::Unclassified;
  }
  return NeutronInteractionChannel::None;
}

const char* NeutronInteractionChannelName(NeutronInteractionChannel channel)
{
  switch (channel) {
    case NeutronInteractionChannel::Elastic:
      return "elastic";
    case NeutronInteractionChannel::Inelastic:
      return "inelastic";
    case NeutronInteractionChannel::Capture:
      return "capture";
    case NeutronInteractionChannel::Fission:
      return "fission";
    case NeutronInteractionChannel::OtherHadronic:
      return "other_hadronic";
    case NeutronInteractionChannel::TransportExit:
      return "transport_exit";
    case NeutronInteractionChannel::Unclassified:
      return "unclassified";
    case NeutronInteractionChannel::None:
      return "none";
  }
  return "unclassified";
}

std::string NeutronInteractionCsvHeader()
{
  return "event_id,interaction_index,track_id,parent_id,step_number,"
         "volume_name,volume_roles,process_name,channel,x_cm,y_cm,z_cm,"
         "time_ns,pre_kinetic_energy_keV,post_kinetic_energy_keV,"
         "local_energy_deposit_keV,secondary_count,secondaries\n";
}

std::string NeutronInteractionCsvRow(const NeutronInteractionRecord& record)
{
  std::ostringstream row;
  row << std::setprecision(12)
      << record.eventId << ','
      << record.interactionIndex << ','
      << record.trackId << ','
      << record.parentId << ','
      << record.stepNumber << ','
      << EscapeCsv(record.volumeName) << ','
      << EscapeCsv(record.volumeRoles) << ','
      << EscapeCsv(record.processName) << ','
      << NeutronInteractionChannelName(record.channel) << ','
      << record.xCm << ',' << record.yCm << ',' << record.zCm << ','
      << record.timeNs << ','
      << record.preKineticEnergyKeV << ','
      << record.postKineticEnergyKeV << ','
      << record.localEnergyDepositKeV << ','
      << record.secondaries.size() << ','
      << EscapeCsv(EncodeSecondaries(record.secondaries)) << '\n';
  return row.str();
}

std::string NeutronInteractionCsvFilename(const std::string& runLabel)
{
  if (runLabel.empty()) {
    return "neutron_interactions.csv";
  }

  std::string safeLabel = runLabel;
  std::transform(safeLabel.begin(), safeLabel.end(), safeLabel.begin(),
                 [](unsigned char c) {
                   return std::isalnum(c) || c == '_' || c == '-' ? c : '_';
                 });
  return "neutron_interactions_" + safeLabel + ".csv";
}

bool InitializeNeutronInteractionCsv(const std::string& filename)
{
  std::lock_guard<std::mutex> lock(interactionCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::trunc);
  if (!output) {
    return false;
  }
  output << NeutronInteractionCsvHeader();
  return output.good();
}

bool AppendNeutronInteractionCsv(
    const std::vector<NeutronInteractionRecord>& records,
    const std::string& filename)
{
  if (records.empty()) {
    return true;
  }

  std::ostringstream payload;
  for (const auto& record : records) {
    payload << NeutronInteractionCsvRow(record);
  }

  std::lock_guard<std::mutex> lock(interactionCsvMutex);
  std::ofstream output(filename, std::ios::out | std::ios::app);
  if (!output) {
    return false;
  }
  output << payload.str();
  return output.good();
}
