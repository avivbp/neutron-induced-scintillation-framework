#ifndef NeutronInteraction_h
#define NeutronInteraction_h 1

#include <cstddef>
#include <string>
#include <vector>

enum class NeutronInteractionChannel {
  None,
  Elastic,
  Inelastic,
  Capture,
  Fission,
  OtherHadronic,
  TransportExit,
  Unclassified,
};

struct NeutronSecondaryTruth {
  std::string particle;
  int pdgCode = 0;
  double kineticEnergyKeV = 0.0;
};

struct NeutronInteractionRecord {
  int eventId = 0;
  std::size_t interactionIndex = 0;
  int trackId = 0;
  int parentId = 0;
  int stepNumber = 0;
  std::string volumeName;
  std::string volumeRoles;
  bool isFiducialLAr = false;
  std::string processName;
  NeutronInteractionChannel channel = NeutronInteractionChannel::None;
  double xCm = 0.0;
  double yCm = 0.0;
  double zCm = 0.0;
  double timeNs = 0.0;
  double preKineticEnergyKeV = 0.0;
  double postKineticEnergyKeV = 0.0;
  double localEnergyDepositKeV = 0.0;
  std::vector<NeutronSecondaryTruth> secondaries;
};

NeutronInteractionChannel ClassifyNeutronInteraction(
    const std::string& processName, bool isHadronic, bool exitsWorld);
const char* NeutronInteractionChannelName(NeutronInteractionChannel channel);

std::string NeutronInteractionCsvHeader();
std::string NeutronInteractionCsvRow(const NeutronInteractionRecord& record);
std::string NeutronInteractionCsvFilename(const std::string& runLabel);

bool InitializeNeutronInteractionCsv(
    const std::string& filename = "neutron_interactions.csv");
bool AppendNeutronInteractionCsv(
    const std::vector<NeutronInteractionRecord>& records,
    const std::string& filename = "neutron_interactions.csv");

#endif
