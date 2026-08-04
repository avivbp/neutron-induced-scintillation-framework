#include "NeutronInteraction.hh"

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
  Require(ClassifyNeutronInteraction("hadElastic", true, false) ==
              NeutronInteractionChannel::Elastic,
          "hadElastic was not classified as elastic");
  Require(ClassifyNeutronInteraction("neutronInelastic", true, false) ==
              NeutronInteractionChannel::Inelastic,
          "neutronInelastic was not classified as inelastic");
  Require(ClassifyNeutronInteraction("nCapture", true, false) ==
              NeutronInteractionChannel::Capture,
          "nCapture was not classified as capture");
  Require(ClassifyNeutronInteraction("nFission", true, false) ==
              NeutronInteractionChannel::Fission,
          "fission process was not classified as fission");
  Require(ClassifyNeutronInteraction("neutronGeneral", true, false) ==
              NeutronInteractionChannel::OtherHadronic,
          "unknown hadronic process was not retained");
  Require(ClassifyNeutronInteraction("Transportation", false, true) ==
              NeutronInteractionChannel::TransportExit,
          "world exit was not classified as transport_exit");
  Require(ClassifyNeutronInteraction("StepLimiter", false, false) ==
              NeutronInteractionChannel::Unclassified,
          "unknown physical process was not retained");
  Require(ClassifyNeutronInteraction("Transportation", false, false) ==
              NeutronInteractionChannel::None,
          "ordinary geometry crossing was recorded as an interaction");

  NeutronInteractionRecord record;
  record.eventId = 4;
  record.interactionIndex = 2;
  record.trackId = 8;
  record.parentId = 1;
  record.stepNumber = 7;
  record.volumeName = "active,lar";
  record.volumeRoles = "active_lar|fiducial_lar";
  record.processName = "neutronInelastic";
  record.channel = NeutronInteractionChannel::Inelastic;
  record.xCm = 1.25;
  record.timeNs = 3.5;
  record.preKineticEnergyKeV = 1000.0;
  record.postKineticEnergyKeV = 750.0;
  record.localEnergyDepositKeV = 10.0;
  record.secondaries.push_back({"gamma", 22, 240.0});

  const auto row = NeutronInteractionCsvRow(record);
  Require(row.find("\"active,lar\"") != std::string::npos,
          "CSV text field was not escaped");
  Require(row.find(",inelastic,") != std::string::npos,
          "CSV channel name is missing");
  Require(row.find(",1,gamma:22:240\n") != std::string::npos,
          "secondary truth was not serialized");
  Require(NeutronInteractionCsvFilename("") == "neutron_interactions.csv",
          "default truth filename changed");
  Require(NeutronInteractionCsvFilename("C 01/test") ==
              "neutron_interactions_C_01_test.csv",
          "run label was not made filename-safe");

  std::cout << "Neutron interaction tests passed\n";
  return EXIT_SUCCESS;
}
