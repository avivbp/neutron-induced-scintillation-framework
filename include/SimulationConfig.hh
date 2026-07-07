#ifndef SIMULATION_CONFIG_HH
#define SIMULATION_CONFIG_HH

#include <string>

// Lightweight configuration record for future C++-level parameterization.
// The first package version uses Python to translate YAML input into Geant4
// macro commands.  This struct is provided so the remaining hard-coded detector
// and optical quantities can later be moved into C++ without changing the user
// facing config format.
struct SimulationConfig {
  double innerRadiusCm = 5.0;
  double innerHeightCm = 10.0;
  double outerRadiusCm = 17.5;
  double outerHeightCm = 70.0;

  double larScintYieldPerMeV = 51300.0;
  double larYieldScale = 0.01;
  double larAbsLengthCm = 150.0;
  double larFastTimeNs = 7.0;
  double larSlowTimeNs = 1500.0;
  double larFastFraction = 1.0;
  double birksConstantMmPerMeV = 0.03;

  double tpbEfficiency = 0.266;
  double tpbWLSTimeNs = 1.7;
  double tpbMeanNumberPhotons = 1.0;

  int topPMTs = 4;
  int bottomPMTs = 4;
  int sipmRows = 0;
  double pmtPDE = 0.266;
  double sipmPDE = 0.117;
  double pmtReflectivity = 0.10;
  double sipmReflectivity = 0.10;

  double tofMinNs = 40.0;
  double tofMaxNs = 50.0;

  std::string outputDir = "output/default";
  std::string configId = "default";
};

#endif
