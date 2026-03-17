cmake_minimum_required(VERSION 3.14)
project(SND_sim)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple)

dd4hep_add_plugin(SiPadDetector_plugin SOURCES detector_plugin/SiPadDetector.cpp)
target_link_libraries(SiPadDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(SiTargetDetector_plugin SOURCES detector_plugin/SiTargetDetector.cpp)
target_link_libraries(SiTargetDetector_plugin DD4hep::DDCore)

gaudi_add_module(SND_reco
  SOURCES
    gaudi_source/GeV2MIPConversion.cpp
    gaudi_source/BasicDigitizer.cpp
    gaudi_source/ContributionInspector.cpp
    gaudi_source/DelayTagger.cpp
    gaudi_source/EventMerger.cpp
    gaudi_source/EventWindowSplitter.cpp
    gaudi_source/EventShuffler.cpp
    gaudi_source/EDM4HEP2RNTuple.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podioRootIO
    DD4hep::DDCore
    ROOT::ROOTNTuple
)

install(TARGETS SiPadDetector_plugin SiTargetDetector_plugin LIBRARY DESTINATION lib)
install(TARGETS SND_reco LIBRARY DESTINATION lib)
