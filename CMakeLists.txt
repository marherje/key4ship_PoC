cmake_minimum_required(VERSION 3.14)
project(SND_sim)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple)

dd4hep_add_plugin(SiPadDetector_plugin SOURCES plugins/SiPadDetector.cpp)
target_link_libraries(SiPadDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(SiTargetDetector_plugin SOURCES plugins/SiTargetDetector.cpp)
target_link_libraries(SiTargetDetector_plugin DD4hep::DDCore)

gaudi_add_module(SND_reco
  SOURCES
    reco/GeV2MIPConversion.cpp
    reco/BasicDigitizer.cpp
    reco/ContributionInspector.cpp
    reco/DelayTagger.cpp
    reco/EventMerger.cpp
    reco/EventWindowSplitter.cpp
    reco/EventShuffler.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podioRootIO
    DD4hep::DDCore
)

install(TARGETS SiPadDetector_plugin SiTargetDetector_plugin LIBRARY DESTINATION lib)
install(TARGETS SND_reco LIBRARY DESTINATION lib)
