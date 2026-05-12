cmake_minimum_required(VERSION 3.14)
project(SND_sim)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple)
find_package(Acts 44.3.0 REQUIRED COMPONENTS Core PluginDD4hep)
find_package(k4ActsTracking REQUIRED)

dd4hep_add_plugin(SiPadDetector_plugin SOURCES detector_plugin/SiPadDetector.cpp)
target_link_libraries(SiPadDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(SiTargetDetector_plugin SOURCES detector_plugin/SiTargetDetector.cpp)
target_link_libraries(SiTargetDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(MTCDetector_plugin SOURCES detector_plugin/MTCDetector.cpp)
target_link_libraries(MTCDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(SND_SciFiAction_plugin SOURCES simulation/ddg4/SND_SciFiAction.cpp)
target_link_libraries(SND_SciFiAction_plugin DD4hep::DDG4 DD4hep::DDCore)
set_target_properties(SND_SciFiAction_plugin PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

gaudi_add_module(SND_reco
  SOURCES
    gaudi_source/GeV2MIPConversion.cpp
    gaudi_source/BasicDigitizer.cpp
    gaudi_source/SciFiDigitizer.cpp
    gaudi_source/ContributionInspector.cpp
    gaudi_source/EventWindowSplitter.cpp
    gaudi_source/EventShuffler.cpp
    gaudi_source/EDM4HEP2RNTuple.cpp
    gaudi_source/ACTSGeoSvc.cpp
    gaudi_source/ACTSProtoTracker.cpp
    gaudi_source/SiTargetMeasConverter.cpp
    gaudi_source/SiPadMeasConverter.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podioRootIO
    DD4hep::DDCore
    DD4hep::DDRec
    ROOT::ROOTNTuple
    Acts::Core
    Acts::PluginDD4hep
    k4ActsTracking::k4ActsTracking
)

 set_target_properties(SND_reco PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN      TRUE
)
target_include_directories(SND_reco PRIVATE
    /cvmfs/sw.hsf.org/key4hep/releases/2026-02-01/x86_64-almalinux9-gcc14.2.0-opt/acts/44.3.0-37xe73/include
)

install(TARGETS SiPadDetector_plugin SiTargetDetector_plugin MTCDetector_plugin SND_SciFiAction_plugin LIBRARY DESTINATION lib)
install(TARGETS SND_reco LIBRARY DESTINATION lib)
