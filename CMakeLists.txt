cmake_minimum_required(VERSION 3.14)
project(SND_sim)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)

dd4hep_add_plugin(SiPixelDetector_plugin SOURCES plugins/SiPixelDetector.cpp)
target_link_libraries(SiPixelDetector_plugin DD4hep::DDCore)

dd4hep_add_plugin(SiTargetDetector_plugin SOURCES plugins/SiTargetDetector.cpp)
target_link_libraries(SiTargetDetector_plugin DD4hep::DDCore)

gaudi_add_module(SND_reco
  SOURCES
    reco/GeV2MIPConversion.cpp
    reco/BasicDigitizer.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
)

install(TARGETS SiPixelDetector_plugin SiTargetDetector_plugin LIBRARY DESTINATION lib)
install(TARGETS SND_reco LIBRARY DESTINATION lib)
