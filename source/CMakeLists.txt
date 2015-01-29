# - Top level CMakeLists.txt for SuperNEMO gamma tracking studies

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(snemo_gamma_tracking_efficiency)

# Use C++11
set(CMAKE_CXX_FLAGS "-W -Wall -std=c++14")

# - Third party
find_package(Falaise 1.0.0 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR} ${Falaise_INCLUDE_DIRS})

add_library(snemo_gamma_tracking_efficiency SHARED
  snemo_gamma_tracking_efficiency_module.h snemo_gamma_tracking_efficiency_module.cc)

target_link_libraries(snemo_gamma_tracking_efficiency ${Falaise_LIBRARIES})

install(FILES
  ${PROJECT_BINARY_DIR}/libsnemo_gamma_tracking_efficiency${CMAKE_SHARED_LIBRARY_SUFFIX}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE
  )

# end of CMakeLists.txt