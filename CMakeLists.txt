cmake_minimum_required(VERSION 3.14)
project (processing-cgal VERSION 0.0.1)

add_definitions(-DGF_PLUGIN_NAME=\"${CMAKE_PROJECT_NAME}\")

find_package(geoflow REQUIRED)

# CGAL
find_package(LASlib)
find_package(glm)
find_package(CGAL QUIET COMPONENTS Core)

add_library(gfp_cgal MODULE
  src/cgal_nodes.cpp
  src/isolines.cpp
  src/linesimp.cpp
  src/tinsimp.cpp
  plugin.cpp
)
target_include_directories( gfp_cgal PRIVATE src ${glm_DIR}/.. ${LASLIB_INCLUDE_DIR})
target_link_libraries( gfp_cgal
  geoflow-core
  LASlib
  CGAL::CGAL CGAL::CGAL_Core
)
set_target_properties(gfp_cgal PROPERTIES CXX_STANDARD 17)
if (MSVC)
  target_compile_options(gfp_cgal PRIVATE "/Zc:__cplusplus")
endif()

install(TARGETS gfp_cgal DESTINATION ${geoflow_PLUGIN_FOLDER})