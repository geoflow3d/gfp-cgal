cmake_minimum_required(VERSION 3.14)
project (processing-cgal VERSION 0.0.1)

add_definitions(-DGF_PLUGIN_NAME=\"${CMAKE_PROJECT_NAME}\")

find_package(geoflow REQUIRED)

# CGAL
find_package(LASlib)
find_package(glm)
find_package(CGAL QUIET COMPONENTS Core)

add_library(gfp_cgal MODULE
  cgal_nodes.cpp
  isolines.cpp
  linesimp.cpp
  tinsimp.cpp
  plugin.cpp
)
target_link_libraries( gfp_cgal
  geoflow-core
  LASlib
  CGAL::CGAL CGAL::CGAL_Core
)
set_target_properties(gfp_cgal PROPERTIES CXX_STANDARD 17)
install(TARGETS gfp_cgal DESTINATION ${geoflow_PLUGIN_FOLDER})