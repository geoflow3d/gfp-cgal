cmake_minimum_required(VERSION 3.14)
project (processing-cgal VERSION 0.0.1)

find_package(geoflow REQUIRED)

# CGAL
find_package(LASlib)
find_package(glm)
find_package(CGAL QUIET COMPONENTS Core)

set(GF_PLUGIN_NAME ${PROJECT_NAME})
set(GF_PLUGIN_TARGET_NAME "gfp_cgal")
set(GF_PLUGIN_REGISTER ${PROJECT_SOURCE_DIR}/register.hpp)
geoflow_create_plugin(
  src/cgal_nodes.cpp
  src/isolines.cpp
  src/linesimp.cpp
  src/tinsimp.cpp
)
target_include_directories( gfp_cgal PRIVATE src ${glm_DIR}/.. ${LASLIB_INCLUDE_DIR})
target_link_libraries( gfp_cgal
  geoflow-core
  LASlib
  CGAL::CGAL CGAL::CGAL_Core
)
if (MSVC)
  add_definitions(-DNOMINMAX)
  target_compile_options(gfp_cgal PRIVATE "/Zc:__cplusplus")
endif()