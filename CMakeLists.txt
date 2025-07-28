cmake_minimum_required(VERSION 3.25)
project(PennerFeature)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-fPIC)
add_compile_options(-Wno-register)

if(PROJECT_IS_TOP_LEVEL)
  # Options for the project
  option(USE_PYBIND "Generate pybindings" ON)
  option(USE_MULTIPRECISION "Use high precision floating point" OFF)
  option(USE_MULTIPRECISION_OVERLAY "Use high precision floating point for overlay" ON)

  option(USE_SUITESPARSE "Use suitesparse for QR factorization" ON)
  # TODO Need to disable most features for mpfr, including suitesparse
  option(USE_COMISO "Use Comiso for field generation" ON)

  option(RENDER_TEXTURE "Render results" ON)
  option(ENABLE_VISUALIZATION "Generate viewers for visualization" ON)

  option(CHECK_VALIDITY "Check validity pre and post conditions" OFF)
  option(BUILD_CURVATURE_METRIC_TESTS "Build tests" ON)
  option(USE_EMBREE "Use embree" ON)
  option(USE_WARNINGS "Compile with warnings" ON)

  # Set libigl options
  option(LIBIGL_PREDICATES "Use Predicates" ON)
  set ( SUITESPARSE_ENABLE_PROJECTS "suitesparse_config;cholmod;spqr;umfpack" )
  option ( SUITESPARSE_USE_CUDA OFF )
endif()


list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# add core Penner methods
add_subdirectory(ext/penner-optimization)

if(NOT ENABLE_VISUALIZATION)
    option(POLYSCOPE_BACKEND_OPENGL3_GLFW OFF)
endif()

# Install executables to bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add executables
add_subdirectory(src/app)

