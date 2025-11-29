cmake_minimum_required(VERSION 3.15)
project(IOccultCalc VERSION 1.0.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fortran settings for OrbFit wrapper
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)

# Options
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Find required packages
find_package(CURL REQUIRED)
find_package(LibXml2 REQUIRED)

# OpenMP support (macOS requires libomp from Homebrew)
if(APPLE)
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()
find_package(OpenMP REQUIRED)

# Add IOC_GaiaLib subdirectory
add_subdirectory(external/IOC_GaiaLib)

# Add IOC_SPICE subdirectory
add_subdirectory(external/IOC_SPICE)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CURL_INCLUDE_DIRS}
    ${LIBXML2_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/cspice/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
    ${CMAKE_CURRENT_SOURCE_DIR}/external/IOC_GaiaLib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/IOC_SPICE/include
)

# Source files
set(SOURCES
    src/astdys_client.cpp
    # gaia_client.cpp and gaia_cache.cpp now provided by IOC_GaiaLib via adapter
    src/gaia_adapter.cpp
    src/orbital_elements.cpp
    src/ephemeris.cpp
    src/occultation_predictor.cpp
    src/chebyshev_approximation.cpp
    src/chebyshev_detector.cpp
    src/kml_exporter.cpp
    src/occult4_xml.cpp
    src/prediction_report.cpp
    src/preston_parser.cpp
    src/coordinates.cpp
    src/time_utils.cpp
    src/observation.cpp
    src/mpc_client.cpp
    src/orbit_determination.cpp
    src/orbit_fitter.cpp
    src/orbit_determination_gauss.cpp
    src/force_model.cpp
    src/orbfit_force_model.cpp
    src/rkf78_integrator.cpp
    src/jpl_ephemeris.cpp
    src/orbit_propagator.cpp
    src/jpl_horizons_client.cpp
    src/jpl_de_downloader.cpp
    src/spk_reader.cpp
    src/spice_spk_reader.cpp
    src/jpleph.cpp
    src/ra15_integrator.cpp
    src/initial_orbit.cpp
    src/orbfit_wrapper.cpp
    src/iers_data.cpp
    src/config_manager.cpp
    src/asteroid_filter.cpp
    src/data_manager.cpp
    src/asteroid_database.cpp
    src/refraction.cpp
    src/topocentric.cpp
    src/nutation.cpp
    src/planetary_aberration.cpp
    src/cubic_spline.cpp
    src/parallel_algorithms.cpp
    src/output_manager.cpp
    src/asteroid_list_reader.cpp
    src/gaia_catalog_preloader.cpp
)

# Fortran wrapper for OrbFit
set(FORTRAN_SOURCES
    src/orbfit_c_wrapper.f90
)

# Header files
set(HEADERS
    include/ioccultcalc/astdys_client.h
    include/ioccultcalc/gaia_client.h
    include/ioccultcalc/gaia_cache.h
    include/ioccultcalc/orbital_elements.h
    include/ioccultcalc/ephemeris.h
    include/ioccultcalc/occultation_predictor.h
    include/ioccultcalc/kml_exporter.h
    include/ioccultcalc/occult4_xml.h
    include/ioccultcalc/prediction_report.h
    include/ioccultcalc/preston_parser.h
    include/ioccultcalc/coordinates.h
    include/ioccultcalc/time_utils.h
    include/ioccultcalc/types.h
    include/ioccultcalc/observation.h
    include/ioccultcalc/mpc_client.h
    include/ioccultcalc/orbit_determination.h
    include/ioccultcalc/force_model.h
    include/ioccultcalc/jpl_ephemeris.h
    include/ioccultcalc/orbit_propagator.h
    include/ioccultcalc/jpl_horizons_client.h
    include/ioccultcalc/jpl_de_downloader.h
    include/ioccultcalc/spk_reader.h
    include/ioccultcalc/spice_spk_reader.h
    include/ioccultcalc/jpleph.h
    include/ioccultcalc/jpl_int.h
    include/ioccultcalc/get_bin.h
    include/ioccultcalc/orbfit_wrapper.h
    include/ioccultcalc/numerical_integrator.h
    include/ioccultcalc/relativistic_corrections.h
    include/ioccultcalc/star_catalog.h
    include/ioccultcalc/asteroid_shape.h
    include/ioccultcalc/uncertainty_propagation.h
    include/ioccultcalc/ra15_integrator.hpp
    include/ioccultcalc/initial_orbit.h
    include/ioccultcalc/config_manager.h
    include/ioccultcalc/asteroid_filter.h
    include/ioccultcalc/data_manager.h
    include/ioccultcalc/asteroid_database.h
    include/refraction.h
    include/topocentric.h
    include/nutation.h
    include/planetary_aberration.h
    include/cubic_spline.h
    include/parallel_algorithms.h
    include/ioccultcalc/gaia_catalog_preloader.h
)

# OrbFit integration via Fortran wrapper with bind(C)
option(USE_ORBFIT "Enable OrbFit integration" OFF)

if(USE_ORBFIT)
    set(ORBFIT_PATH "/Users/michelebigi/Astro/OrbFit" CACHE PATH "Path to OrbFit installation")
    
    # Find Fortran compiler
    if(NOT CMAKE_Fortran_COMPILER)
        message(FATAL_ERROR "Fortran compiler not found! Install gfortran.")
    endif()
    
    # OrbFit libraries
    set(ORBFIT_LIBS
        ${ORBFIT_PATH}/src/lib/libgauss.a
        ${ORBFIT_PATH}/src/lib/libprop.a
        ${ORBFIT_PATH}/src/lib/libsuit.a
        ${ORBFIT_PATH}/src/lib/libmoid.a
    )
    
    # Find Fortran runtime libraries
    execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.dylib
        OUTPUT_VARIABLE GFORTRAN_LIB
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT EXISTS "${GFORTRAN_LIB}")
        execute_process(
            COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.a
            OUTPUT_VARIABLE GFORTRAN_LIB
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    
    # Find quadmath (for quad-precision)
    execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libquadmath.dylib
        OUTPUT_VARIABLE QUADMATH_LIB
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT EXISTS "${QUADMATH_LIB}")
        execute_process(
            COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libquadmath.a
            OUTPUT_VARIABLE QUADMATH_LIB
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    
    message(STATUS "OrbFit integration ENABLED")
    message(STATUS "  OrbFit path: ${ORBFIT_PATH}")
    message(STATUS "  gfortran lib: ${GFORTRAN_LIB}")
    message(STATUS "  quadmath lib: ${QUADMATH_LIB}")
    
    # Compile Fortran wrapper
    add_library(orbfit_wrapper_fortran STATIC ${FORTRAN_SOURCES})
    target_include_directories(orbfit_wrapper_fortran PUBLIC
        ${ORBFIT_PATH}/src/include
        ${ORBFIT_PATH}/src/propag    # Per propag_state.mod
        ${ORBFIT_PATH}/src/prelim    # Per prelim modules
        ${ORBFIT_PATH}/src/suit      # Per utility modules
    )
    
    set(ORBFIT_WRAPPER_LIB orbfit_wrapper_fortran)
    set(ORBFIT_ALL_LIBS ${ORBFIT_WRAPPER_LIB} ${ORBFIT_LIBS} ${GFORTRAN_LIB} ${QUADMATH_LIB})
else()
    message(STATUS "OrbFit integration DISABLED (use -DUSE_ORBFIT=ON to enable)")
    set(ORBFIT_ALL_LIBS "")
endif()

# Create library
add_library(ioccultcalc STATIC ${SOURCES} ${HEADERS})

target_link_libraries(ioccultcalc
    CURL::libcurl
    ${LIBXML2_LIBRARIES}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/cspice/lib/cspice.a
    ${ORBFIT_ALL_LIBS}
    ioc_gaialib
    ioc_spice
    OpenMP::OpenMP_CXX
)

target_include_directories(ioccultcalc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Se AstDyn è abilitato, aggiungi dipendenze Eigen e Boost
if(USE_ASTDYN AND ASTDYN_LIBRARIES)
    find_package(Eigen3 QUIET)
    find_package(Boost QUIET)
    
    if(Eigen3_FOUND)
        target_include_directories(ioccultcalc PRIVATE ${EIGEN3_INCLUDE_DIR})
        target_link_libraries(ioccultcalc Eigen3::Eigen)
    endif()
    
    if(Boost_FOUND)
        target_include_directories(ioccultcalc PRIVATE ${Boost_INCLUDE_DIRS})
        target_link_libraries(ioccultcalc ${Boost_LIBRARIES})
    endif()
endif()

# Examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Tools
add_subdirectory(tools)

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
install(TARGETS ioccultcalc
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ioccultcalc
    DESTINATION include
)
