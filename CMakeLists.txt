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
find_package(SQLite3 REQUIRED)

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

# Add IOC_Config subdirectory
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/IOC_Config/CMakeLists.txt")
    # Disable examples and tests when used as submodule to avoid name conflicts
    set(BUILD_EXAMPLES_OLD ${BUILD_EXAMPLES})
    set(BUILD_TESTS_OLD ${BUILD_TESTS})
    set(BUILD_EXAMPLES OFF)
    set(BUILD_TESTS OFF)
    add_subdirectory(external/IOC_Config)
    set(BUILD_EXAMPLES ${BUILD_EXAMPLES_OLD})
    set(BUILD_TESTS ${BUILD_TESTS_OLD})
    set(IOC_CONFIG_FOUND TRUE)
    message(STATUS "IOC_Config found and will be used")
else()
    set(IOC_CONFIG_FOUND FALSE)
    message(WARNING "IOC_Config not found, using built-in config_manager")
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CURL_INCLUDE_DIRS}
    ${LIBXML2_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/cspice/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
    ${CMAKE_CURRENT_SOURCE_DIR}/external/IOC_GaiaLib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/ITALOccultLibrary/italoccultlibrary/include
)

# Add IOC_Config include if available
if(IOC_CONFIG_FOUND)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/IOC_Config/include)
endif()

# Source files
set(SOURCES
    src/astdys_client.cpp
    src/gaia_client.cpp
    src/gaia_cache.cpp
    src/orbital_elements.cpp
    src/ephemeris.cpp
    src/occultation_predictor.cpp
    # src/chebyshev_approximation.cpp  # RIMOSSO: Usa quello di ITALOccultLibrary
    src/chebyshev_detector.cpp
    #src/chebyshev_occultation_manager.cpp
    src/phase1_candidate_screening.cpp
    src/phase2_occultation_geometry.cpp
    external/ITALOccultLibrary/italoccultlibrary/src/eq1_parser.cpp
    external/ITALOccultLibrary/italoccultlibrary/src/astdyn_wrapper.cpp
    external/ITALOccultLibrary/italoccultlibrary/src/orbital_conversions.cpp
    external/ITALOccultLibrary/italoccultlibrary/src/chebyshev_approximation.cpp
    external/ITALOccultLibrary/italoccultlibrary/src/chebyshev_rkf78_propagation.cpp
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
    # src/orbfit_force_model.cpp
    src/rkf78_integrator.cpp
    src/jpl_ephemeris.cpp
    src/orbit_propagator.cpp
    src/ra15_integrator.cpp
    src/jpl_horizons_client.cpp
    src/jpl_de_downloader.cpp
    src/spk_reader.cpp
    src/spice_spk_reader.cpp
    src/jpleph.cpp
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
    src/astdyn_propagation_helper.cpp
    src/occultation_search_astdyn.cpp
    src/astdyn_interface.cpp
    src/allnum_database.cpp
    # src/gaia_catalog_preloader.cpp
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
    # include/ioccultcalc/ra15_integrator.hpp  # File mancante
    include/ioccultcalc/initial_orbit.h
    include/ioccultcalc/config_manager.h
    include/ioccultcalc/asteroid_filter.h
    include/ioccultcalc/data_manager.h
    include/ioccultcalc/asteroid_database.h
    include/ioccultcalc/allnum_database.h
    include/refraction.h
    include/topocentric.h
    include/nutation.h
    include/planetary_aberration.h
    include/cubic_spline.h
    include/parallel_algorithms.h
    include/phase1_candidate_screening.h
    include/phase2_occultation_geometry.h
    # include/ioccultcalc/gaia_catalog_preloader.h  # File mancante
)

# RIMUOVI questa riga - è codice C++, non CMake
# using AstDySClient = AstDysClient;

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

# Verifica esistenza librerie opzionali
set(CSPICE_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cspice/lib/cspice.a")
set(ASTDYN_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a")

if(EXISTS ${CSPICE_LIB_PATH})
    message(STATUS "Found CSPICE: ${CSPICE_LIB_PATH}")
    set(CSPICE_LIB ${CSPICE_LIB_PATH})
else()
    message(WARNING "CSPICE not found at ${CSPICE_LIB_PATH} - some features may be unavailable")
    set(CSPICE_LIB "")
endif()

if(EXISTS ${ASTDYN_LIB_PATH})
    message(STATUS "Found AstDyn: ${ASTDYN_LIB_PATH}")
    set(ASTDYN_LIB ${ASTDYN_LIB_PATH})
else()
    message(WARNING "AstDyn not found at ${ASTDYN_LIB_PATH} - test_astdyn_integration will fail")
    message(STATUS "To build AstDyn: cd external/ITALOccultLibrary/astdyn && ./build.sh")
    set(ASTDYN_LIB "")
endif()

target_link_libraries(ioccultcalc
    CURL::libcurl
    ${LIBXML2_LIBRARIES}
    SQLite::SQLite3
    ${CSPICE_LIB}
    ${ASTDYN_LIB}
    ${ORBFIT_ALL_LIBS}
    ioc_gaialib
    OpenMP::OpenMP_CXX
)

# Link IOC_Config if available
if(IOC_CONFIG_FOUND)
    target_link_libraries(ioccultcalc ioc_config_static)
    target_compile_definitions(ioccultcalc PRIVATE IOC_CONFIG_AVAILABLE)
endif()

# Rimuovi il blocco if separato per CSPICE se esiste

target_include_directories(ioccultcalc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Trova Eigen3 (richiesto da ChebyshevOccultationManager e AstDyn)
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    target_include_directories(ioccultcalc PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ioccultcalc Eigen3::Eigen)
    message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
else()
    message(WARNING "Eigen3 not found - ChebyshevOccultationManager requires Eigen3")
endif()

# Se AstDyn è abilitato, aggiungi anche dipendenze Boost
if(USE_ASTDYN AND ASTDYN_LIBRARIES)
    find_package(Boost QUIET)
    
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

# # Test 17030 position
# add_executable(test_17030_position test_17030_position.cpp)
# target_link_libraries(test_17030_position PRIVATE ioccultcalc cspice)

# # Test 17030 effemeridi geocentriche (usa AstDyn wrapper)
# add_executable(test_17030_geocentric 
#     test_17030_geocentric.cpp
#     external/ITALOccultLibrary/italoccultlibrary/src/chebyshev_rkf78_propagation.cpp
#     external/ITALOccultLibrary/italoccultlibrary/src/astdyn_wrapper.cpp
#     external/ITALOccultLibrary/italoccultlibrary/src/eq1_parser.cpp
#     external/ITALOccultLibrary/italoccultlibrary/src/orbital_conversions.cpp
# )
# target_include_directories(test_17030_geocentric PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/italoccultlibrary/include
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/include
# )
# target_link_libraries(test_17030_geocentric PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test 17030 geocentrico usando direttamente AstDyn Propagator
# add_executable(test_17030_geocentric_astdyn 
#     test_17030_geocentric_astdyn.cpp
# )
# target_include_directories(test_17030_geocentric_astdyn PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/include
# )
# target_link_libraries(test_17030_geocentric_astdyn PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test 17030 geocentrico caricando elementi da archivio
# add_executable(test_17030_geocentric_from_archive 
#     test_17030_geocentric_from_archive.cpp
# )
# target_include_directories(test_17030_geocentric_from_archive PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/include
# )
# target_link_libraries(test_17030_geocentric_from_archive PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test 17030 effemeridi con Chebyshev polynomials
# add_executable(test_17030_chebyshev_ephemerides 
#     test_17030_chebyshev_ephemerides.cpp
# )
# target_include_directories(test_17030_chebyshev_ephemerides PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/include
# )
# target_link_libraries(test_17030_chebyshev_ephemerides PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test validazione Corridor API per 17030
# add_executable(test_17030_corridor_validation
#     test_17030_corridor_validation.cpp
# )
# target_link_libraries(test_17030_corridor_validation PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test occultazione specifica 17030 × UCAC4 552-011427
# add_executable(test_17030_occultation_ucac4
#     test_17030_occultation_ucac4_552_011427.cpp
# )
# target_link_libraries(test_17030_occultation_ucac4 PRIVATE
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     Eigen3::Eigen
# )

# # Test FASE 1: Screening stelle candidate con path ad alta risoluzione
# add_executable(test_17030_phase1_screening
#     test_17030_phase1_candidate_screening.cpp
# )
# target_link_libraries(test_17030_phase1_screening PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test classe Phase1CandidateScreening
# add_executable(test_phase1_class
#     test_phase1_class.cpp
# )
# target_link_libraries(test_phase1_class PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test Chebyshev RA/Dec + Grid Search Optimization
# add_executable(test_chebyshev_ar_dec
#     test_chebyshev_ar_dec_optimization.cpp
# )
# target_link_libraries(test_chebyshev_ar_dec PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test comparativo con 16, 32, 48 punti di controllo
# add_executable(test_chebyshev_comparison
#     test_chebyshev_comparison.cpp
# )
# target_link_libraries(test_chebyshev_comparison PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
# )

# # Test Phase 2 con orbital fitting
# add_executable(test_phase2_orbital_fitting
#     test_phase2_orbital_fitting.cpp
# )
# target_link_libraries(test_phase2_orbital_fitting PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
#     CURL::libcurl
# )

# # Test integrazione Phase 1 → Phase 2
# add_executable(test_phase1_phase2_integration
#     test_phase1_phase2_integration.cpp
# )
# target_link_libraries(test_phase1_phase2_integration PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
#     CURL::libcurl
# )

# # Test integrazione Phase 1 → Phase 2 (collegato)
# add_executable(test_phase1_phase2_linked
#     test_phase1_phase2_linked.cpp
# )
# target_link_libraries(test_phase1_phase2_linked PRIVATE
#     ioccultcalc
#     ioc_gaialib
#     ${CMAKE_SOURCE_DIR}/external/ITALOccultLibrary/astdyn/build/src/libastdyn.a
#     /opt/homebrew/lib/libboost_filesystem.dylib
#     /opt/homebrew/lib/libboost_atomic.dylib  
#     /opt/homebrew/lib/libboost_program_options.dylib
#     /opt/homebrew/lib/libboost_date_time.dylib
#     Eigen3::Eigen
#     CURL::libcurl
# )

