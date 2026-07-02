cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(fos-parameterization
        VERSION 1.0.0
        LANGUAGES Fortran CXX
        DESCRIPTION "Fourier-over-Spheroid (FoS) nuclear shape parameterization"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        GCC-Compiler-Options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/GCC-Compiler-Options.git
        GIT_TAG 1.4.1
)
FetchContent_MakeAvailable(GCC-Compiler-Options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        Fortran-Foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/Fortran-Foundations.git
        GIT_TAG 2.3.0
)
FetchContent_MakeAvailable(Fortran-Foundations)

# --- Preprocessor & build type ---
set(CMAKE_Fortran_PREPROCESS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags via GCC-Compiler-Options ---
include(GCCCompilerOptions/FortranCompilerOptions)
create_fortran_library_interface(TARGET fos_parameterization_flags)

# --- Source files ---
set(FOS_PARAM_SOURCES
        src/fos_parameterization_mod.f08
)

# --- Static library target ---
add_library(fos_parameterization STATIC ${FOS_PARAM_SOURCES})
target_link_libraries(fos_parameterization
        PUBLIC fortran_foundations
        PRIVATE fos_parameterization_flags
)
target_compile_definitions(fos_parameterization PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(fos_parameterization PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
target_include_directories(fos_parameterization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(FosParameterization::fos_parameterization ALIAS fos_parameterization)

# --- Shared library target (for FFI consumers: Python ctypes, etc.) ---
add_library(fos_parameterization_shared SHARED ${FOS_PARAM_SOURCES})
target_link_libraries(fos_parameterization_shared
        PUBLIC fortran_foundations
        PRIVATE fos_parameterization_flags
)
target_compile_definitions(fos_parameterization_shared PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(fos_parameterization_shared PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared
        OUTPUT_NAME fos_parameterization
)
target_include_directories(fos_parameterization_shared
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(FosParameterization::fos_parameterization_shared
        ALIAS fos_parameterization_shared)
