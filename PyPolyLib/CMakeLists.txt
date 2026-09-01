cmake_minimum_required(VERSION 3.18)

project(pypolylib LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)


# ----------------------------------------------------------------------
# find PolyLib using pkg-config (looks for polylibgmp[-uninstalled].pc)
#
# If POLYLIB_USE_UNINSTALLED is specified, look for the uninstalled lib
# in POLYLIB_ROOT and in POLYLIB_ROOT/build
# The user can manually set PKG_CONFIG_PATH to search other directories
#
# Otherwise, look for an installed PolyLib.
# ----------------------------------------------------------------------

option(POLYLIB_USE_UNINSTALLED
    "Use the uninstalled PolyLib build"
    OFF
)
set(POLYLIB_ROOT
    "${CMAKE_CURRENT_SOURCE_DIR}/../build"
    CACHE PATH
    "PolyLib source/build or installation directory"
)

if(POLYLIB_USE_UNINSTALLED)
    set(ENV{PKG_CONFIG_PATH}
        "${POLYLIB_ROOT}/build:${POLYLIB_ROOT}:$ENV{PKG_CONFIG_PATH}"
    )
    pkg_check_modules(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp-uninstalled
    )
else()
    pkg_check_modules(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp
    )
endif()


# message(STATUS "PolyLib include directory: ${POLYLIB_INCLUDE_DIR}")


# ----------------------------------------------------------------------
# Python extension
# ----------------------------------------------------------------------

pybind11_add_module(_core
    src/_core.cpp
)
target_link_libraries(_core PRIVATE
    PkgConfig::POLYLIB
)
# build in the right place
set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/pypolylib"
)


set(PYTHON_SOURCES
    __init__.py
    lbl.py
    lbl_io.py
    lbl_plot.py
)


add_custom_target(pypolylib_python ALL
COMMAND ${CMAKE_COMMAND} -E copy_directory
"${CMAKE_CURRENT_SOURCE_DIR}/src/python"
"${CMAKE_BINARY_DIR}/pypolylib"
)
add_dependencies(pypolylib_python _core)
