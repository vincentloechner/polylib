cmake_minimum_required(VERSION 3.18)

project(pypolylib LANGUAGES CXX)

find_package(PkgConfig REQUIRED)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

# ----------------------------------------------------------------------
# find PolyLib using pkg-config (look for polylibgmp[-uninstalled].pc)
#
# If POLYLIB_USE_UNINSTALLED is ON, look for polylibgmp-uninstalled.pc in
# POLYLIB_DIR (defaut: ../build).
# Else, look for polylibgmp.pc in POLYLIB_DIR and standard paths.
# The user can also manually set PKG_CONFIG_PATH to search other directories.
# ----------------------------------------------------------------------

option(POLYLIB_USE_UNINSTALLED
    "Use the uninstalled PolyLib build (pkg-config file polylibgmp-uninstalled.pc)"
    OFF
)
set(POLYLIB_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../build"
    CACHE PATH
    "PolyLib source/build or installation directory"
)

if(POLYLIB_USE_UNINSTALLED)
    set(ENV{PKG_CONFIG_PATH}
        "${POLYLIB_DIR}:$ENV{PKG_CONFIG_PATH}"
    )
    pkg_check_modules(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp-uninstalled
    )
else()
    # TODO: add POLYLIB_DIR to the path if user-specified
    pkg_check_modules(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp
    )
endif()


message(STATUS "PolyLib include dirs: ${POLYLIB_INCLUDE_DIRS}")
message(STATUS "PolyLib libraries: ${POLYLIB_LIBRARIES}")
message(STATUS "PolyLib libraries dirs: ${POLYLIB_LIBRARY_DIRS}")


# ----------------------------------------------------------------------
# Python extension
# ----------------------------------------------------------------------
pybind11_add_module(_core
    src/_core.cpp
)


# make sure to use the given polylib include dir,
# in case an older version of polylib is installed system-wide:
target_include_directories(_core BEFORE PRIVATE
    ${POLYLIB_INCLUDE_DIRS}
)
target_link_directories(_core PRIVATE
    ${POLYLIB_LIBRARY_DIRS}
)
target_link_libraries(_core PRIVATE
    ${POLYLIB_LIBRARIES}
    PkgConfig::GMP
)
# instead of this (cmake uses the system directories first by default)
# target_link_libraries(_core PRIVATE
#     PkgConfig::POLYLIB
#     PkgConfig::GMP
# )

# move to the right place
set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/pypolylib"
)

# move the python sources:
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
