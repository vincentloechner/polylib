cmake_minimum_required(VERSION 3.18)

project(pypolylib LANGUAGES CXX)

find_package(PkgConfig REQUIRED)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

# ----------------------------------------------------------------------
# find PolyLib using pkg-config
# If POLYLIB_USE_UNINSTALLED is ON, look for polylibgmp-uninstalled.pc in
# POLYLIB_DIR, .., ../build and standard paths.
# Else, look for polylibgmp.pc or polylibgmp-uninstalled.pc in POLYLIB_DIR
# and standard paths.
# ----------------------------------------------------------------------

option(POLYLIB_USE_UNINSTALLED
    "Force use of the uninstalled PolyLib build tree"
    OFF
)
set(POLYLIB_DIR
    ""
    CACHE PATH
    "PolyLib build or installation directory"
)

if(POLYLIB_USE_UNINSTALLED)
    set(ENV{PKG_CONFIG_PATH}
        "${POLYLIB_DIR}:..:../build:$ENV{PKG_CONFIG_PATH}"
    )
    pkg_check_modules(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp-uninstalled
    )
else()
    if(POLYLIB_DIR)
        set(ENV{PKG_CONFIG_PATH}
            "${POLYLIB_DIR}:$ENV{PKG_CONFIG_PATH}")
    endif()
    pkg_search_module(POLYLIB REQUIRED IMPORTED_TARGET
        polylibgmp
        polylibgmp-uninstalled
    )
endif()
# tell the user which .pc file was found:
pkg_get_variable(POLYLIB_PC_DIR ${POLYLIB_MODULE_NAME} pcfiledir)
message(STATUS "Using POLYLIB via pkg-config module from: \
${POLYLIB_PC_DIR}/${POLYLIB_MODULE_NAME}.pc")


# DEBUG: (to be removed)
message(STATUS "DEBUG: PolyLib include dirs: ${POLYLIB_INCLUDE_DIRS}")
message(STATUS "DEBUG: PolyLib libraries: ${POLYLIB_LIBRARIES}")
message(STATUS "DEBUG: PolyLib libraries dirs: ${POLYLIB_LIBRARY_DIRS}")


# ----------------------------------------------------------------------
# Python extension
# ----------------------------------------------------------------------
pybind11_add_module(_core
    src/_core.cpp
)


# configure polylib explicitly: make sure to use the given polylib include
# dir, in case an older version of polylib is installed system-wide
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
