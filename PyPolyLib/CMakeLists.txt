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
    # need this on some MacOS or else the system directories will still come
    # first, and if you have an old polylib installed it will fail...
    set_property(TARGET PkgConfig::POLYLIB PROPERTY SYSTEM FALSE)
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
# message(STATUS "DEBUG: PolyLib include dirs: ${POLYLIB_INCLUDE_DIRS}")
# message(STATUS "DEBUG: PolyLib libraries: ${POLYLIB_LIBRARIES}")
# message(STATUS "DEBUG: PolyLib libraries dirs: ${POLYLIB_LIBRARY_DIRS}")
# message(STATUS "DEBUG: PolyLib LDFLAGS: ${POLYLIB_LDFLAGS}")


# ----------------------------------------------------------------------
# Python extension build
# ----------------------------------------------------------------------
pybind11_add_module(_core
    src/_core.cpp
)

# make sure to use the uninstalled header files of polylib, if set:
target_include_directories(_core BEFORE PRIVATE
    ${POLYLIB_INCLUDE_DIRS})

target_link_libraries(_core PRIVATE
    PkgConfig::POLYLIB
    PkgConfig::GMP
)

# ----------------------------------------------------------------------
# Python package
# ----------------------------------------------------------------------
file(GLOB PYTHON_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/python/*.py"
)

install(TARGETS _core
    LIBRARY DESTINATION pypolylib
)

install(FILES
    ${PYTHON_SOURCES}
    DESTINATION pypolylib
)
