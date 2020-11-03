#=============================================================================
#
#     Program: b1map-sim
#     Author: Alessandro Arduino <a.arduino@inrim.it>
#
#  MIT License
#
#  Copyright (c) 2020  Alessandro Arduino
#  Istituto Nazionale di Ricerca Metrologica (INRiM)
#  Strada delle cacce 91, 10135 Torino
#  ITALY
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#=============================================================================

cmake_minimum_required(VERSION 3.16)

project(b1map-sim VERSION 0.1.0 LANGUAGES C CXX)

# Disallow in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR
        "FATAL: Out-of-source build required. Create a separate directory for build files.")
endif()

# Set some CMAKE_ variables
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Configure files
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file(${PROJECT_SOURCE_DIR}/include/b1map/version.h.in
    ${PROJECT_BINARY_DIR}/include/b1map/version.h
    @ONLY)

# Find packages
# hdf5
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 COMPONENTS C CXX REQUIRED)
include_directories("${HDF5_INCLUDE_DIRS}")

# Set the applications
add_subdirectory(src)

# ----- Installing -----
include(InstallRequiredSystemLibraries)

set(CMAKE_INSTALL_MFC_LIBRARIES ON)

install(FILES ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
    DESTINATION bin)

install(TARGETS b1map-sim RUNTIME DESTINATION bin)

# ----- Packaging -----

set(CPACK_PACKAGE_NAME "b1map-sim")
set(CPACK_PACKAGE_VENDOR "Alessandro Arduino")
set(CPACK_PACKAGE_CONTACT "Alessandro Arduino")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_DESCRIPTION "A generator of virtual B1-mapping images.")
set(CPACK_PACKAGE_HOMEPAGE_URL "")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")

set(CPACK_WIX_UPGRADE_GUID "4F21507D-EDB1-4864-82A0-E078CD590B03")
configure_file(${PROJECT_SOURCE_DIR}/LICENSE ${PROJECT_BINARY_DIR}/LICENSE.txt COPYONLY)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/LICENSE.txt")
set(CPACK_WIX_CMAKE_PACKAGE_REGISTRY "b1map-sim")

set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /.vscode
    /.vs
    /*out*
    /*build*
    /*install*
)

if(WIN32)
    set(CPACK_GENERATOR "WIX")
endif()

include(CPack)
