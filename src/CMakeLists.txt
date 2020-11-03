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

set(B1MAPSIM_SRC
    b1mapping.cc
    body.cc
    image.cc
    main.cc
    sequences.cc
    util.cc
    version.cc
    io/io_hdf5.cc
    io/io_toml.cc
    io/io_util.cc)

add_executable(b1map-sim ${B1MAPSIM_SRC})

target_include_directories(b1map-sim
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/tinytoml/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
        $<INSTALL_INTERFACE:share/b1map/extern/tinytoml/include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src)

if(${CMAKE_TOOLCHAIN_FILE} MATCHES "vcpkg")
    target_link_libraries(b1map-sim PUBLIC ${HDF5_LIBRARIES})
else()
    target_link_libraries(b1map-sim PUBLIC hdf5 hdf5_cpp)
endif()

target_compile_features(b1map-sim PUBLIC cxx_std_11)

set_property(TARGET b1map-sim
    PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
