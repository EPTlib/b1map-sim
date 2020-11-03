/*****************************************************************************
*
*     Program: b1map-sim
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020  Alessandro Arduino
*  Istituto Nazionale di Ricerca Metrologica (INRiM)
*  Strada delle cacce 91, 10135 Torino
*  ITALY
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*
*****************************************************************************/

#ifndef B1MAPSIM_UTIL_H_
#define B1MAPSIM_UTIL_H_

#include <complex>
#include <functional>
#include <numeric>
#include <string>

/// Number of spatial dimensions.
constexpr int NDIM = 3;
/// Pi.
constexpr double PI = 3.14159265358979323846;
/// Speed of light [m/s].
constexpr double C0 = 299792458.0;
/// Vacuum permeability [H/m].
constexpr double MU0 = 4.0e-7*PI;
/// Vacuum permittivity [F/m].
constexpr double EPS0 = 1.0/MU0/C0/C0;
/// Gyromagnetic ratio of proton [rad/ms/T].
constexpr double GAMMA = 267522.18744;

/**
 * Return the license boilerplate as a string.
 * 
 * @return the license boilerplate.
 */
const std::string LicenseBoilerplate();

/**
 * Compute the sum of all the elements in a container.
 * 
 * @tparam T container typename.
 * 
 * @param v container of elements.
 * 
 * @return the sum of all the elements in `v'.
 */
template <typename T>
inline typename T::value_type Sum(const T &v) {
    using type = typename T::value_type;
    type result = 0;
    for (auto it = v.begin(); it<v.end(); ++it) {
        if (*it == *it) {
            result += *it;
        }
    }
    return result;
}
/**
 * 
 */
template <typename T>
inline typename T::value_type Avg(const T &v) {
    using type = typename T::value_type;
    type result = 0;
    int num = 0;
    for (auto it = v.begin(); it<v.end(); ++it) {
        if (*it == *it) {
            result += *it;
            ++num;
        }
    }
    return result/num;
}
/**
 * Compute the products of all the elements in a container.
 * 
 * @tparam T container typename.
 * 
 * @param v container of elements.
 * 
 * @return the products of all the elements in `v'.
 */
template <typename T>
inline typename T::value_type Prod(const T &v) {
    using type = typename T::value_type;
    type result = 1;
    for (auto it = v.begin(); it<v.end(); ++it) {
        if (*it == *it) {
            result *= *it;
        }
    }
    return result;
}
/**
 * 
 */
template <typename T>
inline typename T::value_type Norm2(const T &v) {
    using type = typename T::value_type;
    type result = 0.0;
    for (auto it = v.begin(); it<v.end(); ++it) {
        if (*it == *it) {
            result += std::abs(*it)*std::abs(*it);
        }
    }
    return result;
}
/**
 * 
 */
template <typename T>
inline typename T::value_type DiffNorm2(const T &v, const T &u) {
    using type = typename T::value_type;
    type result = 0.0;
    for (auto itv = v.begin(), itu = u.begin(); itv<v.end(); ++itv, ++itu) {
        if (*itv == *itv && *itu == *itu) {
            type tmp = *itv-*itu;
            result += std::abs(tmp)*std::abs(tmp);
        }
    }
    return result;
}

/**
 * Translates from multi-index to index assuming the first index the fastest.
 * 
 * @tparam T,U iterator typenames.
 * 
 * @param ii multi-index to the voxel.
 * @param nn number of voxels in each direction.
 * 
 * @return the single index to the voxel.
 * 
 * Arguments `ii' and `nn' must have the methods `begin', `end', `size' and
 * `operator[]' (any sequence containers from STL works fine).
 * 
 * @deprecated
 */
template <typename T, typename U>
int MultiIdxToIdx(const T &ii, const U &nn);
/**
 * Translates from index to multi-index assuming the first index the fastest.
 * 
 * @tparam T,U iterator typenames.
 * 
 * @param[out] ii multi-index to the voxel.
 * @param[in] idx single index to the voxel
 * @param[in] nn number of voxels in each direction.
 * 
 * Arguments `ii' and `nn' must have the methods `begin', `end', `size' and
 * `operator[]' (any sequence containers from STL works fine).
 * 
 * @decrecated
 */
template <typename T, typename U>
void IdxToMultiIdx(T &ii, int idx, const U &nn);

// ---------------------------------------------------------------------------
// -------------------------  Implementation detail  -------------------------
// ---------------------------------------------------------------------------

// Return the license boilerplate as a string.
inline const std::string LicenseBoilerplate() {
    const std::string boilerplate = "MIT License\n"
        "Copyright (c) 2020  Alessandro Arduino\n"
        "Istituto Nazionale di Ricerca Metrologica (INRiM)\n"
        "\n"
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
        "SOFTWARE.\n";
    return boilerplate;
}

// Multi-index to index
template <typename T, typename U>
int MultiIdxToIdx(const T &ii, const U &nn) {
    assert(ii.size()==nn.size());
    auto n_dim = ii.size();
    int idx = 0;
    for (decltype(n_dim) d = n_dim; d>0; --d) {
        idx = ii[d-1] + idx*nn[d-1];
    }
    return idx;
}
// Index to multi-index
template <typename T, typename U>
void IdxToMultiIdx(T &ii, int idx, const U &nn) {
    assert(ii.size()==nn.size());
    auto n_dim = ii.size();
    for (decltype(n_dim) d = 0; d<n_dim; ++d) {
        ii[d] = idx%nn[d];
        idx = idx/nn[d];
    }
    return;
}

#endif  // B1MAPSIM_UTIL_H_
