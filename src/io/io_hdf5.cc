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

#include "b1map/io/io_hdf5.h"

#include <algorithm>

using namespace b1map;
using namespace b1map::io;

namespace {

    /**
     * Provide the uri, given url and urn.
     * 
     * @param url url
     * @param urn urn
     * 
     * @return uri
     */
    inline std::string URI(const std::string &url, const std::string &urn) {
        std::string uri("/"+url+"/"+urn);
        size_t pos = uri.find("//");
        while (pos != std::string::npos) {
            uri.replace(pos, 2, "/");
            pos = uri.find("//");
        }
        return uri;
    }

    /**
     * Create an hdf5 group, given an url.
     * 
     * @param file hdf5 file
     * @param url url
     * 
     * @return hdf5 group
     */
    inline H5::Group CreateGroup(const H5::H5File &file, const std::string &url) {
        size_t depth = 0;
        size_t snip = url.find_first_of("/", 0) ? 0 : 1;
        size_t snap = url.find_first_of("/", snip);
        std::string subpath = url.substr(snip,snap-snip);
        H5::Group group;
        while (!subpath.empty()) {
            try {
                group = depth ? group.openGroup(subpath) : file.openGroup(subpath);
            } catch (const H5::Exception&) {
                group = depth ? group.createGroup(subpath) : file.createGroup(subpath);
            }
            snip = ++snap;
            snap = url.find_first_of("/",snip);
            subpath = url.substr(snip,snap-snip);
            depth++;
        }
        return group;
    }

    // HDF5 types traits
    template <typename T>
    struct HDF5Types;
    // traits specialisations
    template <>
    struct HDF5Types<size_t> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_ULONG;
        }
    };
    template<>
    struct HDF5Types<double> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_DOUBLE;
        }
    };
    template<>
    struct HDF5Types<float> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_FLOAT;
        }
    };
    template<>
    struct HDF5Types<int> {
        static const H5::DataType Type()  {
            return H5::PredType::NATIVE_INT;
        }
    };
    template<>
    struct HDF5Types<long> {
        static const H5::DataType Type() {
            return H5::PredType::NATIVE_LONG;
        }
    };

}  //

// IOh5 constructor
IOh5::
IOh5(const std::string &fname, const Mode mode) :
    fname_(fname), mode_(mode) {
    H5::Exception::dontPrint();
    switch (mode_) {
        case Mode::In:
            file_ = H5::H5File(fname_, H5F_ACC_RDONLY);
            break;
        case Mode::Out:
            file_ = H5::H5File(fname_, H5F_ACC_TRUNC);
            break;
        case Mode::Append:
            try {
                file_ = H5::H5File(fname_, H5F_ACC_RDWR);
            } catch (const H5::FileIException&) {
                file_ = H5::H5File(fname_, H5F_ACC_TRUNC);
            }
            break;
    }
    return;
}

// IOh5 destructor
IOh5::
~IOh5() {
    file_.close();
    return;
}

// IOh5 read dataset
template <typename T>
State IOh5::
ReadDataset(Image<T> *img, const std::string &url, const std::string &urn) {
    H5::Exception::dontPrint();
    try {
        // locate the dataset
        H5::DataSet dset = file_.openDataSet(URI(url,urn));
        H5::DataSpace dspace = dset.getSpace();
        // read the data
        std::vector<hsize_t> dims(dspace.getSimpleExtentNdims());
        size_t ndim = dspace.getSimpleExtentDims(dims.data(),NULL);
        std::vector<int> nn(dims.size());
        std::reverse_copy(dims.begin(),dims.end(),nn.begin());
        *img = Image<T>(nn);
        dset.read(img->GetData().data(),::HDF5Types<T>::Type());
    } catch (const H5::FileIException&) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException&) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException&) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException&) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

// IOh5 write dataset
template <typename T>
State IOh5::
WriteDataset(const Image<T> &img, const std::string &url, const std::string &urn) const {
    H5::Exception::dontPrint();
    try {
        // open or create the group
        H5::Group group;
        try {
            group = file_.openGroup(url);
        } catch (const H5::Exception&) {
            group = CreateGroup(file_,url);
        }
        // create the dataset
        int n_dim = img.GetNDim();
        std::vector<hsize_t> dims(n_dim);
        std::reverse_copy(img.GetSize().begin(),img.GetSize().end(),dims.begin());
        H5::DataSpace dspace(n_dim,dims.data());
        H5::DataType dtype(::HDF5Types<T>::Type());
        H5::DataSet dset;
        try {
            dset = group.createDataSet(urn,dtype,dspace);
        } catch (const H5::Exception&) {
            dset = group.openDataSet(urn);
        }
        // write the data in the dataset
        dset.write(img.GetData().data(),dtype,dspace);
    } catch (const H5::FileIException&) {
        return State::HDF5FileException;
    } catch (const H5::GroupIException&) {
        return State::HDF5FileException;
    } catch (const H5::DataSetIException&) {
        return State::HDF5DatasetException;
    } catch (const H5::DataSpaceIException&) {
        return State::HDF5DataspaceException;
    } catch (const H5::DataTypeIException&) {
        return State::HDF5DatatypeException;
    }
    return State::Success;
}

// Template specialisations
// ReadDataset
template State IOh5::ReadDataset<size_t>(Image<size_t> *img, const std::string &url, const std::string &urn);
template State IOh5::ReadDataset<float>(Image<float> *img, const std::string &url, const std::string &urn);
template State IOh5::ReadDataset<double>(Image<double> *img, const std::string &url, const std::string &urn);
template State IOh5::ReadDataset<int>(Image<int> *img, const std::string &url, const std::string &urn);
template State IOh5::ReadDataset<long>(Image<long> *img, const std::string &url, const std::string &urn);
// WriteDataset
template State IOh5::WriteDataset<size_t>(const Image<size_t> &img, const std::string &url, const std::string &urn) const;
template State IOh5::WriteDataset<float>(const Image<float> &img, const std::string &url, const std::string &urn) const;
template State IOh5::WriteDataset<double>(const Image<double> &img, const std::string &url, const std::string &urn) const;
template State IOh5::WriteDataset<int>(const Image<int> &img, const std::string &url, const std::string &urn) const;
template State IOh5::WriteDataset<long>(const Image<long> &img, const std::string &url, const std::string &urn) const;
