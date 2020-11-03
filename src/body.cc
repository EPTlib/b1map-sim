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

#include "b1map/body.h"

#include "b1map/io/io_hdf5.h"
#include "b1map/io/io_toml.h"

namespace b1map {

void ReadConfig(const io::IOtoml &file, std::pair<std::string,std::string> *arg);
template <typename T> void LoadMap(Image<T> *map, const std::string &address);

// Body constructors
Body::
Body() :
	materials_(), t1_(), t2star_() {
	return;
}
Body::
Body(const std::string &fname) {
	std::pair<std::string,std::string> mat_addr; mat_addr.second = "body.materials";
	std::pair<std::string,std::string> rho_addr; rho_addr.second = "body.proton-density";
	std::pair<std::string,std::string> t1_addr; t1_addr.second = "body.longitudinal-relaxation";
	std::pair<std::string,std::string> t2star_addr; t2star_addr.second = "body.transverse-relaxation";
	// read the configuration file
	io::IOtoml io_toml(fname,io::Mode::In);
	ReadConfig(io_toml,&mat_addr);
	ReadConfig(io_toml,&rho_addr);
	ReadConfig(io_toml,&t1_addr);
	ReadConfig(io_toml,&t2star_addr);
	// load the maps
	LoadMap(&materials_,mat_addr.first);
	LoadMap(&rho_,rho_addr.first);
	LoadMap(&t1_,t1_addr.first);
	LoadMap(&t2star_,t2star_addr.first);
	return;
}

// Getters
Image<int>& Body::
GetMaterials() {
	return materials_;
}
const Image<int>& Body::
GetMaterials() const {
	return materials_;
}
Image<double>& Body::
GetRho() {
	return rho_;
}
const Image<double>& Body::
GetRho() const {
	return rho_;
}
Image<double>& Body::
GetT1() {
	return t1_;
}
const Image<double>& Body::
GetT1() const {
	return t1_;
}
Image<double>& Body::
GetT2Star() {
	return t2star_;
}
const Image<double>& Body::
GetT2Star() const {
	return t2star_;
}

// Read an argumento from the configuration file
void ReadConfig(const io::IOtoml &file, std::pair<std::string,std::string> *arg) {
	io::IOError error = file.GetValue<std::string>(arg->first,arg->second);
	if (error!=io::IOError::Success) {
		throw std::runtime_error(io::ToString(error)+" '"+arg->second+"'");
	}
	return;
}

// Load a map from the h5 file
template <typename T>
void LoadMap(Image<T> *map, const std::string &address) {
	std::string fname;
	std::string uri;
	io::GetAddress(address,fname,uri);
	io::IOh5 io_h5(fname,io::Mode::In);
	io::State state = io_h5.ReadDataset(map,"/",uri);
	if (state!=io::State::Success) {
		throw std::runtime_error(io::ToString(state)+" '"+address+"'");
	}
	return;
}

}  // namespace b1map
