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

#ifndef B1MAPSIM_BODY_H_
#define B1MAPSIM_BODY_H_

#include "b1map/body.h"

#include <string>

#include "b1map/image.h"

namespace b1map {

/**
 * Class for the imaged body description.
 */
class Body {
	public:
		/**
		 * Default constructor.
		 */
		Body();
		/**
		 * Standard constructor.
		 * 
		 * @param fname Address of the .toml file to read.
		 */
		Body(const std::string &fname);
		/**
		 * Get a reference to the material codes.
		 * 
		 * @return a reference to the materials image.
		 */
		Image<int>& GetMaterials();
		/**
		 * Get a constant reference to the material codes.
		 * 
		 * @return a constant reference to the materials image.
		 */
		const Image<int>& GetMaterials() const;
		/**
		 * Get a reference to the proton density list.
		 * 
		 * @return a reference to the proton density list.
		 */
		Image<double>& GetRho();
		/**
		 * Get a constant reference to the proton density list.
		 * 
		 * @return a constant reference to the proton density list.
		 */
		const Image<double>& GetRho() const;
		/**
		 * Get a reference to the T1 list.
		 * 
		 * @return a reference to the T1 list.
		 */
		Image<double>& GetT1();
		/**
		 * Get a constant reference to the T1 list.
		 * 
		 * @return a constant reference to the T1 list.
		 */
		const Image<double>& GetT1() const;
		/**
		 * Get a reference to the T2star list.
		 * 
		 * @return a reference to the T2star list.
		 */
		Image<double>& GetT2Star();
		/**
		 * Get a constant reference to the T2star list.
		 * 
		 * @return a constant reference to the T2star list.
		 */
		const Image<double>& GetT2Star() const;
	private:
		/// 3D image of the material codes.
		Image<int> materials_;
		/// List of the proton densities.
		Image<double> rho_;
		/// List of the longitudinal relaxation times.
		Image<double> t1_;
		/// List of the transverse relaxation times.
		Image<double> t2star_;
};

}  // namespace b1map

#endif  // B1MAPSIM_BODY_H_
