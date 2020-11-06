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

#ifndef B1MAPSIM_B1MAPPING_H_
#define B1MAPSIM_B1MAPPING_H_

#include <array>

#include "b1map/body.h"
#include "b1map/image.h"
#include "b1map/util.h"

namespace b1map {

/**
 * Abstract interface representing any B1-mapping method.
 */
class B1Mapping {
    public:
        /**
         * Constructor.
         */
        B1Mapping();
        /**
         * Virtual destructor.
         */
        virtual ~B1Mapping() = 0;
        /**
         * Abstract method performing the b1-mapping.
		 * 
		 * @param alpha_est Pointer to the flip-angle estimate destination.
		 * @param sigma Standard deviation of the noise in the images.
         */
        virtual void Run(Image<double> *alpha_est, const double sigma) = 0;
		/**
		 * 
		 */
		Image<std::complex<double> >& GetImg(const int d);
	protected:
		/// Complex-valued MRI images.
		std::array<Image<std::complex<double> >,2> imgs;
};

/**
 * Implementation of the double-angle B1-mapping method.
 */
class DoubleAngle : public B1Mapping {
    public:
        /**
         * Constructor.
		 *
		 * @param alpha_nom Nominal flip-angle in radian.
		 * @param TR Repetition time in millisecond.
		 * @param TE Echo time in millisecond.
		 * @param b1p Complex-valued B1+ distribution in tesla.
		 * @param b1m Complex-valued B1- distribution.
		 * @param spoiling Spoiling coefficient for transverse magnetization:
		 *     1 is ideal spoiling; 0 is no spoiling.
		 * @param body Physical description of the imaging body.
         */
        DoubleAngle(const double alpha_nom, const double TR, const double TE,
			const Image<std::complex<double> > &b1p,
			const Image<std::complex<double> > &b1m, const double spoiling,
			const Body &body);
        /**
         * Virtual destructor.
         */
        virtual ~DoubleAngle();
        /**
         * Abstract method performing the b1-mapping.
		 * 
		 * @param alpha_est Pointer to the flip-angle estimate destination.
		 * @param sigma Standard deviation of the noise in the images.
         */
        virtual void Run(Image<double> *alpha_est, const double sigma);
};

/**
 * Implementation of the actual flip-angle B1-mapping method.
 */
class ActualFlipAngle : public B1Mapping {
    public:
        /**
         * Constructor.
		 *
		 * @param alpha_nom Nominal flip-angle in radian.
		 * @param TR1,TR2 Repetition times in millisecond.
		 * @param TE Echo time in millisecond.
		 * @param b1p Complex-valued B1+ distribution in tesla.
		 * @param b1m Complex-valued B1- distribution.
		 * @param spoiling Spoiling coefficient for transverse magnetization:
		 *     1 is ideal spoiling; 0 is no spoiling.
		 * @param body Physical description of the imaging body.
         */
        ActualFlipAngle(const double alpha_nom, const double TR1,
			const double TR2, const double TE,
			const Image<std::complex<double> > &b1p,
			const Image<std::complex<double> > &b1m, const double spoiling,
			const Body &body);
        /**
         * Virtual destructor.
         */
        virtual ~ActualFlipAngle();
        /**
         * Abstract method performing the b1-mapping.
		 * 
		 * @param alpha_est Pointer to the flip-angle estimate destination.
		 * @param sigma Standard deviation of the noise in the images.
         */
        virtual void Run(Image<double> *alpha_est, const double sigma);
	private:
		/// Ratio of the repetition times
		double TRratio;
};

/**
 * Implementation of the Bloch-Siegert shift B1-mapping method.
 */
class BlochSiegertShift : public B1Mapping {
    public:
        /**
         * Constructor.
		 *
		 * @param alpha_nom Nominal flip-angle in radian.
		 * @param TR Repetition time in millisecond.
		 * @param TE Echo time in millisecond.
		 * @param bss_offres Off-resonance frequency of the Bloch-Siegert pulse in
    	 *     radian per millisecond.
    	 * @param bss_length Length of the Bloch-Siegert pulse in millisecond.
		 * @param b1p Complex-valued B1+ distribution in tesla.
		 * @param b1m Complex-valued B1- distribution.
		 * @param spoiling Spoiling coefficient for transverse magnetization:
		 *     1 is ideal spoiling; 0 is no spoiling.
		 * @param body Physical description of the imaging body.
         */
        BlochSiegertShift(const double alpha_nom, const double TR,
			const double TE, const double bss_offres, const double bss_length,
			const Image<std::complex<double> > &b1p,
			const Image<std::complex<double> > &b1m, const double spoiling,
			const Body &body);
        /**
         * Virtual destructor.
         */
        virtual ~BlochSiegertShift();
        /**
         * Abstract method performing the b1-mapping.
		 * 
		 * @param alpha_est Pointer to the flip-angle estimate destination.
		 * @param sigma Standard deviation of the noise in the images.
         */
        virtual void Run(Image<double> *alpha_est, const double sigma);
	private:
		/// 
		double Kbs;
};

/**
 * 
 */
double ComputeSigma(const std::array<Image<std::complex<double> >,2> &imgs,
	const double noise);

/**
 * 
 */
void AddNoise(std::array<Image<std::complex<double> >,2> *imgs_noise,
	const std::array<Image<std::complex<double> >,2> &imgs,
	const double sigma);

}  // namespace b1map

#endif  // B1MAPSIM_B1MAPPING_H_
