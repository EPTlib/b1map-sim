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

#include "b1map/body.h"
#include "b1map/image.h"
#include "b1map/util.h"

namespace b1map {

/**
 * Compute a flip-angle estimate provided by a double-angle method with the
 * provided operative parameters.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Pointer to the complex-valued MRI images destination.
 * @param alpha_nom Nominal flip-angle in radian.
 * @param TR Repetition time in millisecond.
 * @param TE Echo time in millisecond.
 * @param b1p Complex-valued B1+ distribution in tesla.
 * @param b1m Complex-valued B1- distribution.
 * @param spoiling Spoiling coefficient for transverse magnetization:
 *     1 is ideal spoiling; 0 is no spoiling.
 * @param body Physical description of the imaging body.
 */
void DoubleAngle(Image<double> *alpha_est, Image<std::complex<double> > *img1,
	Image<std::complex<double> > *img2, const double alpha_nom, const double TR,
	const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body);
/**
 * Compute a flip-angle estimate provided by a double-angle method applied on
 * the provided images.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Complex-valued MRI images for double-angle b1-mapping.
 */
void DoubleAngleAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2);

/**
 * Compute a flip-angle estimate provided by an actual flip-angle method with
 * the provided operative parameters.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Pointer to the complex-valued MRI images destination.
 * @param alpha_nom Nominal flip-angle in radian.
 * @param TR1,TR2 Repetition times in millisecond.
 * @param TE Echo time in millisecond.
 * @param b1p Complex-valued B1+ distribution in tesla.
 * @param b1m Complex-valued B1- distribution.
 * @param spoiling Spoiling coefficient for transverse magnetization:
 *     1 is ideal spoiling; 0 is no spoiling.
 * @param body Physical description of the imaging body.
 */
void ActualFlipAngle(Image<double> *alpha_est,
	Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
	const double alpha_nom, const double TR1, const double TR2,
	const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body);
/**
 * Compute a flip-angle estimate provided by an actual flip-angle method
 * applied on the provided images.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Complex-valued MRI images for actual flip-angle b1-mapping.
 * @param TRratio Ratio of the repetition times of the two images.
 */
void ActualFlipAngleAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2,
	const double TRratio);

/**
 * Compute a flip-angle estimate provided by a Bloch-Siegert shift method with
 * the provided operative parameters.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Pointer to the complex-valued MRI images destination.
 * @param alpha_nom Nominal flip-angle in radian.
 * @param TR Repetition time in millisecond.
 * @param TE Echo time in millisecond.
 * @param bss_offres Off-resonance frequency of the Bloch-Siegert pulse in
 *     radian per second.
 * @param bss_length Length of the Bloch-Siegert pulse in millisecond.
 * @param b1p Complex-valued B1+ distribution in tesla.
 * @param b1m Complex-valued B1- distribution.
 * @param spoiling Spoiling coefficient for transverse magnetization:
 *     1 is ideal spoiling; 0 is no spoiling.
 * @param body Physical description of the imaging body.
 */
void BlochSiegertShift(Image<double> *alpha_est,
	Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
	const double alpha_nom, const double TR, const double TE,
	const double bss_offres, const double bss_length,
	const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body);
/**
 * Compute a flip-angle estimate provided by a Bloch-Siegert shift method
 * applied on the provided images.
 * 
 * @param alpha_est Pointer to the flip-angle estimate destination.
 * @param img1,img2 Complex-valued MRI images for Bloch-Siegert b1-mapping.
 * @param Kbs Phase coefficient.
 */
void BlochSiegertShiftAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2,
	const double Kbs);

}  // namespace b1map

#endif  // B1MAPSIM_B1MAPPING_H_
