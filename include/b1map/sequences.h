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

#ifndef B1MAPSIM_SEQUENCES_H_
#define B1MAPSIM_SEQUENCES_H_

#include "b1map/body.h"
#include "b1map/image.h"
#include "b1map/util.h"

namespace b1map {

/**
 * Generate a complex-valued MRI image acquired by a GRE sequence with the
 * provided operative parameters.
 * 
 * @param img Pointer to the image destination.
 * @param alpha_nom Nominal flip-angle in radian.
 * @param TR Repetition time in millisecond.
 * @param TE Echo time in millisecond.
 * @param b1p Complex-valued B1+ distribution in tesla.
 * @param b1m Complex-valued B1- distribution.
 * @param spoiling Spoiling coefficient for transverse magnetization:
 *     1 is ideal spoiling; 0 is no spoiling.
 * @param body Physical description of the imaging body.
 */
void GREImage(Image<std::complex<double> > *img, const double alpha_nom,
	const double TR, const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body);

/**
 * Generate two complex-valued MRI images acquired by interleaved GRE sequences
 * with the provided operative parameters, as used by the actual flip-angle
 * b1-mapping method.
 * 
 * @param img1,img2 Pointers to the image destinations.
 * @param alpha_nom Nominal flip-angle in radian.
 * @param TR1,TR2 Repetition times in millisecond.
 * @param TE Echo time in millisecond.
 * @param b1p Complex-valued B1+ distribution in tesla.
 * @param b1m Complex-valued B1- distribution.
 * @param spoiling Spoiling coefficient for transverse magnetization:
 *     1 is ideal spoiling; 0 is no spoiling.
 * @param body Physical description of the imaging body.
 */
void AFIImage(Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
	const double alpha_nom, const double TR1, const double TR2, const double TE,
	const Image<std::complex<double> > &b1p, const Image<std::complex<double> > &b1m,
	const double spoiling, const Body &body);

/**
 * Generate a complex-valued MRI images acquired by a GRE sequence with the
 * provided operative parameters, as used by the Bloch-Siegert shift
 * b1-mapping method.
 * 
 * @param img Pointer to the image destination.
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
void BSSImage(Image<std::complex<double> > *img, const double alpha_nom,
	const double TR, const double TE, const double bss_offres,
	const double bss_length, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body);

/**
 * Evaluate the actual flip-angle distribution.
 * 
 * @param alpha Pointer to the destination.
 * @param b1p Complex-valued B1+ distribution.
 * @param alpha_nom Nominal flip-angle in radian.
 */
void EvalAlpha(Image<double> *alpha, const Image<std::complex<double> > &b1p,
	const double alpha_nom);

/**
 * Apply the rotation due to an RF pulse with given flip-angle.
 * 
 * @param mx,my,mz Pointers to the magnetization components.
 * @param alpha Actual flip-angle distribution.
 */
void RFPulse(Image<double> *mx, Image<double> *my, Image<double> *mz,
	const Image<double> &alpha);

/**
 * Apply the spin relaxations.
 * 
 * @param mx,my,mz Pointers to the magnetization components.
 * @param e1 Longitudinal relaxation coefficients.
 * @param e2 Transverse relaxation coefficients.
 * @param mat Material codes.
 */
void Relax(Image<double> *mx, Image<double> *my, Image<double> *mz,
	const Image<double> &e1, const Image<double> &e2, const Image<int> &mat);

/**
 * 
 */
bool IsSteadyState(const Image<double> &mx, const Image<double> &my,
	const Image<double> &mx_old, const Image<double> &my_old);
/**
 * 
 */
bool IsSteadyState(const Image<double> &m, const Image<double> &m_old);

}  // namespace b1map

#endif  // B1MAPSIM_SEQUENCES_H_
