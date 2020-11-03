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

#include "b1map/b1mapping.h"

#include <iostream>

#include "b1map/sequences.h"

namespace b1map {

// Double-angle b1-mapping
void DoubleAngle(Image<double> *alpha_est, Image<std::complex<double> > *img1,
	Image<std::complex<double> > *img2, const double alpha_nom, const double TR,
	const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	// noiseless GRE images
	GREImage(img1,alpha_nom,TR,TE,b1p,b1m,spoiling,body);
	GREImage(img2,2.0*alpha_nom,TR,TE,b1p,b1m,spoiling,body);
	// estimate the flip-angle
	DoubleAngleAlpha(alpha_est,*img1,*img2);
	return;
}
void DoubleAngleAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2) {
	*alpha_est = Image<double>(img1.GetSize(0),img1.GetSize(1),img1.GetSize(2));
	for (int idx = 0; idx<img1.GetNVox(); ++idx) {
		(*alpha_est)[idx] = std::acos(std::abs(img2[idx])/2.0/std::abs(img1[idx]));
	}
	return;
}

// Actual flip-angle b1-mapping
void ActualFlipAngle(Image<double> *alpha_est,
	Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
	const double alpha_nom, const double TR1, const double TR2,
	const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	// noiseless AFI images
	AFIImage(img1,img2,alpha_nom,TR1,TR2,TE,b1p,b1m,spoiling,body);
	// estimate the flip-angle
	ActualFlipAngleAlpha(alpha_est,*img1,*img2,TR2/TR1);
	return;
}
void ActualFlipAngleAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2,
	const double TRratio) {
	*alpha_est = Image<double>(img1.GetSize(0),img1.GetSize(1),img1.GetSize(2));
	for (int idx = 0; idx<img1.GetNVox(); ++idx) {
		double tmp = std::abs(img2[idx])/std::abs(img1[idx]);
		(*alpha_est)[idx] = std::acos((TRratio*tmp-1.0)/(TRratio-tmp));
	}
	return;
}

// Bloch-Siegert shift b1-mapping
void BlochSiegertShift(Image<double> *alpha_est,
	Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
	const double alpha_nom, const double TR, const double TE,
	const double bss_offres, const double bss_length,
	const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	// noiseless BSS images
	BSSImage(img1,alpha_nom,TR,TE,+bss_offres,bss_length,b1p,b1m,spoiling,body);
	BSSImage(img2,alpha_nom,TR,TE,-bss_offres,bss_length,b1p,b1m,spoiling,body);
	// estimate the flip-angle
	double Kbs = GAMMA*GAMMA*bss_length/2.0/bss_offres;
	BlochSiegertShiftAlpha(alpha_est,*img1,*img2,Kbs);
	return;
}
void BlochSiegertShiftAlpha(Image<double> *alpha_est,
	const Image<std::complex<double> > &img1,
	const Image<std::complex<double> > &img2, const double Kbs) {
	*alpha_est = Image<double>(img1.GetSize(0),img1.GetSize(1),img1.GetSize(2));
	for (int idx = 0; idx<img1.GetNVox(); ++idx) {
		(*alpha_est)[idx] = std::sqrt(std::arg(img1[idx]/img2[idx])/2.0/Kbs);
	}
	return;
}

}  // namespace b1map
