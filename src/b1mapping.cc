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
#include <random>

#include "b1map/sequences.h"

namespace b1map {

// B1Mapping constructor
B1Mapping::
B1Mapping() {
	return;
}
// B1Mapping destructor
B1Mapping::
~B1Mapping() {
	return;
}
// B1Mapping GetImg
Image<std::complex<double> >& B1Mapping::
GetImg(const int d) {
	return imgs[d];
}

// DoubleAngle constructor
DoubleAngle::
DoubleAngle(const double alpha_nom, const double TR, const double TE,
	const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	GREImage(&imgs[0],alpha_nom,TR,TE,b1p,b1m,spoiling,body);
	GREImage(&imgs[1],2.0*alpha_nom,TR,TE,b1p,b1m,spoiling,body);
	return;
}
// DoubleAngle destructor
DoubleAngle::
~DoubleAngle() {
	return;
}
// DoubleAngle Run
void DoubleAngle::
Run(Image<double> *alpha_est, const double sigma) {
	*alpha_est = Image<double>(imgs[0].GetSize(0),imgs[0].GetSize(1),imgs[0].GetSize(2));
	std::array<Image<std::complex<double> >,2>* imgs_noise;
	if (sigma > 0.0) {
		imgs_noise = new std::array<Image<std::complex<double> >,2>();
		AddNoise(imgs_noise,imgs,sigma);
	} else {
		imgs_noise = &imgs;
	}
	for (int idx = 0; idx<imgs[0].GetNVox(); ++idx) {
		(*alpha_est)[idx] = std::acos(std::abs((*imgs_noise)[1][idx])/2.0/std::abs((*imgs_noise)[0][idx]));
	}
	if (sigma>0.0) {
		delete imgs_noise;
	}
	return;
}

// ActualFlipAngle constructor
ActualFlipAngle::
ActualFlipAngle(const double alpha_nom, const double TR1, const double TR2,
	const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	AFIImage(&imgs[0],&imgs[1],alpha_nom,TR1,TR2,TE,b1p,b1m,spoiling,body);
	TRratio = TR2/TR1;
	return;
}
// ActualFlipAngle destructor
ActualFlipAngle::
~ActualFlipAngle() {
	return;
}
// ActualFlipAngle Run
void ActualFlipAngle::
Run(Image<double> *alpha_est, const double sigma) {
	*alpha_est = Image<double> (imgs[0].GetSize(0),imgs[0].GetSize(1),imgs[0].GetSize(2));
	std::array<Image<std::complex<double> >,2>* imgs_noise;
	if (sigma > 0.0) {
		imgs_noise = new std::array<Image<std::complex<double> >,2>();
		AddNoise(imgs_noise,imgs,sigma);
	} else {
		imgs_noise = &imgs;
	}
	for (int idx = 0; idx<imgs[0].GetNVox(); ++idx) {
		double tmp = std::abs((*imgs_noise)[1][idx])/std::abs((*imgs_noise)[0][idx]);
		(*alpha_est)[idx] = std::acos((TRratio*tmp-1.0)/(TRratio-tmp));
	}
	if (sigma>0.0) {
		delete imgs_noise;
	}
	return;
}

// BlochSiegertShift constructor
BlochSiegertShift::
BlochSiegertShift(const double alpha_nom, const double TR, const double TE,
	const double bss_offres, const double bss_length,
	const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	BSSImage(&imgs[0],alpha_nom,TR,TE,+bss_offres,bss_length,b1p,b1m,spoiling,body);
	BSSImage(&imgs[1],alpha_nom,TR,TE,-bss_offres,bss_length,b1p,b1m,spoiling,body);
	Kbs = GAMMA*GAMMA*bss_length/2.0/bss_offres;
	return;
}
// BlochSiegertShift destructor
BlochSiegertShift::
~BlochSiegertShift() {
	return;
}
// BlochSiegertShift run
void BlochSiegertShift::
Run(Image<double> *alpha_est, const double sigma) {
	*alpha_est = Image<double>(imgs[0].GetSize(0),imgs[0].GetSize(1),imgs[0].GetSize(2));
	std::array<Image<std::complex<double> >,2>* imgs_noise;
	if (sigma > 0.0) {
		imgs_noise = new std::array<Image<std::complex<double> >,2>();
		AddNoise(imgs_noise,imgs,sigma);
	} else {
		imgs_noise = &imgs;
	}
	for (int idx = 0; idx<imgs[0].GetNVox(); ++idx) {
		(*alpha_est)[idx] = std::sqrt(std::arg((*imgs_noise)[0][idx]/(*imgs_noise)[1][idx])/2.0/Kbs);
	}
	if (sigma>0.0) {
		delete imgs_noise;
	}
	return;
}

// TRxPhaseGRE constructor
TRxPhaseGRE::
TRxPhaseGRE(const double alpha_nom, const double TR, const double TE,
	const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	GREImage(&imgs[0],alpha_nom,TR,TE,b1p,b1m,spoiling,body);
	return;
}
// TRxPhaseGRE destructor
TRxPhaseGRE::
~TRxPhaseGRE() {
	return;
}
// TRxPhaseGRE Run
void TRxPhaseGRE::
Run(Image<double> *alpha_est, const double sigma) {
	*alpha_est = Image<double>(imgs[0].GetSize(0),imgs[0].GetSize(1),imgs[0].GetSize(2));
	Image<std::complex<double> >* img_noise;
	if (sigma > 0.0) {
		img_noise = new Image<std::complex<double> >;
		AddNoise(img_noise,imgs[0],sigma);
	} else {
		img_noise = &(imgs[0]);
	}
	for (int idx = 0; idx<imgs[0].GetNVox(); ++idx) {
		(*alpha_est)[idx] = std::arg((*img_noise)[idx]);
	}
	if (sigma>0.0) {
		delete img_noise;
	}
	return;
}

// Noise utils
double ComputeSigma(const std::array<Image<std::complex<double> >,2> &imgs,
	const double noise) {
	std::array<double,2> sigma{0.0,0.0};
	for (int d = 0; d<2; ++d) {
		Image<double> tmp(imgs[0].GetSize(0),imgs[0].GetSize(1),imgs[0].GetSize(2));
		for (int idx = 0; idx<imgs[d].GetNVox(); ++idx) {
			tmp[idx] = std::abs(imgs[d][idx]);
		}
		sigma[d] = Avg(tmp.GetData())*noise;
	}
	return (sigma[0]+sigma[1])/2.0;
}
void AddNoise(std::array<Image<std::complex<double> >,2> *imgs_noise,
	const std::array<Image<std::complex<double> >,2> &imgs,
	const double sigma) {
	std::random_device generator;
	std::normal_distribution<double> distribution(0.0,sigma);
	for (int d = 0; d<2; ++d) {
		(*imgs_noise)[d] = Image<std::complex<double> >(imgs[d].GetSize(0),imgs[d].GetSize(1),imgs[d].GetSize(2));
		for (int idx = 0; idx<imgs[d].GetNVox(); ++idx) {
			std::complex<double> tmp(distribution(generator),distribution(generator));
			(*imgs_noise)[d][idx] = imgs[d][idx]+tmp;
		}
	}
	return;
}
void AddNoise(Image<std::complex<double> > *img_noise,
	const Image<std::complex<double> > &img,
	const double sigma) {
	std::random_device generator;
	std::normal_distribution<double> distribution(0.0,sigma);
	(*img_noise) = Image<std::complex<double> >(img.GetSize(0),img.GetSize(1),img.GetSize(2));
	for (int idx = 0; idx<img.GetNVox(); ++idx) {
		std::complex<double> tmp(distribution(generator),distribution(generator));
		(*img_noise)[idx] = img[idx]+tmp;
	}
	return;
}

}  // namespace b1map
