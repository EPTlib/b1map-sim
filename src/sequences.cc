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

#include "b1map/sequences.h"

#include <iostream>

namespace b1map {

// GRE image
void GREImage(Image<std::complex<double> > *img, const double alpha_nom,
	const double TR, const double TE, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	// initialize the result
	*img = Image<std::complex<double> >(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	// shortcut variables
	const Image<int> &mat = body.GetMaterials();
	const Image<double> &rho = body.GetRho();
	const Image<double> &t1 = body.GetT1();
	const Image<double> &t2star = body.GetT2Star();
	// actual flip-angle distribution
	Image<double> alpha;
	EvalAlpha(&alpha,b1p,alpha_nom);
	// relaxation coefficients
	int n_mat = body.GetT1().GetNVox();
	Image<double> e1(n_mat);
	Image<double> e2(n_mat);
	for (int id_mat = 0; id_mat<n_mat; ++id_mat) {
		e1[id_mat] = std::exp(-TR/t1[id_mat]);
		e2[id_mat] = std::exp(-TR/t2star[id_mat])*(1.0-spoiling);
	}
	// solve Bloch equations
	Image<double> my_old(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> mx(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> my(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> mz(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	std::fill(my_old.GetData().begin(),my_old.GetData().end(),0.0);
	std::fill(mx.GetData().begin(),mx.GetData().end(),0.0);
	std::fill(my.GetData().begin(),my.GetData().end(),0.0);
	std::fill(mz.GetData().begin(),mz.GetData().end(),1.0);
	while (true) {
		RFPulse(&mx,&my,&mz,alpha);
		if (IsSteadyState(my,my_old)) {
			break;
		}
		std::copy(my.GetData().begin(),my.GetData().end(),my_old.GetData().begin());
		Relax(&mx,&my,&mz,e1,e2,mat);
	}
	// synthesize the image
	for (int idx = 0; idx<img->GetNVox(); ++idx) {
		(*img)[idx] = std::complex<double>(mx[idx],my[idx]) *
			rho[mat[idx]] *
			std::exp(-TE/t2star[mat[idx]]) *
			b1m[idx]*b1p[idx]/std::abs(b1p[idx]);
	}
	return;
}

// AFI image
void AFIImage(Image<std::complex<double> > *img1, Image<std::complex<double> > *img2,
 	const double alpha_nom, const double TR1, const double TR2, const double TE,
	const Image<std::complex<double> > &b1p, const Image<std::complex<double> > &b1m,
	const double spoiling, const Body &body) {
	// initialize the result
	*img1 = Image<std::complex<double> >(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	*img2 = Image<std::complex<double> >(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	// shortcut variables
	const Image<int> &mat = body.GetMaterials();
	const Image<double> &rho = body.GetRho();
	const Image<double> &t1 = body.GetT1();
	const Image<double> &t2star = body.GetT2Star();
	// actual flip-angle distribution
	Image<double> alpha;
	EvalAlpha(&alpha,b1p,alpha_nom);
	// relaxation coefficients
	int n_mat = body.GetT1().GetNVox();
	Image<double> e11(n_mat);
	Image<double> e12(n_mat);
	Image<double> e21(n_mat);
	Image<double> e22(n_mat);
	for (int id_mat = 0; id_mat<n_mat; ++id_mat) {
		e11[id_mat] = std::exp(-TR1/t1[id_mat]);
		e12[id_mat] = std::exp(-TR1/t2star[id_mat])*(1.0-spoiling);
		e21[id_mat] = std::exp(-TR2/t1[id_mat]);
		e22[id_mat] = std::exp(-TR2/t2star[id_mat])*(1.0-spoiling);
	}
	// solve Bloch equations
	Image<double> mx(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> my(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> mz(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> m1(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> m2(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> m1_old(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> m2_old(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	std::fill(mx.GetData().begin(),mx.GetData().end(),0.0);
	std::fill(my.GetData().begin(),my.GetData().end(),0.0);
	std::fill(mz.GetData().begin(),mz.GetData().end(),1.0);
	std::fill(m1.GetData().begin(),m1.GetData().end(),0.0);
	std::fill(m2.GetData().begin(),m2.GetData().end(),0.0);
	std::fill(m1_old.GetData().begin(),m1_old.GetData().end(),0.0);
	std::fill(m2_old.GetData().begin(),m2_old.GetData().end(),0.0);
	while (true) {
		RFPulse(&mx,&my,&mz,alpha);
		std::copy(my.GetData().begin(),my.GetData().end(),m1.GetData().begin());
		Relax(&mx,&my,&mz,e11,e12,mat);
		RFPulse(&mx,&my,&mz,alpha);
		std::copy(my.GetData().begin(),my.GetData().end(),m2.GetData().begin());
		if (IsSteadyState(m1,m1_old) && IsSteadyState(m2,m2_old)) {
			break;
		}
		std::copy(m1.GetData().begin(),m1.GetData().end(),m1_old.GetData().begin());
		std::copy(m2.GetData().begin(),m2.GetData().end(),m2_old.GetData().begin());
		Relax(&mx,&my,&mz,e21,e22,mat);
	}
	// synthesize the image
	for (int idx = 0; idx<img1->GetNVox(); ++idx) {
		std::complex<double> tmp = rho[mat[idx]] *
			std::exp(-TE/t2star[mat[idx]]) *
			b1m[idx]*b1p[idx]/std::abs(b1p[idx]);
		(*img1)[idx] = std::complex<double>(0.0,m1[idx])*tmp;
		(*img2)[idx] = std::complex<double>(0.0,m2[idx])*tmp;
	}
	return;
}

// BSS image
void BSSImage(Image<std::complex<double> > *img, const double alpha_nom,
	const double TR, const double TE, const double bss_offres,
	const double bss_length, const Image<std::complex<double> > &b1p,
	const Image<std::complex<double> > &b1m, const double spoiling,
	const Body &body) {
	// initialize the result
	*img = Image<std::complex<double> >(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	// shortcut variables
	const Image<int> &mat = body.GetMaterials();
	const Image<double> &rho = body.GetRho();
	const Image<double> &t1 = body.GetT1();
	const Image<double> &t2star = body.GetT2Star();
	// actual flip-angle distribution
	Image<double> alpha;
	EvalAlpha(&alpha,b1p,alpha_nom);
	// relaxation coefficients
	int n_mat = body.GetT1().GetNVox();
	Image<double> e1(n_mat);
	Image<double> e2(n_mat);
	for (int id_mat = 0; id_mat<n_mat; ++id_mat) {
		e1[id_mat] = std::exp(-TR/t1[id_mat]);
		e2[id_mat] = std::exp(-TR/t2star[id_mat])*(1.0-spoiling);
	}
	// Bloch-Siegert angle
	double bss_angle = bss_offres*bss_length;
	Image<double> phi(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	for (int idx = 0; idx<b1p.GetNVox(); ++idx) {
		phi[idx] = std::sqrt(GAMMA*std::abs(b1p[idx])*GAMMA*std::abs(b1p[idx]) + bss_offres*bss_offres)*bss_length;
	}
	// solve Bloch equations
	Image<double> mx_old(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> my_old(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> mx(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> my(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	Image<double> mz(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	std::fill(mx_old.GetData().begin(),mx_old.GetData().end(),0.0);
	std::fill(my_old.GetData().begin(),my_old.GetData().end(),0.0);
	std::fill(mx.GetData().begin(),mx.GetData().end(),0.0);
	std::fill(my.GetData().begin(),my.GetData().end(),0.0);
	std::fill(mz.GetData().begin(),mz.GetData().end(),1.0);
	while (true) {
		RFPulse(&mx,&my,&mz,alpha);
		for (int idx = 0; idx<mx.GetNVox(); ++idx) {
			double angle_x = GAMMA*std::abs(b1p[idx])*bss_length;
			double tmpx = -bss_angle/phi[idx]*my[idx]*std::sin(phi[idx]/2.0);
			double tmpy = (bss_angle*mx[idx]-angle_x*mz[idx])/phi[idx]*std::sin(phi[idx]/2.0);
			double tmpz = angle_x*my[idx]/phi[idx]*std::sin(phi[idx]/2.0);
			mx[idx] += 2.0*std::cos(phi[idx]/2.0)*tmpx - 2.0*bss_angle/phi[idx]*tmpy*std::sin(phi[idx]/2.0);
			my[idx] += 2.0*std::cos(phi[idx]/2.0)*tmpy + 2.0*(bss_angle*tmpx - angle_x*tmpz)/phi[idx]*std::sin(phi[idx]/2.0);
			mz[idx] += 2.0*std::cos(phi[idx]/2.0)*tmpz + 2.0*angle_x/phi[idx]*tmpy*std::sin(phi[idx]/2.0);
			tmpx = mx[idx];
			mx[idx] =   std::cos(bss_angle)*tmpx + std::sin(bss_angle)*my[idx];
			my[idx] = - std::sin(bss_angle)*tmpx + std::cos(bss_angle)*my[idx];
		}
		if (IsSteadyState(mx,my,mx_old,my_old)) {
			break;
		}
		std::copy(mx.GetData().begin(),mx.GetData().end(),mx_old.GetData().begin());
		std::copy(my.GetData().begin(),my.GetData().end(),my_old.GetData().begin());
		Relax(&mx,&my,&mz,e1,e2,mat);
	}
	// synthesize the image
	for (int idx = 0; idx<img->GetNVox(); ++idx) {
		(*img)[idx] = std::complex<double>(mx[idx],my[idx]) *
			rho[mat[idx]] *
			std::exp(-TE/t2star[mat[idx]]) *
			b1m[idx]*b1p[idx]/std::abs(b1p[idx]);
	}
	return;
}

// Evaluate the actual flip-angle
void EvalAlpha(Image<double> *alpha, const Image<std::complex<double> > &b1p, const double alpha_nom) {
	*alpha = Image<double>(b1p.GetSize(0),b1p.GetSize(1),b1p.GetSize(2));
	for (int idx = 0; idx<alpha->GetNVox(); ++idx) {
		(*alpha)[idx] = std::abs(b1p[idx]);
	}
	double avg = Avg(alpha->GetData());
	for (int idx = 0; idx<alpha->GetNVox(); ++idx) {
		(*alpha)[idx] *= alpha_nom/avg;
	}
	return;
}

// RF pulse
void RFPulse(Image<double> *mx, Image<double> *my, Image<double> *mz,
	const Image<double> &alpha) {
	for (int idx = 0; idx<mx->GetNVox(); ++idx) {
		double tmp = (*my)[idx];
		(*my)[idx] = std::cos(alpha[idx])*tmp - std::sin(alpha[idx])*(*mz)[idx];
		(*mz)[idx] = std::sin(alpha[idx])*tmp + std::cos(alpha[idx])*(*mz)[idx];
	}
	return;
}

// Relaxation
void Relax(Image<double> *mx, Image<double> *my, Image<double> *mz,
	const Image<double> &e1, const Image<double> &e2, const Image<int> &mat) {
	for (int idx = 0; idx<mx->GetNVox(); ++idx) {
		(*mx)[idx] *= e2[mat[idx]];
		(*my)[idx] *= e2[mat[idx]];
		(*mz)[idx] = 1.0 + e1[mat[idx]]*((*mz)[idx]-1.0);
	}
	return;
}

// Is steady-state?
bool IsSteadyState(const Image<double> &mx, const Image<double> &my,
	const Image<double> &mx_old, const Image<double> &my_old) {
	double res = DiffNorm2(mx.GetData(),mx_old.GetData()) +
		DiffNorm2(my.GetData(),my_old.GetData());
	double ref = Norm2(mx_old.GetData()) + Norm2(my_old.GetData());
	bool isss = std::sqrt(res) < 1e-10*std::sqrt(ref);
	return isss;
}
bool IsSteadyState(const Image<double> &m, const Image<double> &m_old) {
	double res = DiffNorm2(m.GetData(),m_old.GetData());
	double ref = Norm2(m_old.GetData());
	bool isss = (std::sqrt(res) < 1e-10*std::sqrt(ref));
	return isss;
}

}  // namespace b1map
