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

#include <chrono>
#include <exception>
#include <iostream>
#include <memory>
#include <regex>
#include <utility>
#include <vector>

#include "b1map/io/io_hdf5.h"
#include "b1map/io/io_toml.h"

#include "b1map/b1mapping.h"
#include "b1map/body.h"
#include "b1map/sequences.h"
#include "b1map/version.h"

#include "main.h"

#define LOADMANDATORY(what,io_toml,data,T) { \
    io::IOError MACRO_error = io_toml->what<T>(data.first,data.second); \
    if (MACRO_error!=io::IOError::Success) { \
        string msg = "FATAL ERROR in config file: "+ToString(MACRO_error)+" '"+data.second+"'"; \
        throw runtime_error(msg); \
    } \
}
#define LOADMANDATORYDATA(io_toml,data) LOADMANDATORY(GetValue,io_toml,data,decltype(data.first))
#define LOADMANDATORYLIST(io_toml,data) LOADMANDATORY(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADOPTIONAL(what,io_toml,data,T) { \
    io::IOError MACRO_error = io_toml->what<T>(data.first,data.second); \
    if (WALL && MACRO_error!=io::IOError::Success) { \
        cout<<"WARNING in config file: "+ToString(MACRO_error)+" '"+data.second+"'"<<endl; \
    } \
}
#define LOADOPTIONALDATA(io_toml,data) LOADOPTIONAL(GetValue,io_toml,data,decltype(data.first))
#define LOADOPTIONALLIST(io_toml,data) LOADOPTIONAL(GetArrayOf,io_toml,data,decltype(data.first)::value_type)

#define LOADMAP(map,addr) { \
    string MACRO_fname; \
    string MACRO_uri; \
    io::GetAddress(addr,MACRO_fname,MACRO_uri); \
    io::IOh5 MACRO_ifile(MACRO_fname,io::Mode::In); \
    io::State MACRO_iostate = MACRO_ifile.ReadDataset(&map,"/",MACRO_uri); \
    if (MACRO_iostate!=io::State::Success) { \
        string msg = "FATAL ERROR: "+ToString(MACRO_iostate)+" '"+addr+"'"; \
        throw runtime_error(msg); \
    } \
}
#define SAVEMAP(map,addr) { \
    string MACRO_fname; \
    string MACRO_uri; \
    string MACRO_url; \
    string MACRO_urn; \
    io::GetAddress(addr,MACRO_fname,MACRO_uri); \
    auto const pos = MACRO_uri.find_last_of("/"); \
    MACRO_url = MACRO_uri.substr(0,pos+1); \
    MACRO_urn = MACRO_uri.substr(pos+1); \
    io::IOh5 MACRO_ofile(MACRO_fname,io::Mode::Append); \
    io::State MACRO_iostate = MACRO_ofile.WriteDataset(map,MACRO_url,MACRO_urn); \
    if (MACRO_iostate!=io::State::Success) { \
        string msg = "FATAL ERROR: "+ToString(MACRO_iostate)+" '"+addr+"'"; \
        throw runtime_error(msg); \
    } \
}

using namespace std;
using namespace b1map;

template <class T> using cfgdata = pair<T,string>;
template <class T> using cfglist = pair<array<T,NDIM>,string>;

void SaveComplexMap(Image<complex<double> > img,string addr);

int main(int argc, char **argv) {
    auto start = chrono::system_clock::now();
    // opening boilerplate
    cout<<project::str()<<" ("<<build::str()<<") ["<<compiler::str()<<"]\n"<<endl;
    cout<<LicenseBoilerplate()<<endl;
    // check the number of input
    if (argc<2) {
        cout<<"Usage example: "<<project::str()<<" <config file>"<<endl;
        return -1;
    }
    // load the config file
    unique_ptr<io::IOtoml> io_toml;
    try {
        io_toml.reset(new io::IOtoml(string(argv[1]),io::Mode::In));
    } catch(const ios_base::failure &e) {
        cout<<"FATAL ERROR in config file: "<<e.what()<<endl;
        return 1;
    }
    // declare the input variables
    //   mandatory input
    cfgdata<string> title; title.second = "title";
    cfgdata<int> method; method.second = "method";
    cfglist<int> nn; nn.second = "mesh.size";
    cfglist<double> dd; dd.second = "mesh.step";
    cfgdata<string> body_addr; body_addr.second = "input.body";
    cfgdata<string> txsens_addr; txsens_addr.second = "input.tx-sensitivity";
    cfgdata<string> txphase_addr; txphase_addr.second = "input.tx-phase";
    cfgdata<string> est_addr("","output.alpha-estimate");
    //   optional input
    cfgdata<string> rxsens_addr("","input.rx-sensitivity");
    cfgdata<string> rxphase_addr("","input.rx-phase");
    cfgdata<string> imgs_addr("","output.intermediate-images");
    cfgdata<int> samples(1,"montecarlo.samples");
    cfgdata<double> noise(0.0,"montecarlo.noise");
    // load the input data
    try {
        //   title
        LOADMANDATORYDATA(io_toml,title);
        LOADMANDATORYDATA(io_toml,method);
        //   mesh
        LOADMANDATORYLIST(io_toml,nn);
        LOADMANDATORYLIST(io_toml,dd);
        //   input
        LOADMANDATORYDATA(io_toml,body_addr);
        LOADMANDATORYDATA(io_toml,txsens_addr);
        LOADMANDATORYDATA(io_toml,txphase_addr);
        LOADOPTIONALDATA(io_toml,rxsens_addr);
        LOADOPTIONALDATA(io_toml,rxphase_addr);
        //   output
        LOADMANDATORYDATA(io_toml,est_addr);
        LOADOPTIONALDATA(io_toml,imgs_addr);
        LOADOPTIONALDATA(io_toml,samples);
        LOADOPTIONALDATA(io_toml,noise);
    } catch (const runtime_error &e) {
        cout<<e.what()<<endl;
        return 1;
    }
    cout<<endl;
    // check the provided data
    B1MapMethod b1map_method = static_cast<B1MapMethod>(method.first);
    bool thereis_b1m = (rxsens_addr.first!="" && rxphase_addr.first!="");
    bool thereis_imgs = imgs_addr.first!="";
    bool thereis_noise = noise.first>0;
    //   B1-mapping method
    if (method.first<0||b1map_method>=B1MapMethod::END) {
        cout<<"FATAL ERROR in config file: Wrong data format '"<<method.second<<"'"<<endl;
        return 1;
    }
    //   noise
    if (samples.first<1) {
        cout<<"FATAL ERROR in config file: Wrong data format '"<<samples.second<<"'"<<endl;
        return 1;
    }
    if (!thereis_noise&&samples.first>1) {
        cout<<"WARNING in config file: Without noise the number of samples is set equal to 1"<<endl;
        samples.first = 1;
    }
    // report the readen values
    cout<<"  "<<title.first<<"\n";
    cout<<"\n  Method: ("<<method.first<<") "<<ToString(b1map_method)<<"\n";
    cout<<"\n  Mesh size: ["<<nn.first[0]<<", "<<nn.first[1]<<", "<<nn.first[2]<<"]\n";
    cout<<"  Mesh step: ["<<dd.first[0]<<", "<<dd.first[1]<<", "<<dd.first[2]<<"] m\n";
    cout<<"\n  Number of Monte Carlo samples: "<<samples.first<<"\n";
    cout<<"  Additive noise: "<<noise.first*100.0<<" %\n";
    cout<<"\n  Body details addr.: '"<<body_addr.first<<"'\n";
    cout<<"\n  Tx sensitivity addr.: '"<<txsens_addr.first<<"'\n";
    cout<<"  Tx phase addr.: '"<<txphase_addr.first<<"'\n";
    cout<<"\n  Rx sensitivity addr.: '"<<rxsens_addr.first<<"'\n";
    cout<<"  Rx phase addr.: '"<<rxphase_addr.first<<"'\n";
    cout<<"\n  Output estimate addr.: '"<<est_addr.first<<"'\n";
    cout<<"  Output intermediate images addr.: '"<<imgs_addr.first<<"'\n";
    cout<<endl;
    // load the body details
    Body body;
    try {
        body = Body(body_addr.first);
    } catch (const runtime_error &e) {
        cout<<e.what()<<endl;
        return 1;
    }
    // load b1p and b1m
    Image<complex<double> > b1p(nn.first[0],nn.first[1],nn.first[2]);
    Image<complex<double> > b1m(nn.first[0],nn.first[1],nn.first[2]);
    {
        cout<<"Loading Tx sensitivity and phase:\n"<<flush;
        Image<double> txsens(nn.first[0],nn.first[1],nn.first[2]);
        try {
            LOADMAP(txsens,txsens_addr.first);
        } catch (const runtime_error &e) {
            cout<<e.what()<<endl;
            return 1;
        }
        cout<<"  '"<<txsens_addr.first<<"'\n"<<flush;
        Image<double> txphase(nn.first[0],nn.first[1],nn.first[2]);
        try {
            LOADMAP(txphase,txphase_addr.first);
        } catch (const runtime_error &e) {
            cout<<e.what()<<endl;
            return 1;
        }
        cout<<"  '"<<txphase_addr.first<<"'\n"<<flush;
        for (int idx = 0; idx<b1p.GetNVox(); ++idx) {
            b1p[idx] = txsens[idx]*exp(complex<double>(0.0,txphase[idx]));
        }
    }
    if (thereis_b1m) {
        cout<<"Loading Rx sensitivity and phase:\n"<<flush;
        Image<double> rxsens(nn.first[0],nn.first[1],nn.first[2]);
        try {
            LOADMAP(rxsens,rxsens_addr.first);
        } catch (const runtime_error &e) {
            return 1;
        }
        cout<<"  '"<<rxsens_addr.first<<"'\n"<<flush;
        Image<double> rxphase(nn.first[0],nn.first[1],nn.first[2]);
        try {
            LOADMAP(rxphase,rxphase_addr.first);
        } catch (const runtime_error &e) {
            cout<<e.what()<<endl;
            return 1;
        }
        cout<<"  '"<<rxphase_addr.first<<"'\n"<<flush;
        for (int idx = 0; idx<b1m.GetNVox(); ++idx) {
            b1m[idx] = rxsens[idx]*exp(complex<double>(0.0,rxphase[idx]));
        }
    } else {
        for (int idx = 0; idx<b1m.GetNVox(); ++idx) {
            b1m[idx] = 1.0;
        }
    }
    // load the method parameters and run the method
    Image<double> alpha_est;
    // declare common parameters
    cfgdata<double> alpha_nom; alpha_nom.second="parameter.alpha-nominal";
    cfgdata<double> TR; TR.second="parameter.TR";
    cfgdata<double> TE; TE.second="parameter.TE";
    cfgdata<double> spoiling(1.0,"parameter.spoiling");
    // load the common parameters
    try {
        LOADMANDATORYDATA(io_toml,alpha_nom);
        LOADMANDATORYDATA(io_toml,TR);
        LOADMANDATORYDATA(io_toml,TE);
        LOADOPTIONALDATA(io_toml,spoiling);
    } catch (const runtime_error &e) {
        cout<<e.what()<<endl;
        return 1;
    }
    cout<<endl;
    // check the common parameters
    if (TR.first<0.0) {
        cout<<"FATAL ERROR in config file: Negative '"<<TR.second<<"'"<<endl;
        return 1;
    }
    if (TE.first<0.0) {
        cout<<"FATAL ERROR in config file: Negative '"<<TE.second<<"'"<<endl;
        return 1;
    }
    if (spoiling.first<0.0||spoiling.first>1.0) {
        cout<<"FATAL ERROR in config file: Out of range '"<<spoiling.second<<"'"<<endl;
        return 1;
    }
    // set-up the B1-mapping method
    std::unique_ptr<B1Mapping> b1mapping;
    switch (b1map_method) {
        case B1MapMethod::DA: {
            // report the parameters
            cout<<"Parameters:\n";
            cout<<"  Nominal flip-angle: "<<alpha_nom.first<<" rad\n";
            cout<<"  Repetition time (TR): "<<TR.first<<" ms\n";
            cout<<"  Echo time (TE): "<<TE.first<<" ms\n";
            cout<<"  Spoiling coefficient: "<<spoiling.first<<"\n";
            cout<<endl;
            // initialise the method
            b1mapping.reset(new DoubleAngle(alpha_nom.first,TR.first,TE.first,b1p,b1m,spoiling.first,body));
            break;
        }
        case B1MapMethod::AFI: {
            // declare the special parameters
            cfgdata<double> TRratio; TRratio.second="parameter.TRratio";
            // load the special parameters
            try {
                LOADMANDATORYDATA(io_toml,TRratio);
            } catch (const runtime_error &e) {
                cout<<e.what()<<endl;
                return 1;
            }
            cout<<endl;
            // check the special parameters
            if (TRratio.first<0.0) {
                cout<<"FATAL ERROR in config file: Negative '"<<TRratio.second<<"'"<<endl;
                return 1;
            }
            // report the parameters
            double TR2 = TR.first*TRratio.first;
            cout<<"Parameters:\n";
            cout<<"  Nominal flip-angle: "<<alpha_nom.first<<" rad\n";
            cout<<"  Repetition times (TRs): "<<TR.first<<" ms, "<<TR2<<" ms\n";
            cout<<"  Echo time (TE): "<<TE.first<<" ms\n";
            cout<<"  Spoiling coefficient: "<<spoiling.first<<"\n";
            cout<<endl;
            // initialise the method
            b1mapping.reset(new ActualFlipAngle(alpha_nom.first,TR.first,TR2,TE.first,b1p,b1m,spoiling.first,body));
            break;
        }
        case B1MapMethod::BSS: {
            // declare the special parameters
            cfgdata<double> bss_offres; bss_offres.second="parameter.bss-offres";
            cfgdata<double> bss_length; bss_length.second="parameter.bss-length";
            // load the special parameters
            try {
                LOADMANDATORYDATA(io_toml,bss_offres);
                LOADMANDATORYDATA(io_toml,bss_length);
            } catch (const runtime_error &e) {
                cout<<e.what()<<endl;
                return 1;
            }
            cout<<endl;
            // check the special parameters
            if (bss_offres.first<0.0) {
                cout<<"FATAL ERROR in config file: Negative '"<<bss_offres.second<<"'"<<endl;
                return 1;
            }
            if (bss_length.first<0.0) {
                cout<<"FATAL ERROR in config file: Negative '"<<bss_length.second<<"'"<<endl;
                return 1;
            }
            // report the parameters
            cout<<"Parameters:\n";
            cout<<"  Nominal flip-angle: "<<alpha_nom.first<<" rad\n";
            cout<<"  Repetition time (TR): "<<TR.first<<" ms\n";
            cout<<"  Echo time (TE): "<<TE.first<<" ms\n";
            cout<<"  Spoiling coefficient: "<<spoiling.first<<"\n";
            cout<<"  BSS off-resonance: "<<bss_offres.first<<" kHz\n";
            cout<<"  BSS pulse length: "<<bss_length.first<<" ms\n";
            cout<<endl;
            // initialise the method
            b1mapping.reset(new BlochSiegertShift(alpha_nom.first,TR.first,TE.first,2.0*PI*bss_offres.first,bss_length.first,b1p,b1m,spoiling.first,body));
            break;
        }
    }
    // save the images
    std::array<Image<complex<double> >,2> imgs;
    imgs[0] = b1mapping->GetImg(0);
    imgs[1] = b1mapping->GetImg(1);
    if (thereis_imgs) {
        try {
                SaveComplexMap(imgs[0],imgs_addr.first+"1");
                SaveComplexMap(imgs[1],imgs_addr.first+"2");
        } catch (const runtime_error &e) {
            cout<<e.what()<<endl;
            return 1;
        }
    }
    // apply the method noiseless
    cout<<"Noiseless B1-mapping..."<<flush;
    b1mapping->Run(&alpha_est,0.0);
    // save the result
    try {
        SAVEMAP(alpha_est,est_addr.first);
    } catch (const runtime_error &e) {
        cout<<e.what()<<endl;
        return 1;
    }
    cout<<"done!\n";
    cout<<endl;
    // apply the Monte Carlo with noisy input
    double sigma = ComputeSigma(imgs,noise.first);
    cout<<"Monte Carlo sampling:\n";
    for (int m = 0; m<samples.first; ++m) {
        cout<<"  MC"<<to_string(m)<<"..."<<flush;
        b1mapping->Run(&alpha_est,sigma);
        try {
            SAVEMAP(alpha_est,est_addr.first+"-MC"+to_string(m));
        } catch (const runtime_error &e) {
            cout<<e.what()<<endl;
            return 1;
        }
        cout<<"done!\n";
    }
    cout<<endl;
    //
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::seconds>(end - start);
    cout<<"Execution ended in "<<elapsed.count()<<" s\n";
    cout<<endl;
    return 0;
}

void SaveComplexMap(Image<complex<double> > img,string addr) {
    Image<double> tmp(img.GetSize(0),img.GetSize(1),img.GetSize(2));
    for (int idx = 0; idx<tmp.GetNVox(); ++idx) {
        tmp[idx] = real(img[idx]);
    }
    string real_addr = addr+"/real";
    SAVEMAP(tmp,real_addr);
    for (int idx = 0; idx<tmp.GetNVox(); ++idx) {
        tmp[idx] = imag(img[idx]);
    }
    string imag_addr = addr+"/imag";
    SAVEMAP(tmp,imag_addr)
    return;
}
