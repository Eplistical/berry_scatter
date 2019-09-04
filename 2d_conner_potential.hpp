#ifndef _2D_CONNER_POTENTIAL_HPP
#define _2D_CONNER_POTENTIAL_HPP

#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"


namespace {
    using std::vector;
    using std::complex;
    using std::pow;

    double param_A = 1.0;
    double param_B = 5.0;
    double param_C = 0.15;
    double param_D = 3.0;
    double param_W = 0.0;

    void output_potential_param() {
        ioer::info("# 2D conner potential [with phi = W * r] parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " C = ", param_C,
                    " D = ", param_D,
                    " W = ", param_W);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 5, 
                "set_potenial_params: potential paramter vector size must be >= 5");
        param_A = params[0];
        param_B = params[1];
        param_C = params[2];
        param_D = params[3];
        param_W = params[4];
    }

    double cal_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        return param_W * pow(x*x + y*y, 0.5);
    }

    vector<double> cal_der_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const double R = pow(x*x + y*y, 0.5);
        vector<double> der_phi(r.size(), 0.0);
        der_phi[0] = param_W * x / R;
        der_phi[1] = param_W * y / R;
        return der_phi;
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > H(4, 0.0);
        H[0+0*2] = tanh(y-param_B) - tanh(y+param_B) + tanh(x) + param_D;
        H[1+1*2] = tanh(x-param_B) - tanh(x+param_B) + tanh(y) + param_D;
        H[0+1*2] = param_C * eip;
        H[1+0*2] = conj(H[0+1*2]);
        return H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const vector<double> der_phi = cal_der_phi(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > nablaHx(4, 0.0);
        nablaHx[0+0*2] = 1.0 - pow(tanh(x), 2);
        nablaHx[1+1*2] = pow(tanh(x+param_B), 2) - pow(tanh(x-param_B), 2);
        nablaHx[0+1*2] = param_C * eip * matrixop::IMAGIZ * der_phi[0];
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        return nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const vector<double> der_phi = cal_der_phi(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > nablaHy(4, 0.0);
        nablaHy[0+0*2] = pow(tanh(y+param_B), 2) - pow(tanh(y-param_B), 2);
        nablaHy[1+1*2] = 1.0 - pow(tanh(y), 2);
        nablaHy[0+1*2] = param_C * eip * matrixop::IMAGIZ * der_phi[1];
        nablaHy[1+0*2] = conj(nablaHy[0+1*2]);
        return nablaHy;
    }

    void cal_info_nume(const vector<double>& r, 
            vector< complex<double> >& Fx, vector< complex<double> >& Fy, 
            vector< complex<double> >& dcx, vector< complex<double> >& dcy, 
            vector<double>& eva, vector< complex<double> >& lastevt 
            )
    {
        vector< complex<double> > evt;
        matrixop::hdiag(cal_H(r), eva, evt);
        // correct phase
        if (not lastevt.empty()) {
            auto tmp = matrixop::matCmat(lastevt, evt, 2);
            for (int j = 0; j < 2; ++j) {
                complex<double> eip = tmp[j+j*2] / abs(tmp[j+j*2]);
                for (int k = 0; k < 2; ++k) {
                    evt[k+j*2] /= eip;
                }
            }
        }
        // F, dc
        dcx = matrixop::matCmatmat(evt, cal_nablaHx(r), evt, 2, 2);
        dcy = matrixop::matCmatmat(evt, cal_nablaHy(r), evt, 2, 2);
        Fx.assign(4, 0.0);
        Fy.assign(4, 0.0);
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                Fx[j+k*2] = -dcx[j+k*2];
                Fy[j+k*2] = -dcy[j+k*2];
                if (j == k) {
                    dcx[j+k*2] = 0.0;
                    dcy[j+k*2] = 0.0;
                }
                else {
                    dcx[j+k*2] /= (eva[k] - eva[j]);
                    dcy[j+k*2] /= (eva[k] - eva[j]);
                }
            }
        }
        lastevt = move(evt);
    }
};

#endif
