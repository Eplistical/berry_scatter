#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/fmtstring.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/matrixop.hpp"

#include "2d_conner_potential.hpp"

vector< complex<double> > lastevt;
vector<double> eva(2);
vector< complex<double> > Fx(4), Fy(4);
vector< complex<double> > dcx(4), dcy(4);
vector< complex<double> > dx_dcx(4), dy_dcy(4);

int main(int argc, char** argv) {
    output_potential_param();

    const int N = 128;
    const double xmin = -15.0, xmax = 15.0;
    vector<double> xarr = linspace(xmin, xmax, N);
    vector<double> yarr = linspace(xmin, xmax, N);

    ioer::tabout("# x", "y", "E0", "E1", "H00", "H11", "FBx", "FBy", "Re(dcx01)", "Im(dcx01)", "Re(dcx10)", "Im(dcx10)");

    for (int i(0); i < N; ++i) {
        const double x = xarr[i];
        for (int j(0); j < N; ++j) {
            const double y = yarr[j];
            vector<double> r { x, y };

            vector<double> v { 0.02, -0.02 };

            auto H = cal_H(r);

            // calc info
            cal_info_nume(r, Fx, Fy, dcx, dcy, eva, lastevt);

            double Fx_berry = 0.0, Fy_berry = 0.0;
            const int s = 0;
            Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
            Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();

            ioer::tabout(x, y, 
                    eva[0], eva[1], 
                    H[0+0*2].real(), H[1+1*2].real(),
                    Fx_berry,
                    Fy_berry,
                    dcx[0+1*2].real(),
                    dcx[0+1*2].imag(),
                    dcx[1+0*2].real(),
                    dcx[1+0*2].imag()
                    );
        }
    }
    return 0;
}
