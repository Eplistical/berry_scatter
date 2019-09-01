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
#include "misc/fft.hpp"
#include "misc/matrixop.hpp"
#include "boost/program_options.hpp"

//#include "2d_flat_potential.hpp"
#include "2d_conner_potential.hpp"

using namespace std;
namespace po = boost::program_options;

double L = 16;
int M = 256;

vector<double> potential_params;
double mass = 1000.0;

double xI = -3.0;
double yI = 0.0;
double sigmax = 1.0;
double sigmay = 1.0;
double kxI = 20.0;
double kyI = 0.0;

double init_s = 0.0;
int Nstep = 4000;
double dt = 0.1;

int output_step = 250;

bool enable_adiab = false;
bool enable_abc = false;

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("L", po::value<double>(&L), "grid para, the grid is [-L/2, L/2]")
        ("M", po::value<int>(&M), "grid number")
        ("mass", po::value<double>(&mass), "mass")
        ("init_x", po::value<double>(&xI), "init x")
        ("init_y", po::value<double>(&yI), "init y")
        ("sigma_x", po::value<double>(&sigmax), "init sigma x")
        ("sigma_y", po::value<double>(&sigmay), "init sigma y")
        ("init_px", po::value<double>(&kxI), "init px")
        ("init_py", po::value<double>(&kyI), "init py")
        ("init_s", po::value<double>(&init_s), "init adiab surface")
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("dt", po::value<double>(&dt), "time step")
        ("enable_adiab", po::value<bool>(&enable_adiab), "output observables on adiabats")
        ("enable_abc", po::value<bool>(&enable_abc), "absorbing boundary condition")
        ("output_step", po::value<int>(&output_step), "output step")
        ;
    po::variables_map vm; 
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    if (not potential_params.empty()) {
        set_potenial_params(potential_params);
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

vector< complex<double> > myfftshift(const vector< complex<double> >& Z) 
{
    vector< complex<double> > rst(M * M);
    for (int k = 0; k < M/2; ++k) {
        for (int j = 0; j < M/2; ++j) {
            rst[k+j*M] = Z[(k+M/2) + (j+M/2)*M];
            rst[(k+M/2)+(j+M/2)*M] = Z[k + j*M];
            rst[(k+M/2)+j*M] = Z[k + (j+M/2)*M];
            rst[k+(j+M/2)*M] = Z[(k+M/2) + j*M];
        }
    }
    return rst;
}

bool check_end( const vector< complex<double> >& psi0, const vector< complex<double> >& psi1, 
                const vector<double>& xarr, const vector<double>& yarr)
{
    return false;
}

void exact() {
    // para
    double dx, dy, dkx, dky;
    vector<double> xarr = linspace(-L/2, L/2, M, dx);
    vector<double> yarr = linspace(-L/2, L/2, M, dy);
    dkx = 2 * M_PI / M / dx;
    dky = 2 * M_PI / M / dy;
    vector<double> kxarr = (arange(M) - M / 2) * dkx;
    vector<double> kyarr = (arange(M) - M / 2) * dky;
    double c0 = sqrt(1 - init_s);
    double c1 = sqrt(init_s);
    // for abc
    double x_absorb_accu0 = 0.0;
    double x_absorb_accu1 = 0.0;
    double y_absorb_accu0 = 0.0;
    double y_absorb_accu1 = 0.0;
    vector<double> reduce_x(M*M, 0.0), reduce_y(M*M, 0.0);
    if (enable_abc) {
        const double U0 = 0.1;
        const double alpha = 0.1;
        for (int k = 0; k < M; ++k) {
            for (int j = 0; j < M; ++j) {
                reduce_x[k+j*M] = (1.0 - U0 / pow(cosh(alpha * j), 2) * dt)
                                    * (1.0 - U0 / pow(cosh(alpha * (M-j)), 2) * dt);
                reduce_y[k+j*M] = (1.0 - U0 / pow(cosh(alpha * k), 2) * dt)
                                    * (1.0 - U0 / pow(cosh(alpha * (M-k)), 2) * dt);
            }
        }
    }
    // construct TU on k grid
    vector< complex<double> > TU(M * M);
    for (int k = 0; k < M; ++k) {
        double ky = kyarr[k];
        for (int j = 0; j < M; ++j) {
            double kx = kxarr[j];
            TU[k+j*M] = exp(-matrixop::IMAGIZ * dt * (kx*kx + ky*ky) / 2.0 / mass);
        }
    }
    TU = myfftshift(TU);
    // construct VU on x grid
    vector< complex<double> > V00(M*M), V01(M*M), V10(M*M), V11(M*M);
    vector< complex<double> > evts00(M*M), evts01(M*M), evts10(M*M), evts11(M*M);
    vector< complex<double> > H00(M*M), H01(M*M), H10(M*M), H11(M*M);
    for (int k = 0; k < M; ++k) {
        double y = yarr[k];
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            vector< complex<double> > H = cal_H(vector<double> {x, y});
            vector<double> eva;
            vector< complex<double> > evt, evamat;
            matrixop::eigh(H, eva, evt);
            // propagation matrix
            evamat.assign(4, 0.0);
            evamat[0+0*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[0]);
            evamat[1+1*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[1]);
            auto tmp = matrixop::matmatmatC(evt, evamat, evt, 2, 2);
            V00[k+j*M] = tmp[0+0*2];
            V01[k+j*M] = tmp[0+1*2];
            V10[k+j*M] = tmp[1+0*2];
            V11[k+j*M] = tmp[1+1*2];
            H00[k+j*M] = H[0+0*2];
            H01[k+j*M] = H[0+1*2];
            H10[k+j*M] = H[1+0*2];
            H11[k+j*M] = H[1+1*2];
            // eva & evt
            if (j == 0 and k == 0) {
                const complex<double> phase0 = evt[0+0*2] / abs(evt[0+0*2]);
                evt[0+0*2] *= conj(phase0);
                evt[1+0*2] *= conj(phase0);
                const complex<double> phase1 = evt[1+1*2] / abs(evt[1+1*2]);
                evt[0+1*2] *= conj(phase1);
                evt[1+1*2] *= conj(phase1);
            }
            else if (k == 0) {
                const complex<double> phase0 = conj(evts00[k+(j-1)*M]) * evt[0+0*2] + conj(evts10[k+(j-1)*M]) * evt[1+0*2];
                evt[0+0*2] *= conj(phase0);
                evt[1+0*2] *= conj(phase0);
                const complex<double> phase1 = conj(evts01[k+(j-1)*M]) * evt[0+1*2] + conj(evts11[k+(j-1)*M]) * evt[1+1*2];
                evt[0+1*2] *= conj(phase1);
                evt[1+1*2] *= conj(phase1);
            }
            else {
                const complex<double> phase0 = conj(evts00[(k-1)+j*M]) * evt[0+0*2] + conj(evts10[(k-1)+j*M]) * evt[1+0*2];
                evt[0+0*2] *= conj(phase0);
                evt[1+0*2] *= conj(phase0);
                const complex<double> phase1 = conj(evts01[(k-1)+j*M]) * evt[0+1*2] + conj(evts11[(k-1)+j*M]) * evt[1+1*2];
                evt[0+1*2] *= conj(phase1);
                evt[1+1*2] *= conj(phase1);
            }
            evts00[k+j*M] = evt[0+0*2];
            evts01[k+j*M] = evt[0+1*2];
            evts10[k+j*M] = evt[1+0*2];
            evts11[k+j*M] = evt[1+1*2];
        }
    }
    // initialize wf on adiabats
    vector< complex<double> > psiad0(M * M, 0.0), psiad1(M * M, 0.0);
    vector< complex<double> > psi0(M * M, 0.0), psi1(M * M, 0.0);
    for (int k = 0; k < M; ++k) {
        double y = yarr[k];
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            psiad0[k+j*M] = c0 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
            psiad1[k+j*M] = c1 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
        }
    }
    double nm = norm(psiad0 | psiad1);
    psiad0 /= nm;
    psiad1 /= nm;
    // convert to diab
    psi0 = evts00 * psiad0 + evts01 * psiad1;
    psi1 = evts10 * psiad0 + evts11 * psiad1;
    nm = norm(psi0 | psi1);
    psi0 /= nm;
    psi1 /= nm;
    // covinience vairables
    vector<int> dim{ M, M };
    // statistics
    double KE = 0.0, PE = 0.0;
    double n0d = 0.0, n1d = 0.0;
    double n0ad = 0.0, n1ad = 0.0;
    // propagate WF
    for (int istep = 0; istep < Nstep; ++istep) {
        // exp(-iVdt/2)
        vector< complex<double> > psi_k0 = V00 * psi0 + V01 * psi1;
        vector< complex<double> > psi_k1 = V10 * psi0 + V11 * psi1;
        // exp(-iTdt)
        psi_k0 = misc::fftn(psi_k0, dim);
        psi_k1 = misc::fftn(psi_k1, dim);
        psi_k0 *= TU;
        psi_k1 *= TU;
        // exp(-iVdt/2)
        psi_k0 = misc::ifftn(psi_k0, dim);
        psi_k1 = misc::ifftn(psi_k1, dim);
        psi0 = V00 * psi_k0 + V01 * psi_k1;
        psi1 = V10 * psi_k0 + V11 * psi_k1;

        // abosrbing boundary condition, notice energy conservation will break if abc is enabled
        if (enable_abc) {
            // reduce x boundary 
            double n0_before_reduce = sum(pow(abs(psi0), 2));
            double n1_before_reduce = sum(pow(abs(psi1), 2));

            psi0 *= reduce_x;
            psi1 *= reduce_x;

            double n0_after_reduce = sum(pow(abs(psi0), 2));
            double n1_after_reduce = sum(pow(abs(psi1), 2));

            x_absorb_accu0 += n0_before_reduce - n0_after_reduce;
            x_absorb_accu1 += n1_before_reduce - n1_after_reduce;

            // reduce y boundary
            n0_before_reduce = sum(pow(abs(psi0), 2));
            n1_before_reduce = sum(pow(abs(psi1), 2));

            psi0 *= reduce_y;
            psi1 *= reduce_y;

            n0_after_reduce = sum(pow(abs(psi0), 2));
            n1_after_reduce = sum(pow(abs(psi1), 2));

            y_absorb_accu0 += n0_before_reduce - n0_after_reduce;
            y_absorb_accu1 += n1_before_reduce - n1_after_reduce;
        }
        // analysis & output
        if (istep % output_step == 0) {
            // header
            if (istep == 0) {
                output_potential_param();
                ioer::info("# EXACT 2D: ");
                if (enable_adiab) {
                    ioer::info("# ADIAB ");
                }
                else {
                    ioer::info("# DIAB ");
                }
                ioer::info("# para: ", " L = ", L, " M = ", M, " mass = ", mass, 
                                       " xI = ", xI, " yI = ", yI, " sigmax = ", sigmax, " sigmay = ", sigmay, " kxI = ", kxI, " kyI = ", kyI, " init_s = ", init_s, " c0 = ", c0, " c1 = ", c1,
                                       " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step);
                ioer::info("# dx = ", dx, " dy = ", dy, " dkx = ", dkx, " dky = ", dky);
                if (enable_adiab) {
                    ioer::tabout('#', "t", "n0ad", "n1ad", "KE", "PE", "Etot");
                }
                else {
                    ioer::tabout('#', "t", "n0d", "n1d", "KE", "PE", "Etot");
                }
            }
            // get psi_k
            psi_k0 = myfftshift(misc::fftn(psi0, dim));
            psi_k1 = myfftshift(misc::fftn(psi1, dim));
            double nm = norm(psi_k0 | psi_k1);
            psi_k0 /= nm;
            psi_k1 /= nm;
            KE = 0.0;
            for (int k = 0; k < M; ++k) {
                for (int j = 0; j < M; ++j) {
                    KE += (pow(abs(psi_k0[k+j*M]), 2) + pow(abs(psi_k1[k+j*M]), 2)) * (kxarr[j]*kxarr[j] + kyarr[k]*kyarr[k]) * 0.5 / mass;
                }
            }
            PE = real(sum(conj(psi0) * H00 * psi0 + conj(psi0) * H01 * psi1 + conj(psi1) * H10 * psi0 + conj(psi1) * H11 * psi1));

            if (enable_adiab) {
                // analysis in adiab
                psiad0 = conj(evts00) * psi0 + conj(evts10) * psi1;
                psiad1 = conj(evts01) * psi0 + conj(evts11) * psi1;
                double nm = norm(psiad0 | psiad1);
                psiad0 /= nm;
                psiad1 /= nm;
                // get psiad fft
                auto psiad_k0 = myfftshift(misc::fftn(psiad0, dim));
                auto psiad_k1 = myfftshift(misc::fftn(psiad1, dim));
                nm = norm(psiad_k0 | psiad_k1);
                psiad_k0 /= nm;
                psiad_k1 /= nm;
                // analysis
                // XXX
            }
            else {
                n0d = sum(pow(abs(psi0), 2));
                n1d = sum(pow(abs(psi1), 2));
                ioer::tabout('#', istep * dt, n0d, n1d, KE, PE, KE+PE, 
                        x_absorb_accu0, x_absorb_accu1,
                        y_absorb_accu0, y_absorb_accu1,
                        x_absorb_accu0 + x_absorb_accu1 + y_absorb_accu0 + y_absorb_accu1
                        );
            }
            // check end
            if (check_end(psi0, psi1, xarr, yarr) == true) {
                ioer::info("# check_end returns true");
                break;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    exact();
    ioer::info("# ", timer::toc());
    return 0;
}
