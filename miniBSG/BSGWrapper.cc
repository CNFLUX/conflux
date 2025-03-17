/// @file BSGWrapper.cc Wrapper for python-callable spectrum calculation

#include "SpectralFunctions.h"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <array>
using std::array;

using namespace bsg;
using namespace bsg::SpectralFunctions;

extern "C" {

double py_phase_space(double W, double W0, double numass = 0) {
    if(W <= 1) return 0;
    double nue = W0 - W - numass;
    if(nue <= numass) return 0;
    return W * nue * sqrt((W*W - 1) * (nue*nue - numass*numass));
}

double py_shape_factor_gamow_teller(double W, double Z, double W0, double R, double A, double b, double c, double d, double Lambda) {
    double M = A * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV;
    double bt = Z/fabs(Z);

    double L = Lambda;

    double C0 = (-1./5*W0*W0*R*R
            +4./9*R*R*(1.-L/20)
            +1./3*W0/M/c*(-bt*b+d)
            +2./5*ALPHA*Z/(M*R*c)*(bt*2*b+d)
            +2./35*ALPHA*Z*W0*R*(1-L)
            -233./630*ALPHA*ALPHA*Z*Z);

    double C1 = (bt*4./3*b/M/c
            +4./9*W0*R*R*(1-L/10)
            -4./7*ALPHA*Z*R*(1-L/10));

    double Cm1 = (-1./3/M/c*(2*bt*b+d)
            -2./45*W0*R*R*(1-L)
            -ALPHA*Z*R/70);

    double C2 = -4./9*R*R*(1-L/10);

    return 1 + C0 + C1*W + Cm1/W + C2*W*W;
}

array<double,7>& get_aNeg(int Z) {

    static map<int, array<double,7>> m_aNeg;
    auto it = m_aNeg.find(Z);
    if(it != m_aNeg.end()) return it->second;

    constexpr double bNeg[7][6] = {
        {0.115, -1.8123, 8.2498, -11.223, -14.854, 32.086},
        {-0.00062, 0.007165, 0.01841, -0.53736, 1.2691, -1.5467},
        {0.02482, -0.5975, 4.84199, -15.3374, 23.9774, -12.6534},
        {-0.14038, 3.64953, -38.8143, 172.137, -346.708, 288.787},
        {0.008152, -1.15664, 49.9663, -273.711, 657.629, -603.703},
        {1.2145, -23.9931, 149.972, -471.298, 662.191, -305.68},
        {-1.5632, 33.4192, -255.133, 938.53, -1641.28, 1095.36}
    };

    array<double,7> aNeg = {};
    double c = 1;
    for(int j = 0; j < 6; j++) {
        c *= ALPHA*Z;
        for(int i = 0; i < 7; i++) aNeg[i] += bNeg[i][j] * c;
    }

    return (m_aNeg[Z] = aNeg);
}


array<double,7>& get_aPos(int Z) {

    static map<int, array<double,7>> m_aPos;
    auto it = m_aPos.find(Z);
    if(it != m_aPos.end()) return it->second;

    constexpr double bPos[7][6] = {
        {0.0701, -2.572, 27.5971, -128.658, 272.264, -214.925},
        {-0.002308, 0.066463, -0.6407, 2.63606, -5.6317, 4.0011},
        {0.07936, -2.09284, 18.4546, -80.9375, 160.838, -124.893},
        {-0.93832, 22.0251, -197.002, 807.188, -1566.61, 1156.33},
        {4.27618, -96.8241, 835.265, -3355.84, 6411.33, -4681.57},
        {-8.2135, 179.086, -1492.13, 5872.54, -11038.7, 7963.47},
        {5.4583, -115.892, 940.831, -3633.92, 6727.63, -4795.05}
    };

    array<double,7> aPos = {};
    double c = 1;
    for(int j = 0; j < 6; j++) {
        c *= ALPHA * Z;
        for(int i = 0; i < 7; i++) aPos[i] += bPos[i][j] * c;
    };

    return (m_aPos[Z] = aPos);
}


double BSG_FermiFunction(double W, double Z, double R) {
    // normalization changed to match python version
    double g = sqrt(1 - ALPHA*ALPHA*Z*Z);
    return FermiFunction(W, Z, R, BETA_MINUS) * 2/(g+1);
}

double BSG_L0Correction(double W, int Z, double R) {
    return L0Correction(W, Z, R, BETA_MINUS, nullptr, get_aNeg(Z).data());
}

double BSG_AtomicScreeningCorrection(double W, int Z) {
    if(Z == 1) return 1.;
    return AtomicScreeningCorrection(W, Z, BETA_MINUS);
}

double  BSG_RecoilCorrection(double W, double W0, int A) {
    return RecoilCorrection(W, W0, A, GAMOW_TELLER, 0);
}

double BSG_QCorrection(double W, double W0, int Z, int A) {
    return QCorrection(W, W0, Z, A, BETA_MINUS, GAMOW_TELLER, 0);
}

double BSG_RadiativeCorrection(double W, double W0, int Z, double R) {
    return RadiativeCorrection(W, W0, Z, R, BETA_MINUS, 1.27, 4.706);
}

double BSG_NeutrinoRadiativeCorrection(double Wv) {
    return NeutrinoRadiativeCorrection(Wv);
}



double gamma_k(int k, int Z) { return sqrt(k*k - ALPHA*ALPHA * Z*Z); }

/// Implementation of the generalized Fermi function F_{k-1} according to Behrens et al.
/*
    :param W: Total electron energy in units of its rest mass
    :param Z: Proton number of daughter
    :param R: nuclear radius in natural units
    :param k: absolute value of kappa
 */
double generalizedFermiFunction(double W, int Z, double R, int k) {
    if(W <= 1) return 0;
    double p = sqrt(W*W - 1.);
    double y = ALPHA*Z*W/p;
    double gk = gamma_k(k, Z);
    double kff = k * gsl_sf_doublefact(2*k-1);
    double prefactor = kff*kff * pow(4.,k) * pow(2*p*R, 2.*(gk-k)) * exp(M_PI*y);

    double lng = 0, arg = 0;
    gsl_sf_result magn;
    gsl_sf_result phase;
    gsl_sf_lngamma_complex_e(gk, y, &magn, &phase);

    double gammas = exp(magn.val) / gsl_sf_gamma(1 + 2*gk);
    return prefactor * gammas * gammas;
}

/// Coulomb function $\lambda_k$ as per Behrens et al.
double lambda_k(double W, int Z, double R, int k) {
    double gk = gamma_k(k, Z);
    double g1 = gamma_k(1, Z);
    return generalizedFermiFunction(W, Z, R, k) / generalizedFermiFunction(W, Z, R, 1) * (k+gk)/(k*(1+g1));
}

/* Unique forbidden shape factor

    :param L: int Spin change
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param Z: Proton number of the final nuclear state
*/
double shape_factor_unique_forbidden(double W, int L, double W0, int Z, double R) {
    if(!(W>1)) return 1;

    double pe = sqrt(W*W - 1.);
    double pnu = W0 - W;
    L = abs(L);
    double C = 0;
    for(int k = 1; k <= L; ++k)
        C += lambda_k(W, Z, R, k) * pow(pe, 2*(k-1)) * pow(pnu, 2*(L-k)) / gsl_sf_fact(2*k-1) / gsl_sf_fact(2*(L-k)+1);
    return gsl_sf_fact(2*L-1)*C;
}


double  BSG_beta_spectrum(
    double W,
    double W0,
    int A,
    int Z,
    double R,
    int L,
    bool isNu
) {
    if(W - 1 < 1e-5) return 0;
    double c = py_phase_space(W, W0);
    if(!c) return c;

    c *= BSG_FermiFunction(W, Z, R);

    // very small
    c *= BSG_RecoilCorrection(W, W0, A);
    c *= BSG_QCorrection(W, W0, Z, A);


    // don't have this in BSG... porting from thecobs .py function
    if(!L) c *= py_shape_factor_gamow_teller(W, Z, W0, R, A, A*4.7, 1.0, 0.0, 0.0);
    else c *= shape_factor_unique_forbidden(W, L, W0, Z, R);

    if(isNu) c *= BSG_NeutrinoRadiativeCorrection(W);
    else c *= BSG_RadiativeCorrection(W, W0, Z, R); // 7% effect, 1% +- 0.5% difference from thecobs

    c *= BSG_L0Correction(W, Z, R);

    c *= BSG_AtomicScreeningCorrection(W, Z);


    if(!(std::isfinite(c)) || c < 0) {
        printf("c = %g, W-1 = %g, W0 = %g, pf = %g, Z = %i R = %g, L = %i\n", c, W-1, W0, py_phase_space(W, W0), Z, R, L);
        throw std::runtime_error("invalid c");
    }

    return c;
}

}
