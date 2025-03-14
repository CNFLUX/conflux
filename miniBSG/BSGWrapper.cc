/// @file BSGWrapper.cc Wrapper for python-callable spectrum calculation

#include "SpectralFunctions.h"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <vector>
using std::vector;
#include <map>
using std::map;

using namespace bsg;
using namespace bsg::SpectralFunctions;

extern "C" {

double py_phase_space(double W, double W0, double numass = 0) {
    if(W <= 1) return 0;
    double nue = W0 - W - numass;
    if(nue <= numass) return 0;
    return W * nue * sqrt((W*W - 1) * (nue*nue - numass*numass));
}

double py_finite_size_L0(double W, double Z, double R) {
    return (
            1
            - ALPHA * Z * W * R
            + 13./60 * ALPHA*ALPHA*Z*Z
            - ALPHA*Z*R/W
            );
}

double py_recoil_gamow_teller(double W, double W0, double A) {
    double M = A * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV;
    double M2 = M*M;

    double Ar0 = -2. * W0 / 3. / M - W0 * W0 / 6. / M2 - 77. / 18. / M2;
    double Ar1 = -2. / 3. / M + 7. * W0 / 9. / M2;
    double Ar2 = 10. / 3. / M - 28. * W0 / 9. / M2;
    double Ar3 = 88. / 9. / M2;

    return (1
            +Ar0
            +Ar1/W
            +Ar2*W
            +Ar3*W*W);
}

double py_recoil_Coulomb_gamow_teller(double W, double W0, double Z, double A) {
    double M = A * NUCLEON_MASS_KEV / ELECTRON_MASS_KEV;
    double p = sqrt(W*W - 1);

    return (1
            -ALPHA*Z*M_PI/M/p
            *(1-1/3*(W0-W)/(3*W)));
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

vector<double>& get_aNeg(int Z) {

    static map<int, vector<double>> m_aNeg;
    auto it = m_aNeg.find(Z);
    if(it != m_aNeg.end()) return it->second;

    double bNeg[7][6];
    bNeg[0][0] = 0.115;
    bNeg[0][1] = -1.8123;
    bNeg[0][2] = 8.2498;
    bNeg[0][3] = -11.223;
    bNeg[0][4] = -14.854;
    bNeg[0][5] = 32.086;
    bNeg[1][0] = -0.00062;
    bNeg[1][1] = 0.007165;
    bNeg[1][2] = 0.01841;
    bNeg[1][3] = -0.53736;
    bNeg[1][4] = 1.2691;
    bNeg[1][5] = -1.5467;
    bNeg[2][0] = 0.02482;
    bNeg[2][1] = -0.5975;
    bNeg[2][2] = 4.84199;
    bNeg[2][3] = -15.3374;
    bNeg[2][4] = 23.9774;
    bNeg[2][5] = -12.6534;
    bNeg[3][0] = -0.14038;
    bNeg[3][1] = 3.64953;
    bNeg[3][2] = -38.8143;
    bNeg[3][3] = 172.1368;
    bNeg[3][4] = -346.708;
    bNeg[3][5] = 288.7873;
    bNeg[4][0] = 0.008152;
    bNeg[4][1] = -1.15664;
    bNeg[4][2] = 49.9663;
    bNeg[4][3] = -273.711;
    bNeg[4][4] = 657.6292;
    bNeg[4][5] = -603.7033;
    bNeg[5][0] = 1.2145;
    bNeg[5][1] = -23.9931;
    bNeg[5][2] = 149.9718;
    bNeg[5][3] = -471.2985;
    bNeg[5][4] = 662.1909;
    bNeg[5][5] = -305.6804;
    bNeg[6][0] = -1.5632;
    bNeg[6][1] = 33.4192;
    bNeg[6][2] = -255.1333;
    bNeg[6][3] = 938.5297;
    bNeg[6][4] = -1641.2845;
    bNeg[6][5] = 1095.358;

    vector<double> aNeg(7);
    for (int i = 0; i < 7; i++) {
        aNeg[i] = 0;
        for (int j = 0; j < 6; j++) {
            aNeg[i] += bNeg[i][j] * std::pow(ALPHA * Z, j + 1);
        }
    }

    return (m_aNeg[Z] = aNeg);
}

double BSG_FermiFunction(double W, double Z, double R) {
    // normalization matching python
    double g = sqrt(1 - ALPHA*ALPHA*Z*Z);
    return FermiFunction(W, Z, R, BETA_MINUS) * 2/(g+1);

    //return FermiFunction(W, Z, R, BETA_MINUS);
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
