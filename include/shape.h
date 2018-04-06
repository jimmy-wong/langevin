#include <string>
#include <gsl/gsl_vector.h>
#include <cmath>
#include <chrono>

using namespace std;
class shape{
public:
    shape() = default;
    shape(double l, double r, double z, double c, double s):
            _para_l(l), _para_r(r), _para_z(z), _para_c(c), _para_s(s){}
    double get_l(){ return _para_l;}
    double get_r(){ return _para_r;}
    double get_z(){ return _para_z;}
    double get_c(){ return _para_c;}
    double get_s(){ return _para_s;}
    double get_Rcn() { return _Rcn; }
    double get_excited_energy() { return _excited_energy;}
    double get_level_density() { return _para_a;}
    double get_density(){return _density;}
    double get_average_v(){return _average_v;}
    int* get_steps(){ return _steps;}
    void set_l(double l){ _para_l = l;}
    void set_r(double r){ _para_r = r;}
    void set_z(double z){ _para_z = z;}
    void set_c(double c){ _para_c = c;}
    void set_s(double s){ _para_s = s;}
    void set_excited_energy(double excited_energy){_excited_energy = excited_energy+_ground_state_energy;}
    void set_ground_energy(double ground_energy){_ground_state_energy = ground_energy;}
    void set_level_density() { _para_a = _Acn*(1.0+3.114*pow(_Acn,-1/3) + 5.626 * pow(_Acn,-2./3.))/14.61;}
    void set_steps(int steps[5]);
    void grid(double * starting_point, double * step_length, int* grid);
    void efficiency();
    double grid_energy(double* storation, int* step);
    double AH(gsl_vector* generalized_coordinates);
    double Rho(double zeta);
    double RhoDerivative(double zeta, char label);
    double RhoDDerivative(double zeta, char label_i, char label_j);
    double IntegrateRhoDerivative(double zeta, char label);
    double IntegrateRhoDDerivative(double zeta, char label_i, char label_j);
    double IntegrateRhoDerivativeLowHigh(double z_high, double z_low, char label_i);
    double IntegrateRhoDDerivativeLowHigh(double z_high, double z_low, char label_i, char label_j);
    double CenterOfMassDerivative(char side, char label);
private:
    // 1u = 931.5MeV, U236 mass = 236.045566201u=279876.444916MeV
    // \rho_m = U236 mass/(4\pi/3*_Rcn^3)=106.5378621MeV fm^-3
    // v_F = \sqrt[3]{9\pi Acn/4}/R_cn*\hbar/(U236 mass/236)
    double _Acn = 236.;
    double _Rcn = 1.27808*pow(_Acn,1./3.);
    double _para_l, _para_r, _para_z, _para_c, _para_s;
    double _a0 = 0., _a1 = 0., _a2 = 0., _a3 = 0., _a4 = 0.;
    double _density = 106.5378621; // unit MeV fm^-3
    double _average_v = 0.3161252227942924*4./3.; // unit 1
    double _excited_energy = 0, _ground_state_energy = 0.;
    double _para_a;// level density parameter: a
    int _steps[5];
};
double RhoShape(double zeta, void* params);
//  A_i
double A(shape shape, double z_high, double z_low, char label_i);
//  \frac{\partial A_i}{\partial z}==A_i'
double ADerivativeZ(shape shape, double z_high, double z_low, char label_i);
//  \frac{\partial A_i}{\partial q_l}
double ADerivativeQ(shape shape, double z_high, double z_low, char label_i, char label_l);
//  \frac{\partial A_i'}{\partial q_l}
double ADDerivativeZQ(shape shape, double z_high, double z_low, char label_i, char label_l);
