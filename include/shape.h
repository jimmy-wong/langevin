#include <string>
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
    int* get_steps(){ return _steps;}
    void set_l(double l){ _para_l = l;}
    void set_r(double r){ _para_r = r;}
    void set_z(double z){ _para_z = z;}
    void set_c(double c){ _para_c = c;}
    void set_s(double s){ _para_s = s;}
    void set_steps(int steps[5]);
    int* grid(double * starting_point, double * step_length);
    void efficiency();
    double grid_energy(double* storation, int* step, int* steps);
    double Rho(double zeta);
    double RhoDerivative(double zeta, char label);
    double A_Block(double z, char label);
    double CenterOfMassDerivative(char side, char label);
    double get_density(){return _density;}
    double get_average_v(){return _average_v;}
private:
    double _para_l, _para_r, _para_z, _para_c, _para_s;
    double _a0 = 0, _a1 = 0, _a2 = 0, _a3 = 0, _a4 = 0;
    double _density = 0;
    double _average_v = 0;
    int _steps[5];
};
double RhoShape(double zeta, void* params);
double A_derivative(shape shape, double z_high, double z_low, char label);
double A_dderivative(shape shape, double z_high, double z_low, char label);
double A_para(shape shape, double z_high, double z_low, char label);