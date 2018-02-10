#include "../include/shape.h"
#include <string>
#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace std;

void shape::set_steps(int *steps) {
    for(int i=0; i<5; i++){
        _steps[i] = steps[i];
    }
}
int* shape::grid(double * starting_point, double * step_length){
    double shape_para[5]={_para_l,_para_r,_para_z,_para_c,_para_s};
    int lower_limit[5];
    for(int i=0;i<5;i++){
        lower_limit[i] = floor((shape_para[i]-starting_point[i])/step_length[i]);
    }
    return lower_limit;
}
void shape::efficiency() {
    // shape_para unit: Rcn, not fm
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;
    // auxiliary variable
    double term0;
    double l2,l3,l4,l5,l6,l7,z2,z3,z4,z5,z6;

    l2=l*l;  l3=l2*l;  l4=l2*l2;  l5=l3*l2;  l6=l3*l3;  l7=l3*l4;
    z2=z*z;  z3=z2*z;  z4=z2*z2;  z5=z3*z2;  z6=z3*z3;

    term0 = l2-z2;
    _a0 = r*r/term0;
    _a1 = 2.*_a0*z/term0;
    _a2 = (c*r+_a0+2.*_a1*z)/term0;
    _a3 = (-7*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z -
            165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) +
            410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) -
            312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
            136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
            140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
            (pow(l,5)*pow((l2 - z2),3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245 *z6));
    _a4 = (-7*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z -
            60*pow(l,8)*z2 - 12*c*pow(l,11)*r*z2 + 90*pow(l,9)*pow(r,2)*z2 +
            140*l6*s*z3 + 270*l6*z4 + 50*c*pow(l,9)*r*z4 -
            330*l7*pow(r,2)*z4 + 210*l4*s*z5 - 300*l4*z6 -
            76*c*l7*r*z6 + 350*l5*pow(r,2)*z6 - 420*l2*s*pow(z,7) +
            105*l2*pow(z,8) + 35*c*l5*r*pow(z,8) + 175*s*pow(z,9)))/
            (l5*pow((l2 - z2),3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
}
double shape::grid_energy(double* storation, int* step, int* steps) {
    int index;
    index = step[0]*steps[1]*steps[2]*steps[3]*steps[4]+
            step[1]*steps[2]*steps[3]*steps[4]+
            step[2]*steps[3]*steps[4]+
            step[3]*steps[4]+
            step[4];
    return storation[index];
}
double shape::Rho(double x){
    double result;
    double l = _para_l;
    double z = _para_z;
    double s = _para_s;
    result = (pow(l,2) - pow(s + x,2))*(_a0 + _a1*(s + x - z) + _a2*pow(s + x - z,2) + _a3*pow(s + x - z,3) + _a4*pow(s + x - z,4));
    return result;
}

double shape::RhoDerivative(double x, const char label) {
    double result;
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;
    switch (label) {
        case 'x' :
            result = (pow(l, 2) - pow(s + x, 2)) *
                 (_a1 + 2 * _a2 * (s + x - z) + 3 * _a3 * pow(s + x - z, 2) + 4 * _a4 * pow(s + x - z, 3)) -
                 2 * (s + x) * (_a0 + _a1 * (s + x - z) + _a2 * pow(s + x - z, 2) + _a3 * pow(s + x - z, 3) +
                                _a4 * pow(s + x - z, 4));
        case 'l' :
            result = (pow(l,2) - pow(s + x,2))*((-8*l*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) -
                                                    (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2) -
                                                    (7*pow(s + x - z,4)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z -
                                                                           480*pow(l,7)*pow(z,2) - 132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) +
                                                                           840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) + 450*c*pow(l,8)*r*pow(z,4) -
                                                                           2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                           532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) +
                                                                           175*c*pow(l,4)*r*pow(z,8)))/
                                                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                  245*pow(z,6))) - (7*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z +
                                                                                                                                            792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) + 320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) +
                                                                                                                                            2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) + 612*c*pow(l,8)*r*pow(z,5) -
                                                                                                                                            2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                                                                                            952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) +
                                                                                                                                            350*c*pow(l,4)*r*pow(z,9)))/
                                                                                                                     (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                   245*pow(z,6))) + (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                                                                                                                                                       (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                                                                                                                                        270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                                                                                                                                                        210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                                                                                                                                                        350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                                                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                         245*pow(z,6),2)) + (42*pow(s + x - z,4)*
                                                                                                                                                                                                                                                               (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                                                                                                                                                                12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                                                                                                                                                                                                                270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                                                                                                                                                                                                                                210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                                                                                                                                                                                                                                350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                                                                                                                                                                35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                                                                                                                                                              (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                            245*pow(z,6))) + (35*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) -
                                                                                                                                                                                                                                                                                                                                                       105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) - 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) +
                                                                                                                                                                                                                                                                                                                                                       140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                       330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                                                                                                                                                                                                                                                                                                       76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) +
                                                                                                                                                                                                                                                                                                                                                       105*pow(l,2)*pow(z,8) + 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                                                                                                                                                                                                                               (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                             245*pow(z,6))) + (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                                                                                                                                                                                                                                                                                                                                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                                                                                                                                                                                                                                                                                                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) -
                                                                                                                                                                                                                                                                                                                                                                                                  360*pow(l,4)*pow(z,7) - 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) -
                                                                                                                                                                                                                                                                                                                                                                                                  315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                                                                                                                                                                                                                                                                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                   245*pow(z,6),2)) + (42*pow(s + x - z,3)*
                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          360*pow(l,4)*pow(z,7) - 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      245*pow(z,6))) + (35*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       245*pow(z,6))) + (pow(s + x - z,2)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) -
                                                    (2*l*pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                    pow(pow(l,2) - pow(z,2),2)) + 2*l*((2*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) +
                                                                                             pow(r,2)/(pow(l,2) - pow(z,2)) - (7*pow(s + x - z,4)*
                                                                                                                                     (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                                      12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                                                                                      270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                                                                                                      210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                                                                                                      350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                                      35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                  245*pow(z,6))) - (7*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z +
                                                                                                                                                                                                                            72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) +
                                                                                                                                                                                                                            410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) -
                                                                                                                                                                                                                            312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                                                                                                                            136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                                                                                                                                                            140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                                                                                                     (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                   245*pow(z,6))) + (pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                                                                                                                                                                                                                           pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
        case 'r' :
            result = (pow(l,2) - pow(s + x,2))*((4*r*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2)) -
                                                    (7*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                           50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                                           35*c*pow(l,5)*pow(z,8)))/
                                                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    - (7*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                             624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (pow(s + x - z,2)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/
                                                      (pow(l,2) - pow(z,2)));
        case 'z' :
            result = (pow(l,2) - pow(s + x,2))*((8*pow(r,2)*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) +
                                                    (2*pow(r,2)*(s + x - z))/pow(pow(l,2) - pow(z,2),2) -
                                                    (7*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z +
                                                                           420*pow(l,6)*s*pow(z,2) + 1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) -
                                                                           1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) - 1800*pow(l,4)*pow(z,5) -
                                                                           456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                           840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                       (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                        270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                        210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                        350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                         245*pow(z,6),2)) - (42*pow(s + x - z,4)*z*
                                                                                                                               (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                                12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                                                                                270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                                                                                                210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                                                                                                350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                                35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                              (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z -
                                                                              60*pow(l,8)*pow(z,2) - 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                              270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) +
                                                                              210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) +
                                                                              350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    - (7*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z +
                                                                             120*pow(l,8)*pow(z,2) - 24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) +
                                                                             340*c*pow(l,9)*r*pow(z,4) - 1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) -
                                                                             2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) + 3920*pow(l,5)*pow(r,2)*pow(z,6) -
                                                                             2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) + 1750*s*pow(z,9)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                       (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                        40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                        68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                        136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                        140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                         245*pow(z,6),2)) - (42*pow(s + x - z,3)*z*
                                                                                                                               (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                                                                40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                                                                68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                                136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                                                                140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                              (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z -
                                                                              165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) +
                                                                              240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) -
                                                                              360*pow(l,4)*pow(z,7) - 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) -
                                                                              315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    + (pow(s + x - z,2)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) +
                                                                           (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) +
                                                    (2*pow(s + x - z,2)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                    pow(pow(l,2) - pow(z,2),2) - (2*(s + x - z)*
                                                                                        (c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                                                       (pow(l,2) - pow(z,2)));
        case 'c' :
            result = (pow(l,2) - pow(s + x,2))*((r*pow(s + x - z,2))/(pow(l,2) - pow(z,2)) -
                                                    (7*pow(s + x - z,4)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                           35*pow(l,5)*r*pow(z,8)))/
                                                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
                                                    - (7*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) -
                                                                             136*pow(l,7)*r*pow(z,7) + 70*pow(l,5)*r*pow(z,9)))/
                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))
            );
        case 's' :
            result = (pow(l,2) - pow(s + x,2))*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) -
                    (7*pow(s + x - z,4)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) +
                            210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) - (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) -
                                                                                                                                             105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) - 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) +
                                                                                                                                             140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) -
                                                                                                                                             330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                                                                                             76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) +
                                                                                                                                             105*pow(l,2)*pow(z,8) + 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                                                     (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                   245*pow(z,6))) - (7*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) -
                                                                                                                                                                                                             90*pow(l,4)*pow(z,6) - 315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                                                                                                                                                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                    245*pow(z,6))) - (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z +
                                                                                                                                                                                                                                                                               72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) +
                                                                                                                                                                                                                                                                               410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) -
                                                                                                                                                                                                                                                                               312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                                                                                                                                                                               136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                                                                                                                                                                                                               140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                                                                                                                                                                       (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                     245*pow(z,6))) + (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                                                                                                                                                                                                                                                                        pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                     2*(s + x)*((2*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2)) -
                                (7*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z -
                                                       60*pow(l,8)*pow(z,2) - 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) +
                                                       140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) -
                                                       330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                       76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) +
                                                       105*pow(l,2)*pow(z,8) + 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                              245*pow(z,6))) - (7*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z +
                                                                                                                        72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) +
                                                                                                                        410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) -
                                                                                                                        312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                        136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                                                        140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                                 (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                               245*pow(z,6))) + (pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                                                                                                                       pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    return result;
}
double shape::RhoDDerivative(double x, char label_i, char label_j) {
    double result;
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;
    if (label_i=='x'&&label_j=='l'){
        result = (pow(l,2) - pow(s + x,2))*((-8*l*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),3) -
                                                (28*pow(s + x - z,3)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                        132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                        450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                        532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                        320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                        612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                        952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (168*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (140*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (126*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (105*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                                                (pow(l,2) - pow(z,2)) - (4*l*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                              pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) -
                 2*(s + x)*((-8*l*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2) -
                            (7*pow(s + x - z,4)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                   132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                   450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                   532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                            /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                   320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                   612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                   952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                            /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                             (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                            (42*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (35*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                             (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                            (42*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (35*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                            (pow(l,2) - pow(z,2)) - (2*l*pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                 pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) +
                 2*l*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) - (28*pow(s + x - z,3)*
                                                                           (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                            12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                            50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                            76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                            35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                          (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='x' && label_j=='r'){
        result = (pow(l,2) - pow(s + x,2))*((4*r*z)/pow(pow(l,2) - pow(z,2),2) -
                                                (28*pow(s + x - z,3)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                        50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                                        35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                        624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 2*(s + x)*((4*r*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2)) -
                            (7*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                   50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                   35*c*pow(l,5)*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                   624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='x'&& label_j=='z'){
        result = (pow(l,2) - pow(s + x,2))*((8*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (2*pow(r,2))/pow(pow(l,2) - pow(z,2),2) -
                                                (28*pow(s + x - z,3)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                        1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                        1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                        840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (168*pow(s + x - z,3)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (84*pow(s + x - z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                        24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                        1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                        3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                        1750*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (126*pow(s + x - z,2)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (42*(s + x - z)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                                                (pow(l,2) - pow(z,2)) + (4*(s + x - z)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                              pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                                                (2*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 2*(s + x)*((8*pow(r,2)*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (2*pow(r,2)*(s + x - z))/pow(pow(l,2) - pow(z,2),2) -
                            (7*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                   1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                   1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                   840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                             (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                            (42*pow(s + x - z,4)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                      12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                      50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                      76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                      35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                   24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                   1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                   3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                   1750*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                             (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                            (42*pow(s + x - z,3)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                      40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                      68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                      136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                      70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                            (pow(l,2) - pow(z,2)) + (2*pow(s + x - z,2)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                 pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                            (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='x'&&label_j=='c'){
        result = (pow(l,2) - pow(s + x,2))*((2*r*(s + x - z))/(pow(l,2) - pow(z,2)) -
                                                (28*pow(s + x - z,3)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                        35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                        70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))) -
                 2*(s + x)*((r*pow(s + x - z,2))/(pow(l,2) - pow(z,2)) -
                            (7*pow(s + x - z,4)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                   35*pow(l,5)*r*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                   70*pow(l,5)*r*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))));
    }
    else if (label_i=='x'&& label_j=='s'){
        result = (pow(l,2) - pow(s + x,2))*((-28*pow(s + x - z,3)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) -
                                                                         420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (84*pow(s + x - z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                        315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*(s + x - z)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 2*(s + x)*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) -
                            (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)))
                 - 2*(s + x)*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) -
                              (7*pow(s + x - z,4)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                              (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                              (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                      12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                      50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                      76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                      35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                              (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                              (7*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                     315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                              (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                              (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                      40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                      68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                      136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                      70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                              (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                              (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)))
                 - 2*((2*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2)) -
                      (7*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                             12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                             50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                             76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                             35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (7*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                             40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                             68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                             136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                             70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                      (pow(l,2) - pow(z,2)));
    }
    else if (label_i=='l'&&label_j=='l'){
        result = (pow(l,2) - pow(s + x,2))*((48*pow(l,2)*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),4) +
                                                (8*pow(l,2)*pow(r,2))/pow(pow(l,2) - pow(z,2),3) - (8*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) -
                                                (2*pow(r,2))/pow(pow(l,2) - pow(z,2),2) - (7*pow(s + x - z,4)*
                                                                                                   (-1350*pow(l,8) + 468*c*pow(l,11)*r + 1980*pow(l,9)*pow(r,2) - 5880*pow(l,6)*s*z - 3360*pow(l,6)*pow(z,2) -
                                                                                                    1320*c*pow(l,9)*r*pow(z,2) + 6480*pow(l,7)*pow(r,2)*pow(z,2) + 4200*pow(l,4)*s*pow(z,3) + 8100*pow(l,4)*pow(z,4) +
                                                                                                    3600*c*pow(l,7)*r*pow(z,4) - 13860*pow(l,5)*pow(r,2)*pow(z,4) + 2520*pow(l,2)*s*pow(z,5) - 3600*pow(l,2)*pow(z,6) -
                                                                                                    3192*c*pow(l,5)*r*pow(z,6) + 7000*pow(l,3)*pow(r,2)*pow(z,6) - 840*s*pow(z,7) + 210*pow(z,8) + 700*c*pow(l,3)*r*pow(z,8)))/
                                                                                                  (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (14*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                  132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                  450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                  532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (84*pow(s + x - z,4)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                        132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                        450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                        532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (70*pow(s + x - z,4)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                        132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                        450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                        532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(-1350*pow(l,8)*s - 5400*pow(l,8)*z + 936*c*pow(l,11)*r*z + 7920*pow(l,9)*pow(r,2)*z -
                                                                       9240*pow(l,6)*s*pow(z,2) + 2240*pow(l,6)*pow(z,3) - 880*c*pow(l,9)*r*pow(z,3) + 12300*pow(l,4)*s*pow(z,4) +
                                                                       7200*pow(l,4)*pow(z,5) + 4896*c*pow(l,7)*r*pow(z,5) - 13104*pow(l,5)*pow(r,2)*pow(z,5) - 1080*pow(l,2)*s*pow(z,6) -
                                                                       4320*pow(l,2)*pow(z,7) - 5712*c*pow(l,5)*r*pow(z,7) + 11200*pow(l,3)*pow(r,2)*pow(z,7) - 630*s*pow(z,8) + 280*pow(z,9) +
                                                                       1400*c*pow(l,3)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (14*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                  320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                  612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                  952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (84*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                        320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                        612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                        952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (70*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                        320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                        612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                        952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,4)*pow(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4),2)*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,4)*(270*pow(l,4) - 756*pow(l,2)*pow(z,2) - 42*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (84*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (70*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (336*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,3)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (378*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (210*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,7)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,3)*pow(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4),2)*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,3)*(270*pow(l,4) - 756*pow(l,2)*pow(z,2) - 42*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (84*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (70*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (336*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,3)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (378*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (210*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,7)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((96*pow(l,2)*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),4) +
                                                                     (8*pow(l,2)*pow(r,2))/pow(pow(l,2) - pow(z,2),3) - (16*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) -
                                                                     (2*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) -
                                                (4*l*pow(s + x - z,2)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                                                pow(pow(l,2) - pow(z,2),2) + (8*pow(l,2)*pow(s + x - z,2)*
                                                                                    (c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),3) -
                                                (2*pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                pow(pow(l,2) - pow(z,2),2)) + 4*l*((-8*l*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) -
                                                                                         (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2) - (7*pow(s + x - z,4)*
                                                                                                                                              (-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                                                                                               132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                                                                                               450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                                                                                               532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                                                                                                             /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                                                         (7*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                                                                320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                                                                612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                                                                952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                                                         /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                                                         (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                                                          (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                         (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                                                         (42*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                                                                 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                                                                 76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                         (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                                                         (35*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                                                                 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                                                                 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                                                                 76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                                                                 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                                                         (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                                                         (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                                                          (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                         (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                                                         (42*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                         (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                                                         (35*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                                                         (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                                                         (pow(s + x - z,2)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                                                                                         (pow(l,2) - pow(z,2)) - (2*l*pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                                                                              pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) +
                 2*((2*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2)) -
                    (7*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                    (7*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                    (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                    (pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                    (pow(l,2) - pow(z,2)));
    }
    else if (label_i=='l'&&label_j=='r'){
        result = (pow(l,2) - pow(s + x,2))*((-16*l*r*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) - (4*l*r)/pow(pow(l,2) - pow(z,2),2) -
                                                (7*pow(s + x - z,4)*(39*c*pow(l,12) + 396*pow(l,10)*r - 132*c*pow(l,10)*pow(z,2) + 1620*pow(l,8)*r*pow(z,2) +
                                                                       450*c*pow(l,8)*pow(z,4) - 4620*pow(l,6)*r*pow(z,4) - 532*c*pow(l,6)*pow(z,6) + 3500*pow(l,4)*r*pow(z,6) +
                                                                       175*c*pow(l,4)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) + 50*c*pow(l,9)*pow(z,4) -
                                                  660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) + 35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                        50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                                        35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                        50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                                        35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(78*c*pow(l,12)*z + 1584*pow(l,10)*r*z - 88*c*pow(l,10)*pow(z,3) + 612*c*pow(l,8)*pow(z,5) -
                                                                       4368*pow(l,6)*r*pow(z,5) - 952*c*pow(l,6)*pow(z,7) + 5600*pow(l,4)*r*pow(z,7) + 350*c*pow(l,4)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) - 624*pow(l,7)*r*pow(z,5) -
                                                  136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                        624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                        624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((-32*l*r*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (4*l*r)/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) -
                                                (2*l*pow(s + x - z,2)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2))
                 + 2*l*((4*r*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2)) -
                        (7*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                               50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                               35*c*pow(l,5)*pow(z,8)))/
                        (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                        (7*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                               624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                        (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                        (pow(s + x - z,2)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if(label_i=='l'&&label_j=='z'){
        result = (pow(l,2) - pow(s + x,2))*((-48*l*pow(r,2)*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),4) -
                                                (8*l*pow(r,2)*(s + x - z))/pow(pow(l,2) - pow(z,2),3) -
                                                (7*pow(s + x - z,4)*(-840*pow(l,7)*s - 960*pow(l,7)*z - 264*c*pow(l,10)*r*z + 1620*pow(l,8)*pow(r,2)*z + 2520*pow(l,5)*s*pow(z,2) +
                                                                       6480*pow(l,5)*pow(z,3) + 1800*c*pow(l,8)*r*pow(z,3) - 9240*pow(l,6)*pow(r,2)*pow(z,3) + 4200*pow(l,3)*s*pow(z,4) -
                                                                       7200*pow(l,3)*pow(z,5) - 3192*c*pow(l,6)*r*pow(z,5) + 10500*pow(l,4)*pow(r,2)*pow(z,5) - 5880*l*s*pow(z,6) +
                                                                       1680*l*pow(z,7) + 1400*c*pow(l,4)*r*pow(z,7)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                  132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                  450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                  532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,4)*z*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                          132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                          450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                          532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                        132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                        450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                        532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(-600*pow(l,9) + 78*c*pow(l,12)*r + 792*pow(l,10)*pow(r,2) - 2640*pow(l,7)*s*z + 960*pow(l,7)*pow(z,2) -
                                                                       264*c*pow(l,10)*r*pow(z,2) + 9840*pow(l,5)*s*pow(z,3) + 7200*pow(l,5)*pow(z,4) + 3060*c*pow(l,8)*r*pow(z,4) -
                                                                       10920*pow(l,6)*pow(r,2)*pow(z,4) - 2160*pow(l,3)*s*pow(z,5) - 10080*pow(l,3)*pow(z,6) - 6664*c*pow(l,6)*r*pow(z,6) +
                                                                       19600*pow(l,4)*pow(r,2)*pow(z,6) - 5040*l*s*pow(z,7) + 2520*l*pow(z,8) + 3150*c*pow(l,4)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                  1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                  1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                  840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                        1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                        1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                        840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                        1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                        1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                        840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                  320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                  612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                  952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,3)*z*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z -
                                                                          1320*pow(l,7)*s*pow(z,2) + 320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) +
                                                                          1440*pow(l,5)*pow(z,5) + 612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) -
                                                                          1440*pow(l,3)*pow(z,7) - 952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) +
                                                                          350*c*pow(l,4)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                        320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                        612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                        952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,4)*(-504*pow(l,3)*z - 168*l*pow(z,3))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,4)*z*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (28*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (35*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (336*pow(s + x - z,4)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (168*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (210*pow(s + x - z,4)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (140*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                  24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                  1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                  3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                  1750*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                        24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                        1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                        3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                        1750*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                        24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                        1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                        3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                        1750*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,3)*(-504*pow(l,3)*z - 168*l*pow(z,3))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,3)*z*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (21*pow(s + x - z,2)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (35*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (336*pow(s + x - z,3)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (126*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (210*pow(s + x - z,3)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (105*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((-96*l*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),4) - (40*l*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),3)))/
                                                (pow(l,2) - pow(z,2)) + (2*pow(s + x - z,2)*z*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) -
                                                                                                     (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/pow(pow(l,2) - pow(z,2),2) -
                                                (2*(s + x - z)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                                                (pow(l,2) - pow(z,2)) - (2*l*pow(s + x - z,2)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) +
                                                                                                     (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/pow(pow(l,2) - pow(z,2),2) -
                                                (8*l*pow(s + x - z,2)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                pow(pow(l,2) - pow(z,2),3) + (4*l*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                                     pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) +
                 2*l*((8*pow(r,2)*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (2*pow(r,2)*(s + x - z))/pow(pow(l,2) - pow(z,2),2) -
                      (7*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                             1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                             1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                             840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                       (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                      (42*pow(s + x - z,4)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (7*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                             24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                             1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                             3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                             1750*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                       (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                        40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                        68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                        136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                        70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                      (42*pow(s + x - z,3)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (pow(s + x - z,2)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                      (pow(l,2) - pow(z,2)) + (2*pow(s + x - z,2)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                           pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                      (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='l'&&label_j=='c'){
        result = (pow(l,2) - pow(s + x,2))*((-2*l*r*pow(s + x - z,2))/pow(pow(l,2) - pow(z,2),2) -
                                                (7*pow(s + x - z,4)*(39*pow(l,12)*r - 132*pow(l,10)*r*pow(z,2) + 450*pow(l,8)*r*pow(z,4) - 532*pow(l,6)*r*pow(z,6) +
                                                                       175*pow(l,4)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) + 35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,4)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                        35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,4)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                        35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(78*pow(l,12)*r*z - 88*pow(l,10)*r*pow(z,3) + 612*pow(l,8)*r*pow(z,5) - 952*pow(l,6)*r*pow(z,7) +
                                                                       350*pow(l,4)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) + 70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                        70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                        70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))) +
                 2*l*((r*pow(s + x - z,2))/(pow(l,2) - pow(z,2)) - (7*pow(s + x - z,4)*
                                                                          (3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) + 35*pow(l,5)*r*pow(z,8)))/
                                                                         (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (7*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                             70*pow(l,5)*r*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))));
    }
    else if (label_i=='l'&&label_j=='s'){
        result = (pow(l,2) - pow(s + x,2))*((-8*l*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),3) -
                                                (7*pow(s + x - z,4)*(-840*pow(l,7)*z + 840*pow(l,5)*pow(z,3) + 840*pow(l,3)*pow(z,5) - 840*l*pow(z,7)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(-150*pow(l,9) - 1320*pow(l,7)*pow(z,2) + 2460*pow(l,5)*pow(z,4) - 360*pow(l,3)*pow(z,6) - 630*l*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (28*pow(s + x - z,3)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                                        132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                                        450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                                        532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,4)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,4)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                                        320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                                        612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                                        952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                                                /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (168*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (140*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                         12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                         50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                         76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                         35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) - 315*pow(l,2)*pow(z,8) +
                                                  175*pow(z,10)))/(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*
                                                                     pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (42*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                        315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (35*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                        315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (126*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (105*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                         40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                         68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                         136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                         70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                                                (pow(l,2) - pow(z,2)) - (4*l*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                              pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) -
                 2*(s + x)*((-8*l*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2) -
                            (7*pow(s + x - z,4)*(-150*pow(l,9) + 39*c*pow(l,12)*r + 198*pow(l,10)*pow(r,2) - 840*pow(l,7)*s*z - 480*pow(l,7)*pow(z,2) -
                                                   132*c*pow(l,10)*r*pow(z,2) + 810*pow(l,8)*pow(r,2)*pow(z,2) + 840*pow(l,5)*s*pow(z,3) + 1620*pow(l,5)*pow(z,4) +
                                                   450*c*pow(l,8)*r*pow(z,4) - 2310*pow(l,6)*pow(r,2)*pow(z,4) + 840*pow(l,3)*s*pow(z,5) - 1200*pow(l,3)*pow(z,6) -
                                                   532*c*pow(l,6)*r*pow(z,6) + 1750*pow(l,4)*pow(r,2)*pow(z,6) - 840*l*s*pow(z,7) + 210*l*pow(z,8) + 175*c*pow(l,4)*r*pow(z,8)))
                            /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(-150*pow(l,9)*s - 600*pow(l,9)*z + 78*c*pow(l,12)*r*z + 792*pow(l,10)*pow(r,2)*z - 1320*pow(l,7)*s*pow(z,2) +
                                                   320*pow(l,7)*pow(z,3) - 88*c*pow(l,10)*r*pow(z,3) + 2460*pow(l,5)*s*pow(z,4) + 1440*pow(l,5)*pow(z,5) +
                                                   612*c*pow(l,8)*r*pow(z,5) - 2184*pow(l,6)*pow(r,2)*pow(z,5) - 360*pow(l,3)*s*pow(z,6) - 1440*pow(l,3)*pow(z,7) -
                                                   952*c*pow(l,6)*r*pow(z,7) + 2800*pow(l,4)*pow(r,2)*pow(z,7) - 630*l*s*pow(z,8) + 280*l*pow(z,9) + 350*c*pow(l,4)*r*pow(z,9)))
                            /(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,4)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                             (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                            (42*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (35*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,3)*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                             (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                            (42*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (35*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*((-16*l*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) - (2*l*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/
                            (pow(l,2) - pow(z,2)) - (2*l*pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                 pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2)) +
                 2*l*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) - (7*pow(s + x - z,4)*
                                                                           (-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                                          (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (7*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                             315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='r'&&label_j=='r'){
        result = (pow(l,2) - pow(s + x,2))*((4*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + 2/(pow(l,2) - pow(z,2)) -
                                                (7*pow(s + x - z,4)*(36*pow(l,11) + 180*pow(l,9)*pow(z,2) - 660*pow(l,7)*pow(z,4) + 700*pow(l,5)*pow(z,6)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(144*pow(l,11)*z - 624*pow(l,7)*pow(z,5) + 1120*pow(l,5)*pow(z,7)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((8*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + 2/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if(label_i=='r'&&label_j=='z'){
        result = (pow(l,2) - pow(s + x,2))*((16*r*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (4*r*(s + x - z))/pow(pow(l,2) - pow(z,2),2) -
                                                (7*pow(s + x - z,4)*(-24*c*pow(l,11)*z + 360*pow(l,9)*r*z + 200*c*pow(l,9)*pow(z,3) - 2640*pow(l,7)*r*pow(z,3) -
                                                                       456*c*pow(l,7)*pow(z,5) + 4200*pow(l,5)*r*pow(z,5) + 280*c*pow(l,5)*pow(z,7)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) + 50*c*pow(l,9)*pow(z,4) -
                                                  660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) + 35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,4)*z*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                          50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) + 35*c*pow(l,5)*pow(z,8)
                                                ))/(pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                        50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) + 35*c*pow(l,5)*pow(z,8)
                                                ))/(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(6*c*pow(l,13) + 144*pow(l,11)*r - 24*c*pow(l,11)*pow(z,2) + 340*c*pow(l,9)*pow(z,4) -
                                                                       3120*pow(l,7)*r*pow(z,4) - 952*c*pow(l,7)*pow(z,6) + 7840*pow(l,5)*r*pow(z,6) + 630*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) - 624*pow(l,7)*r*pow(z,5) -
                                                  136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,3)*z*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                          624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                        624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((32*r*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (20*r*z)/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) +
                                                (2*pow(s + x - z,2)*z*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/
                                                pow(pow(l,2) - pow(z,2),2) - (2*(s + x - z)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/
                                                                                   (pow(l,2) - pow(z,2)));
    }
    else if (label_i=='r'&&label_j=='c'){
        result = (pow(l,2) - pow(s + x,2))*(pow(s + x - z,2)/(pow(l,2) - pow(z,2)) -
                                                (7*pow(s + x - z,4)*(3*pow(l,13) - 12*pow(l,11)*pow(z,2) + 50*pow(l,9)*pow(z,4) - 76*pow(l,7)*pow(z,6) + 35*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(6*pow(l,13)*z - 8*pow(l,11)*pow(z,3) + 68*pow(l,9)*pow(z,5) - 136*pow(l,7)*pow(z,7) +
                                                                       70*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))));
    }
    else if (label_i=='r'&&label_j=='s'){
        result = (pow(l,2) - pow(s + x,2))*((4*r*z)/pow(pow(l,2) - pow(z,2),2) -
                                                (28*pow(s + x - z,3)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                                        50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                                        35*c*pow(l,5)*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                                        624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 2*(s + x)*((4*r*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2)) -
                            (7*pow(s + x - z,4)*(3*c*pow(l,13) + 36*pow(l,11)*r - 12*c*pow(l,11)*pow(z,2) + 180*pow(l,9)*r*pow(z,2) +
                                                   50*c*pow(l,9)*pow(z,4) - 660*pow(l,7)*r*pow(z,4) - 76*c*pow(l,7)*pow(z,6) + 700*pow(l,5)*r*pow(z,6) +
                                                   35*c*pow(l,5)*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(6*c*pow(l,13)*z + 144*pow(l,11)*r*z - 8*c*pow(l,11)*pow(z,3) + 68*c*pow(l,9)*pow(z,5) -
                                                   624*pow(l,7)*r*pow(z,5) - 136*c*pow(l,7)*pow(z,7) + 1120*pow(l,5)*r*pow(z,7) + 70*c*pow(l,5)*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*(c + (8*r*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + (2*r)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='z'&&label_j=='z'){
        result = (pow(l,2) - pow(s + x,2))*((48*pow(r,2)*(s + x - z)*pow(z,3))/pow(pow(l,2) - pow(z,2),4) +
                                                (24*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),3) - (8*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) -
                                                (2*pow(r,2))/pow(pow(l,2) - pow(z,2),2) - (7*pow(s + x - z,4)*
                                                                                                   (-120*pow(l,8) - 24*c*pow(l,11)*r + 180*pow(l,9)*pow(r,2) + 840*pow(l,6)*s*z + 3240*pow(l,6)*pow(z,2) +
                                                                                                    600*c*pow(l,9)*r*pow(z,2) - 3960*pow(l,7)*pow(r,2)*pow(z,2) + 4200*pow(l,4)*s*pow(z,3) - 9000*pow(l,4)*pow(z,4) -
                                                                                                    2280*c*pow(l,7)*r*pow(z,4) + 10500*pow(l,5)*pow(r,2)*pow(z,4) - 17640*pow(l,2)*s*pow(z,5) + 5880*pow(l,2)*pow(z,6) +
                                                                                                    1960*c*pow(l,5)*r*pow(z,6) + 12600*s*pow(z,7)))/
                                                                                                  (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (14*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                  1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                  1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                  840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (84*pow(s + x - z,4)*z*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                          1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                          1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                          840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (56*pow(s + x - z,3)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                        1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                        1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                        840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(-330*pow(l,8)*s + 240*pow(l,8)*z - 48*c*pow(l,11)*r*z + 4920*pow(l,6)*s*pow(z,2) + 4800*pow(l,6)*pow(z,3) +
                                                                       1360*c*pow(l,9)*r*pow(z,3) - 6240*pow(l,7)*pow(r,2)*pow(z,3) - 2700*pow(l,4)*s*pow(z,4) - 15120*pow(l,4)*pow(z,5) -
                                                                       5712*c*pow(l,7)*r*pow(z,5) + 23520*pow(l,5)*pow(r,2)*pow(z,5) - 17640*pow(l,2)*s*pow(z,6) + 10080*pow(l,2)*pow(z,7) +
                                                                       5040*c*pow(l,5)*r*pow(z,7) + 15750*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,4)*pow(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5),2)*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,4)*(-126*pow(l,4) - 252*pow(l,2)*pow(z,2) - 7350*pow(z,4))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (84*pow(s + x - z,4)*z*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (56*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (336*pow(s + x - z,4)*pow(z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z -
                                                                                    60*pow(l,8)*pow(z,2) - 12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) +
                                                                                    270*pow(l,6)*pow(z,4) + 50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) -
                                                                                    300*pow(l,4)*pow(z,6) - 76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) +
                                                                                    105*pow(l,2)*pow(z,8) + 35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (336*pow(s + x - z,3)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (84*pow(s + x - z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (14*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                  24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                  1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                  3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                  1750*s*pow(z,9)))/(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*
                                                                       pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (84*pow(s + x - z,3)*z*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                          24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                          1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                          3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                          1750*s*pow(z,9)))/(pow(l,5)*pow(pow(l,2) - pow(z,2),4)*
                                                                                               (9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (42*pow(s + x - z,2)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                        24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                        1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                        3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                        1750*s*pow(z,9)))/(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*
                                                                                             (9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (14*pow(s + x - z,3)*pow(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5),2)*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),3)) +
                                                (7*pow(s + x - z,3)*(-126*pow(l,4) - 252*pow(l,2)*pow(z,2) - 7350*pow(z,4))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) +
                                                (84*pow(s + x - z,3)*z*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,2)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (336*pow(s + x - z,3)*pow(z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z -
                                                                                    165*pow(l,8)*s*pow(z,2) + 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) +
                                                                                    240*pow(l,6)*pow(z,5) + 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) -
                                                                                    360*pow(l,4)*pow(z,7) - 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) +
                                                                                    140*pow(l,2)*pow(z,9) + 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),5)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                        40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                        68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                        136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                        70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (252*pow(s + x - z,2)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*(s + x - z)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (pow(s + x - z,2)*((96*pow(r,2)*pow(z,4))/pow(pow(l,2) - pow(z,2),4) + (88*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) +
                                                                     (10*pow(r,2))/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) +
                                                (4*pow(s + x - z,2)*z*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                                                pow(pow(l,2) - pow(z,2),2) - (4*(s + x - z)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) +
                                                                                                   (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/(pow(l,2) - pow(z,2)) +
                                                (8*pow(s + x - z,2)*pow(z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                pow(pow(l,2) - pow(z,2),3) + (2*pow(s + x - z,2)*
                                                                                    (c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                                                (8*(s + x - z)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                pow(pow(l,2) - pow(z,2),2) + (2*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                                                                                   (pow(l,2) - pow(z,2)));
    }
    else if (label_i=='z'&&label_j=='c'){
        result = (pow(l,2) - pow(s + x,2))*((2*r*pow(s + x - z,2)*z)/pow(pow(l,2) - pow(z,2),2) - (2*r*(s + x - z))/(pow(l,2) - pow(z,2)) -
                                                (7*pow(s + x - z,4)*(-24*pow(l,11)*r*z + 200*pow(l,9)*r*pow(z,3) - 456*pow(l,7)*r*pow(z,5) + 280*pow(l,5)*r*pow(z,7)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) + 35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,4)*z*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                          35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                        35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(6*pow(l,13)*r - 24*pow(l,11)*r*pow(z,2) + 340*pow(l,9)*r*pow(z,4) - 952*pow(l,7)*r*pow(z,6) +
                                                                       630*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) + 70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,3)*z*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                          70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                        70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))));
    }
    else if (label_i=='z'&&label_j=='s'){
        result = (pow(l,2) - pow(s + x,2))*((8*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (2*pow(r,2))/pow(pow(l,2) - pow(z,2),2) -
                                                (7*pow(s + x - z,4)*(-105*pow(l,8) + 420*pow(l,6)*pow(z,2) + 1050*pow(l,4)*pow(z,4) - 2940*pow(l,2)*pow(z,6) + 1575*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (28*pow(s + x - z,3)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                                        1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                                        1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                                        840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,4)*z*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (7*pow(s + x - z,3)*(-330*pow(l,8)*z + 1640*pow(l,6)*pow(z,3) - 540*pow(l,4)*pow(z,5) - 2520*pow(l,2)*pow(z,7) + 1750*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (28*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                  12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                  50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                  76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                  35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (168*pow(s + x - z,3)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                           12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                           50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                           76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                           35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (84*pow(s + x - z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                                        24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                                        1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                                        3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                                        1750*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) - 315*pow(l,2)*pow(z,8) +
                                                  175*pow(z,10)))/(pow(l,5)*pow(pow(l,2) - pow(z,2),3)*
                                                                     pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (42*pow(s + x - z,3)*z*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                          315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                        315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (21*pow(s + x - z,2)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                 (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                  40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                  68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                  136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                  70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                                                (126*pow(s + x - z,2)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                           40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                           68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                           136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                           70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (42*(s + x - z)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(s + x - z)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                                                (pow(l,2) - pow(z,2)) + (4*(s + x - z)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                              pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                                                (2*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 2*(s + x)*((8*pow(r,2)*(s + x - z)*pow(z,2))/pow(pow(l,2) - pow(z,2),3) + (2*pow(r,2)*(s + x - z))/pow(pow(l,2) - pow(z,2),2) -
                            (7*pow(s + x - z,4)*(-105*pow(l,8)*s - 120*pow(l,8)*z - 24*c*pow(l,11)*r*z + 180*pow(l,9)*pow(r,2)*z + 420*pow(l,6)*s*pow(z,2) +
                                                   1080*pow(l,6)*pow(z,3) + 200*c*pow(l,9)*r*pow(z,3) - 1320*pow(l,7)*pow(r,2)*pow(z,3) + 1050*pow(l,4)*s*pow(z,4) -
                                                   1800*pow(l,4)*pow(z,5) - 456*c*pow(l,7)*r*pow(z,5) + 2100*pow(l,5)*pow(r,2)*pow(z,5) - 2940*pow(l,2)*s*pow(z,6) +
                                                   840*pow(l,2)*pow(z,7) + 280*c*pow(l,5)*r*pow(z,7) + 1575*s*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,4)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                             (-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                              12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                              50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                              76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                              35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                            (42*pow(s + x - z,4)*z*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                      12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                      50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                      76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                      35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(-60*pow(l,10) + 6*c*pow(l,13)*r + 72*pow(l,11)*pow(r,2) - 330*pow(l,8)*s*z + 120*pow(l,8)*pow(z,2) -
                                                   24*c*pow(l,11)*r*pow(z,2) + 1640*pow(l,6)*s*pow(z,3) + 1200*pow(l,6)*pow(z,4) + 340*c*pow(l,9)*r*pow(z,4) -
                                                   1560*pow(l,7)*pow(r,2)*pow(z,4) - 540*pow(l,4)*s*pow(z,5) - 2520*pow(l,4)*pow(z,6) - 952*c*pow(l,7)*r*pow(z,6) +
                                                   3920*pow(l,5)*pow(r,2)*pow(z,6) - 2520*pow(l,2)*s*pow(z,7) + 1260*pow(l,2)*pow(z,8) + 630*c*pow(l,5)*r*pow(z,8) +
                                                   1750*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (7*pow(s + x - z,3)*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                             (-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                              40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                              68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                              136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                              70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6),2)) -
                            (42*pow(s + x - z,3)*z*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                      40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                      68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                      136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                      70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (pow(s + x - z,2)*((16*pow(r,2)*pow(z,3))/pow(pow(l,2) - pow(z,2),3) + (10*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2)))/
                            (pow(l,2) - pow(z,2)) + (2*pow(s + x - z,2)*z*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) +
                                                                                 pow(r,2)/(pow(l,2) - pow(z,2))))/pow(pow(l,2) - pow(z,2),2) -
                            (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)));
    }
    else if (label_i=='c'&&label_j=='c'){
        result = 0.;
    }
    else if (label_i=='c'&&label_j=='s'){
        result = (pow(l,2) - pow(s + x,2))*((2*r*(s + x - z))/(pow(l,2) - pow(z,2)) -
                                                (28*pow(s + x - z,3)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                                        35*pow(l,5)*r*pow(z,8)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (21*pow(s + x - z,2)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                                        70*pow(l,5)*r*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6)))) -
                 2*(s + x)*((r*pow(s + x - z,2))/(pow(l,2) - pow(z,2)) -
                            (7*pow(s + x - z,4)*(3*pow(l,13)*r - 12*pow(l,11)*r*pow(z,2) + 50*pow(l,9)*r*pow(z,4) - 76*pow(l,7)*r*pow(z,6) +
                                                   35*pow(l,5)*r*pow(z,8)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(6*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) + 68*pow(l,9)*r*pow(z,5) - 136*pow(l,7)*r*pow(z,7) +
                                                   70*pow(l,5)*r*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))));
    }
    else{
        result = (pow(l,2) - pow(s + x,2))*((-56*pow(s + x - z,3)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) -
                                                                         420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (84*pow(s + x - z,2)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                                        12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                                        50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                                        76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                                        35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*pow(s + x - z,2)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                                        315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                                                (42*(s + x - z)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                                 40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                                 68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                                 136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                                 70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                                                (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                                                (2*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2))) -
                 4*(s + x)*((2*pow(r,2)*z)/pow(pow(l,2) - pow(z,2),2) -
                            (7*pow(s + x - z,4)*(-105*pow(l,8)*z + 140*pow(l,6)*pow(z,3) + 210*pow(l,4)*pow(z,5) - 420*pow(l,2)*pow(z,7) + 175*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (28*pow(s + x - z,3)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                                    12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                                    50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                                    76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                                    35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (7*pow(s + x - z,3)*(-15*pow(l,10) - 165*pow(l,8)*pow(z,2) + 410*pow(l,6)*pow(z,4) - 90*pow(l,4)*pow(z,6) -
                                                   315*pow(l,2)*pow(z,8) + 175*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                            (21*pow(s + x - z,2)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                                    40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                                    68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                                    136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                                    70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                            (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                            (2*(s + x - z)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/(pow(l,2) - pow(z,2)))
                 - 2*((2*pow(r,2)*(s + x - z)*z)/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2)) -
                      (7*pow(s + x - z,4)*(-15*pow(l,10) + 3*c*pow(l,13)*r + 18*pow(l,11)*pow(r,2) - 105*pow(l,8)*s*z - 60*pow(l,8)*pow(z,2) -
                                             12*c*pow(l,11)*r*pow(z,2) + 90*pow(l,9)*pow(r,2)*pow(z,2) + 140*pow(l,6)*s*pow(z,3) + 270*pow(l,6)*pow(z,4) +
                                             50*c*pow(l,9)*r*pow(z,4) - 330*pow(l,7)*pow(r,2)*pow(z,4) + 210*pow(l,4)*s*pow(z,5) - 300*pow(l,4)*pow(z,6) -
                                             76*c*pow(l,7)*r*pow(z,6) + 350*pow(l,5)*pow(r,2)*pow(z,6) - 420*pow(l,2)*s*pow(z,7) + 105*pow(l,2)*pow(z,8) +
                                             35*c*pow(l,5)*r*pow(z,8) + 175*s*pow(z,9)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) -
                      (7*pow(s + x - z,3)*(-15*pow(l,10)*s - 60*pow(l,10)*z + 6*c*pow(l,13)*r*z + 72*pow(l,11)*pow(r,2)*z - 165*pow(l,8)*s*pow(z,2) +
                                             40*pow(l,8)*pow(z,3) - 8*c*pow(l,11)*r*pow(z,3) + 410*pow(l,6)*s*pow(z,4) + 240*pow(l,6)*pow(z,5) +
                                             68*c*pow(l,9)*r*pow(z,5) - 312*pow(l,7)*pow(r,2)*pow(z,5) - 90*pow(l,4)*s*pow(z,6) - 360*pow(l,4)*pow(z,7) -
                                             136*c*pow(l,7)*r*pow(z,7) + 560*pow(l,5)*pow(r,2)*pow(z,7) - 315*pow(l,2)*s*pow(z,8) + 140*pow(l,2)*pow(z,9) +
                                             70*c*pow(l,5)*r*pow(z,9) + 175*s*pow(z,10)))/
                      (pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) - 245*pow(z,6))) +
                      (pow(s + x - z,2)*(c*r + (4*pow(r,2)*pow(z,2))/pow(pow(l,2) - pow(z,2),2) + pow(r,2)/(pow(l,2) - pow(z,2))))/
                      (pow(l,2) - pow(z,2)));
    }
    return result;
    }
double shape::IntegrateRhoDerivative(double x, char label){
    double l,r,z,c,s;
    double result =0 ;
    l = _para_l;
    r = _para_r;
    z = _para_z;
    c = _para_c;
    s = _para_s;
    switch (label) {
        case 'l' :
            result = (x*(612*pow(l,16)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                         756*pow(l,11)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                          5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                          4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                         28*l*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                          252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                          210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                          420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                          35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                          15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                  21*pow(z,5))) - 300*pow(l,9)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                    280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                    70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                    35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                    84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                         168*pow(z,5))) + 360*pow(l,5)*pow(z,3)*
                                                                                                            (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                             70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                             7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                             28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                            93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                  360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                             2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                  42*x*pow(z,5) + 42*pow(z,6))) - 40*pow(l,7)*z*
                                                                                                                                                      (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                       420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                       14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                       84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                      334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                              1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                       8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                            378*x*pow(z,5) + 63*pow(z,6))) + 60*pow(l,3)*pow(z,5)*
                                                                                                                                                                                                 (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                  280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                  14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                  84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                 59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                       70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                  7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                       414*x*pow(z,5) + 84*pow(z,6))) + 140*pow(l,4)*r*pow(z,6)*
                                                                                                                                                                                                                                            (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                           105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                           35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                             15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                   35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                   7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                         540*pow(l,14)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                 10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                         84*pow(l,6)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                   70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                   140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                   42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                   pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                     6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                          35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                          35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                          7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                          pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                         156*pow(l,12)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                 200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                            c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                               100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                               15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                               5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                               3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                               s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                  396*pow(z,5)))) + 108*pow(l,8)*r*pow(z,2)*
                                                                      (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                             168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                             210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                             14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                             42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                       c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                     420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                     70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                     70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                     105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                          30*pow(z,5)))) + 88*pow(l,10)*r*(-(c*pow(z,2)*
                                                                                                                                 (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                  495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                  30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                  10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                  6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                  6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                       126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                     (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                       245*pow(z,6))) - (x*(54*pow(l,5) - 252*pow(l,3)*pow(z,2) - 42*l*pow(z,4))*
                                                                                           (36*pow(l,17)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                                                                                            105*s*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                                                                                              210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                                                                                              pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                                                                                              140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                                                                                            63*pow(l,12)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                                                                                            5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                                                                                            4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                                                                                            14*pow(l,2)*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                                                                                      252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                                                                                      210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                                                                                      420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                                                                                      35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                                                                                      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                                                                                              21*pow(z,5))) - 30*pow(l,10)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                                                                                                280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                                                                                                70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                                                                                                35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                                                                                                84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                                                                                                7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                                                                                                     168*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                                                                                                                                        (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                                                                                                         70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                                                                                                         7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                                                                                                         28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                                                                                                        93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                              360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                                                                                                         2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                                                                                              42*x*pow(z,5) + 42*pow(z,6))) - 5*pow(l,8)*z*
                                                                                                                                                                                                                                  (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                                                                                                   420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                                                                                                   14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                                                                                                   84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                                                                                                  334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                          1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                                                                                                   8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                        378*x*pow(z,5) + 63*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                                                                                                                                                                                                             (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                                                                                              280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                                                                                              14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                                                                                              84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                                                                                             59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                   70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                                                                                              7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                   414*x*pow(z,5) + 84*pow(z,6))) + 28*pow(l,5)*r*pow(z,6)*
                                                                                                                                                                                                                                                                                                                        (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                       105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                       35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                                                                                                         15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                               35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                               7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                                                                                            36*pow(l,15)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                                                                                   10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                                                                                   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                                                                                            12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                                                                                      70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                                                                                      140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                                                                                      42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                                                                                      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                                                                                        6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                                                                                             35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                                                                                             35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                                                                                             7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                                                                                             pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                                                                                            12*pow(l,13)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                                                                                   200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                                                                                              c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                                                                                                 100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                                                                                                 15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                                                                                                 5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                                                                                                 3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                                                                                                 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                                                                                    396*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                                                                                        (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                                               168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                               210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                                               14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                                               42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                                         c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                                       420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                       70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                                       70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                                       105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                                       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                            30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                                                                                                                  (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                   495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                                                                                   30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                                                                                   10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                                                                                   6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                                                                                   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                                                                                        126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                 450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                                                                                 45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                                                                                 27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                                                                                 9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                                                                                          (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                 245*pow(z,6),2)) - (x*(36*pow(l,17)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                                                                                                                                                                          105*s*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                                                                                                                                                                            210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                                                                                                                                                                            pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                                                                                                                                                                            140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                                                                                                                                                                          63*pow(l,12)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                                                                                                                                                                          5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                                                                                                                                                                          4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                                                                                                                                                                          14*pow(l,2)*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                                                                                                                                                                    252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                                                                                                                                                                    210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                                                                                                                                                                    420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                                                                                                                                                                    35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                                                                                                                                                                    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                                                                                                                                                                            21*pow(z,5))) - 30*pow(l,10)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                              280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                                                                                                                                                                              70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                                                                                                                                                                              35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                                                                                                                                                                              84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                                                                                                                                                                              7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                                                                                                                                                                                   168*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                                                                                                                                                                                                                      (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                                                                                                                                                                                       70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                                                                                                                                                                                       7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                                                                                                                                                                                       28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                                                                                                                                                                                      93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                            360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                                                                                                                                                                                       2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                                                                                                                                                                            42*x*pow(z,5) + 42*pow(z,6))) - 5*pow(l,8)*z*
                                                                                                                                                                                                                                                                                                                (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                 420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                 14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                 84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                                                                                                                                                                                334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                        1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                      378*x*pow(z,5) + 63*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                                                                                                                                                                                                                                                                                           (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                            280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                            14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                            84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                           59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                                                                 70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                                                            7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                 414*x*pow(z,5) + 84*pow(z,6))) + 28*pow(l,5)*r*pow(z,6)*
                                                                                                                                                                                                                                                                                                                                                                                                      (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                                     105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                                     35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                                                                                                                                                                                       15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                             35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                             7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                                                                                                                                                                          36*pow(l,15)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                                                                                                                                                                 10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                                                                                                                                                                          12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                                                                                                                                                                    70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                    140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                                                                                                                                                                    42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                                                                                                                                                                    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                                                                                                                                                                      6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                                                                                                                                                                           35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                                                                                                                                                                           35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                                                                                                                                                                           7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                                                                                                                                                                           pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                                                                                                                                                                          12*pow(l,13)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                                                                                                                                                                 200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                                                                                                                                                                            c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                                                                                                                                                                               100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                                                                                                                                                                               15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                                                                                                                                                                               5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                                                                                                                                                                               3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                                                                                                                                                                               s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                                                                                                                                                                  396*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                                                                                                                                                                      (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                             168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                             210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                                                                                                                             14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                                                                                                                             42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                                                                                                                       c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                     420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                                     70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                                     70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                                                                                                                     105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                                                                                                                     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                                                                                                          30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                                                                                                                                                                                                (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                                 495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                                                                                                                                                                 30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                                                                                                                                                                 10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                                                                                                                                                                 6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                                                                                                                                                                 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                                                                                                                                                                      126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                                                                                               450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                                                                                                                                                               45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                                                                                                                                                               27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                                                                                                                                                               9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                                                                                                                                                                      (2.*pow(l,4)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                       245*pow(z,6))) - (5*x*(36*pow(l,17)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                                                                                                                                                                                                                                                105*s*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                                                                                                                                                                                                                                                  210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                                                                                                                                                                                                                                                  pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                                                                                                                                                                                                                                                  140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                                                                                                                                                                                                                                                63*pow(l,12)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                                                                                                                                                                                                                                                5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                                                                                                                                                                                                                                                4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                                                                                                                                                                                                                                                14*pow(l,2)*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                                                                                                                                                                                                                                          252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                                                                                                                                                                                                                                          210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                                                                                                                                                                                                                                          420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                                                                                                                                                                                                                                          35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                                                                                                                                                                                                                                          15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                                                                                                                                                                                                                                                  21*pow(z,5))) - 30*pow(l,10)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                                                                                                    280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                                                                                                                                                                                                                                                    70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                    35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                    84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                                                                                                                                                                                                                                                         168*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                                                                                                                                                                                                                                                                                            (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                             70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                             7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                             28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                                                                                                                                                                                                                                                            93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                                                  360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                                                                                                                                                                                                                                                             2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                                                                                                                                                                                                                                                  42*x*pow(z,5) + 42*pow(z,6))) - 5*pow(l,8)*z*
                                                                                                                                                                                                                                                                                                                                                                                      (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                       420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                                                       14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                                                       84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                                                                                                                                                                                                                                                      334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                                                                                              1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                                                                                       8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                            378*x*pow(z,5) + 63*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                                                                                                                                                                                                                                                                                                                                                                 (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                  280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                  14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                  84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                 59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                  7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                                                                                       414*x*pow(z,5) + 84*pow(z,6))) + 28*pow(l,5)*r*pow(z,6)*
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                                                                                                                                                                                                                                                36*pow(l,15)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                                                                                                                                                                                                                                       10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                                                                                                                                                                                                                                                12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                                                                                                                                                                                                                                          70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                                                                                          140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                                                                                                                                                                                                                                          42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                                                                                                                                                                                                                                          pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                                                                                                                                                                                                                                            6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                                                                                                                                                                                                                                                 35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                                                                                                                                                                                                                                                 35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                                                                                                                                                                                                                                                 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                                                                                                                                                                                                                                                 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                                                                                                                                                                                                                                                12*pow(l,13)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                                                                                                                                                                                                                                       200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                                                                                                                                                                                                                                                  c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                     100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                                                                                                                                                                                                                                                     15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                                                                                                                                                                                                                                                     5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                                                                                                                                                                                                                                                     3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                                                                                                                                                                                                                                                     s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                                                                                                                                                                                                                                        396*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                                                                                                                                                                                                                                            (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                                                   168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                                                                                                   210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                                                                                                                                                                                                   14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                                                                                                                                                                                                   42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                                                                                                                                                                                             c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                                                           420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                                                                                                           70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                                                                                                           70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                                                                                                                                                                                           105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                                                                                                                                                                                           7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                                                                                                                                                                                30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                                                                                                                                                                                                                                                                      (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                                                                                                       495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                                                                                                                                                                                                                                       30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                       10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                       6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                       6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                            126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                                                                                                                                                                     450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                                                                                                                                                                                                                                     45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                                                     27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                                                     9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                                                                                                                                                                                                                                          (12.*pow(l,6)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                            245*pow(z,6)));
        case 'r' :
            result = (x*(108*pow(l,17)*r - 1080*pow(l,15)*r*pow(z,2) +
                         36*pow(l,17)*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                         420*pow(l,5)*r*pow(z,6)*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                      35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                      7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2))) -
                         36*pow(l,13)*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) - 200*pow(z,4) +
                                           90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) +
                         180*pow(l,9)*r*pow(z,2)*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z +
                                                      84*pow(x,4)*pow(z,2) - 168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                      210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                      14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                      42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) -
                         72*pow(l,7)*r*pow(z,4)*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) +
                                                     105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                     35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                     35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4))) +
                         24*pow(l,11)*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) + 450*pow(x,3)*pow(z,3) -
                                           765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) + 45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) +
                                           45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                           27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                           9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))) +
                         28*pow(l,5)*pow(z,6)*(c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                 105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                 35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                   15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                         35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                         7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                         36*pow(l,15)*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                              10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                              10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                         12*pow(l,7)*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                 70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                 140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                 42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                 pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                   6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                        35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                        35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                        7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                        pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                         12*pow(l,13)*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                              200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                         c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                            100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                            15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                            5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                            3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                            s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                               396*pow(z,5)))) + 12*pow(l,9)*pow(z,2)*
                                                                   (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                          168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                          210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                          14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                          42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                    c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                  420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                  70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                  70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                  105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                       30*pow(z,5)))) + 8*pow(l,11)*(-(c*pow(z,2)*
                                                                                                                           (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                            495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                            30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                            10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                            6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                            6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                 126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                          450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                          45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                          27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                          9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                     (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                       245*pow(z,6)));
        case 'z' :
            result = (x*(36*c*pow(l,17)*r*(-6*s - 3*x + 6*z) + 105*s*pow(z,9)*
                                                                 (-420*pow(s,5) - 210*pow(s,2)*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*pow(s,2)*x*(x - z)*(2*x - z) +
                                                                  70*pow(s,4)*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*pow(z,2)) +
                                                                  140*pow(s,3)*(-10*pow(x,2) + 12*x*z - 3*pow(z,2))) +
                         63*pow(l,12)*(-60*pow(s,3) + 60*pow(s,2)*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*pow(z,2)) +
                                         4*(-20*pow(x,2)*z + 60*x*pow(z,2) - 60*pow(z,3))) +
                         945*s*pow(z,8)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                           210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                           pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                           140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) -
                         14*pow(l,2)*pow(z,7)*(-8820*pow(s,6) - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*pow(z,2) +
                                                   350*pow(x,2)*pow(z,4) + 105*pow(s,5)*(-213*x + 232*z) + 210*pow(s,4)*(-145*pow(x,2) + 240*x*z - 126*pow(z,2)) +
                                                   420*pow(s,2)*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*pow(z,2)) -
                                                   840*pow(s,2)*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                   35*pow(s,3)*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*pow(z,2) + 432*pow(z,3)) +
                                                   15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*pow(z,2) + 336*x*pow(z,3) - 105*pow(z,4))) -
                         30*pow(l,10)*(-420*pow(s,5) - 168*pow(x,4)*z + 840*pow(x,3)*pow(z,2) - 1960*pow(x,2)*pow(z,3) +
                                         2520*x*pow(z,4) - 1680*pow(z,5) + 70*pow(s,4)*(-12*x + 12*z) +
                                         35*pow(s,3)*(-24*pow(x,2) + 24*x*z + 60*pow(z,2)) + 84*pow(s,2)*(-5*pow(x,3) + 45*x*pow(z,2) - 80*pow(z,3)) +
                                         7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*pow(z,2) - 900*x*pow(z,3) + 840*pow(z,4))) -
                         98*pow(l,2)*pow(z,6)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                   252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                   210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                   420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                   35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                   15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                           21*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                             (168*pow(s,6) + 7*pow(s,5)*(87*x - 576*z) + 70*pow(s,4)*(17*pow(x,2) - 132*x*z + 168*pow(z,2)) +
                                                                              7*pow(s,3)*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*pow(z,2) - 2400*pow(z,3)) +
                                                                              28*pow(s,2)*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*pow(z,2) - 720*x*pow(z,3) + 465*pow(z,4)) +
                                                                              7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*pow(z,2) - 1440*pow(x,2)*pow(z,3) + 1305*x*pow(z,4) -
                                                                                   504*pow(z,5)) - 2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*pow(z,2) - 700*pow(x,2)*pow(z,3) +
                                                                                                          210*x*pow(z,4) + 252*pow(z,5)) - 2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) +
                                                                                                                                                  210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) + 42*x*pow(z,5) + 42*pow(z,6))) +
                         180*pow(l,6)*pow(z,2)*(196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) +
                                                    7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                    70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                    7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                    28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                   93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                         360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                    2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                         42*x*pow(z,5) + 42*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                             (-6888*pow(s,6) + 84*pow(s,5)*(-215*x + 372*z) + 280*pow(s,4)*(-92*pow(x,2) + 240*x*z - 198*pow(z,2)) +
                                                                                              14*pow(s,3)*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*pow(z,2) + 3640*pow(z,3)) +
                                                                                              84*pow(s,2)*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*pow(z,2) + 660*x*pow(z,3) - 295*pow(z,4)) +
                                                                                              4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*pow(z,2) - 280*pow(x,2)*pow(z,3) + 294*pow(z,5)) +
                                                                                              7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*pow(z,2) + 3280*pow(x,2)*pow(z,3) - 2070*x*pow(z,4) +
                                                                                                   504*pow(z,5)) + 4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                        70*pow(x,2)*pow(z,4) + 49*pow(z,6))) -
                         5*pow(l,8)*z*(-3276*pow(s,6) + 42*pow(s,5)*(-183*x - 4*z) + 420*pow(s,4)*(-22*pow(x,2) - 12*x*z + 81*pow(z,2)) +
                                         14*pow(s,3)*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*pow(z,2) - 7920*pow(z,3)) +
                                         84*pow(s,2)*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*pow(z,2) - 1800*x*pow(z,3) + 1670*pow(z,4)) +
                                         21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*pow(z,2) - 4320*pow(x,2)*pow(z,3) + 5340*x*pow(z,4) -
                                               3024*pow(z,5)) + 8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*pow(z,2) + 2800*pow(x,2)*pow(z,3) -
                                                                       1890*x*pow(z,4) + 378*pow(z,5)) + 8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) -
                                                                                                                630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 63*pow(z,6))) -
                         5*pow(l,8)*(1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                       420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                       14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                       84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                      334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                              1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                       8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                            378*x*pow(z,5) + 63*pow(z,6))) + 75*pow(l,4)*pow(z,4)*
                                                                                 (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                  280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                  14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                  84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                 59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                       70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                  7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                       414*x*pow(z,5) + 84*pow(z,6))) - 36*pow(l,15)*r*
                                                                                                                            (60*r*z + c*(-20*pow(s,3) - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*pow(z,2) + 80*pow(z,3) + 10*pow(s,2)*(-3*x + 6*z) +
                                                                                                                                         10*s*(-2*pow(x,2) + 6*x*z - 12*pow(z,2)))) +
                         28*pow(l,5)*r*pow(z,6)*(c*pow(z,2)*(-210*pow(s,5) - 210*s*pow(x,3)*(x - z) + 105*pow(s,4)*(-5*x + 2*z) +
                                                                   105*pow(s,2)*pow(x,2)*(-5*x + 4*z) + 35*pow(s,3)*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
                                                     15*r*(-168*pow(s,5) - 420*pow(s,2)*pow(x,2)*(x - z) + 35*pow(s,4)*(-12*x + 6*z) + 70*pow(s,3)*x*(-8*x + 6*z) +
                                                           7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
                                                     2*c*z*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                            105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                            35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2)))) +
                         168*pow(l,5)*r*pow(z,5)*(c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                    105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                    35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                      15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                            35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                            7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                         12*pow(l,13)*r*(3*r*(-540*pow(s,2)*z - 180*pow(x,2)*z + 540*x*pow(z,2) - 800*pow(z,3) +
                                                45*s*(-12*x*z + 24*pow(z,2))) - c*(-42*pow(s,5) - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*pow(z,2) +
                                                                                     680*pow(x,2)*pow(z,3) - 990*x*pow(z,4) + 852*pow(z,5) + 15*pow(s,4)*(-7*x + 22*z) +
                                                                                     5*pow(s,3)*(-28*pow(x,2) + 132*x*z - 240*pow(z,2)) +
                                                                                     3*pow(s,2)*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*pow(z,2) + 680*pow(z,3)) +
                                                                                     s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*pow(z,2) + 2040*x*pow(z,3) - 1980*pow(z,4)))) -
                         12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(-1176*pow(s,5) + 70*pow(s,4)*(-42*x + 30*z) +
                                                                   14*s*x*pow(x - z,2)*(-8*x + 42*z) + 140*pow(s,3)*(-28*pow(x,2) + 30*x*z - 15*pow(z,2)) -
                                                                   28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                   42*pow(s,2)*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*pow(z,2) + 28*pow(z,3)) +
                                                                   pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*pow(z,2) + 392*pow(z,3))) +
                                                     6*r*(-1176*pow(s,5) + 105*pow(s,2)*pow(x - z,2)*(-6*x + 14*z) + 35*pow(s,4)*(-84*x + 90*z) +
                                                          35*pow(s,3)*(-112*pow(x,2) + 180*x*z - 120*pow(z,2)) -
                                                          210*pow(s,2)*(x - z)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                          7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*pow(z,2) + 420*pow(z,3)) +
                                                          pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*pow(z,2) + 980*pow(z,3))) +
                                                     2*c*z*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) + 70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) +
                                                            14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                            140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                            42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                            pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4)))) -
                         48*pow(l,7)*r*pow(z,3)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                   70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                   140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                   42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                   pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                     6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                          35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                          35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                          7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                          pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) +
                         8*pow(l,11)*r*(-(c*pow(z,2)*(-420*pow(s,5) - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*pow(z,2) +
                                                          2160*pow(x,2)*pow(z,3) - 1890*x*pow(z,4) + 1260*pow(z,5) + 30*pow(s,4)*(-35*x + 80*z) +
                                                          10*pow(s,3)*(-140*pow(x,2) + 480*x*z - 594*pow(z,2)) +
                                                          6*pow(s,2)*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*pow(z,2) + 1080*pow(z,3)) +
                                                          6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*pow(z,2) + 1080*x*pow(z,3) - 630*pow(z,4)))) +
                                          3*r*(-1350*pow(s,4)*z - 270*pow(x,4)*z + 1350*pow(x,3)*pow(z,2) - 3060*pow(x,2)*pow(z,3) + 3780*x*pow(z,4) -
                                               3570*pow(z,5) + 45*pow(s,3)*(-60*x*z + 120*pow(z,2)) +
                                               27*pow(s,2)*(-100*pow(x,2)*z + 300*x*pow(z,2) - 340*pow(z,3)) +
                                               9*s*(-150*pow(x,3)*z + 600*pow(x,2)*pow(z,2) - 1020*x*pow(z,3) + 840*pow(z,4))) -
                                          2*c*z*(126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                 495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                 30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                 10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                 6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                      126*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                          (15*r*(-168*pow(s,5) - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*pow(z,2) + 728*pow(x,2)*pow(z,3) -
                                                                                 420*x*pow(z,4) + 294*pow(z,5) + 210*pow(s,4)*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
                                                                                 14*pow(s,3)*(-40*pow(x,2) + 120*x*z - 144*pow(z,2)) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                 42*pow(s,2)*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*pow(z,2) + 52*pow(z,3))) +
                                                                           c*pow(z,2)*(-924*pow(s,5) - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*pow(z,2) + 1260*pow(x,2)*pow(z,3) -
                                                                                         525*x*pow(z,4) + 294*pow(z,5) + 70*pow(s,4)*(-33*x + 42*z) +
                                                                                         70*pow(s,3)*(-44*pow(x,2) + 84*x*z - 72*pow(z,2)) +
                                                                                         105*pow(s,2)*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*pow(z,2) + 36*pow(z,3)) +
                                                                                         7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*pow(z,2) + 540*x*pow(z,3) - 150*pow(z,4))) +
                                                                           2*c*z*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                  420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                  70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                  70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                  105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                       30*pow(z,5)))) + 24*pow(l,9)*r*z*(15*r*
                                                                                                                             (42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                              168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                              210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                              14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                              42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                             c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                           420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                           70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                           70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                           105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                           7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                30*pow(z,5))))))/
                     (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                       245*pow(z,6))) - (x*(-126*pow(l,4)*z - 84*pow(l,2)*pow(z,3) - 1470*pow(z,5))*
                                                                                           (36*pow(l,17)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                                                                                            105*s*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                                                                                              210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                                                                                              pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                                                                                              140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                                                                                            63*pow(l,12)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                                                                                            5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                                                                                            4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                                                                                            14*pow(l,2)*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                                                                                      252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                                                                                      210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                                                                                      420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                                                                                      35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                                                                                      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                                                                                              21*pow(z,5))) - 30*pow(l,10)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                                                                                                280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                                                                                                70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                                                                                                35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                                                                                                84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                                                                                                7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                                                                                                     168*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                                                                                                                                        (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                                                                                                         70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                                                                                                         7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                                                                                                         28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                                                                                                        93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                              360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                                                                                                         2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                                                                                              42*x*pow(z,5) + 42*pow(z,6))) - 5*pow(l,8)*z*
                                                                                                                                                                                                                                  (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                                                                                                   420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                                                                                                   14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                                                                                                   84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                                                                                                  334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                          1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                                                                                                   8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                        378*x*pow(z,5) + 63*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                                                                                                                                                                                                             (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                                                                                              280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                                                                                              14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                                                                                              84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                                                                                             59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                   70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                                                                                              7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                   414*x*pow(z,5) + 84*pow(z,6))) + 28*pow(l,5)*r*pow(z,6)*
                                                                                                                                                                                                                                                                                                                        (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                       105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                       35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                                                                                                         15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                               35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                               7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                                                                                            36*pow(l,15)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                                                                                   10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                                                                                   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                                                                                            12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                                                                                      70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                                                                                      140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                                                                                      42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                                                                                      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                                                                                        6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                                                                                             35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                                                                                             35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                                                                                             7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                                                                                             pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                                                                                            12*pow(l,13)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                                                                                   200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                                                                                              c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                                                                                                 100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                                                                                                 15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                                                                                                 5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                                                                                                 3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                                                                                                 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                                                                                    396*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                                                                                        (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                                               168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                               210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                                               14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                                               42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                                         c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                                       420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                       70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                                       70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                                       105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                                       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                            30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                                                                                                                  (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                   495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                                                                                   30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                                                                                   10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                                                                                   6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                                                                                   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                                                                                        126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                 450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                                                                                 45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                                                                                 27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                                                                                 9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                                                                                          (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*pow(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                 245*pow(z,6),2)) + (x*z*(36*pow(l,17)*r*(3*r + c*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2))) +
                                                                                                                                                                            105*s*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                                                                                                                                                                              210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                                                                                                                                                                              pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                                                                                                                                                                              140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                                                                                                                                                                            63*pow(l,12)*(40*pow(s,4) + 10*pow(s,3)*(7*x - 6*z) + 60*pow(s,2)*(pow(x,2) - x*z - pow(z,2)) +
                                                                                                                                                                                            5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3)) +
                                                                                                                                                                                            4*(pow(x,4) - 10*pow(x,2)*pow(z,2) + 20*x*pow(z,3) - 15*pow(z,4))) -
                                                                                                                                                                            14*pow(l,2)*pow(z,7)*(2520*pow(s,7) + 1260*pow(s,6)*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*pow(z,2) -
                                                                                                                                                                                                      252*pow(x,4)*pow(z,3) + 70*pow(x,2)*pow(z,5) + 105*pow(s,5)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                                                                                                                                                                      210*pow(s,4)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                                                                                                                                                                      420*pow(s,2)*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                                                                                                                                                                      35*pow(s,3)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                                                                                                                                                                      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                                                                                                                                                                              21*pow(z,5))) - 30*pow(l,10)*(84*pow(s,6) + 6*pow(x,6) + 21*pow(s,5)*(11*x - 20*z) - 84*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                                280*pow(x,3)*pow(z,3) - 490*pow(x,2)*pow(z,4) + 504*x*pow(z,5) - 280*pow(z,6) +
                                                                                                                                                                                                                                                70*pow(s,4)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                                                                                                                                                                                35*pow(s,3)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                                                                                                                                                                                84*pow(s,2)*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                                                                                                                                                                                7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) +
                                                                                                                                                                                                                                                     168*pow(z,5))) + 60*pow(l,6)*pow(z,3)*
                                                                                                                                                                                                                                                                        (196*pow(s,7) + 84*pow(s,6)*(7*x + 2*z) + 7*pow(s,5)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                                                                                                                                                                                                                                         70*pow(s,4)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                                                                                                                                                                                                                                         7*pow(s,3)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                                                                                                                                                                                                                                         28*pow(s,2)*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                                                                                                                                                                                                                                                        93*pow(z,5)) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                              360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6)) -
                                                                                                                                                                                                                                                                         2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*pow(z,2) + 210*pow(x,3)*pow(z,3) - 175*pow(x,2)*pow(z,4) +
                                                                                                                                                                                                                                                                              42*x*pow(z,5) + 42*pow(z,6))) - 5*pow(l,8)*z*
                                                                                                                                                                                                                                                                                                                  (1764*pow(s,7) + 252*pow(s,6)*(21*x - 13*z) + 42*pow(s,5)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                   420*pow(s,4)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                   14*pow(s,3)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                   84*pow(s,2)*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                                                                                                                                                                                                                                                                                                                  334*pow(z,5)) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                          1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                   8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*pow(z,2) - 630*pow(x,3)*pow(z,3) + 700*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                        378*x*pow(z,5) + 63*pow(z,6))) + 15*pow(l,4)*pow(z,5)*
                                                                                                                                                                                                                                                                                                                                                             (1176*pow(s,7) + 168*pow(s,6)*(21*x - 41*z) + 84*pow(s,5)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                              280*pow(s,4)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                                                              14*pow(s,3)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                                                              84*pow(s,2)*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                             59*pow(z,5)) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*pow(z,2) + 140*pow(x,3)*pow(z,3) -
                                                                                                                                                                                                                                                                                                                                                                                                   70*pow(x,2)*pow(z,4) + 49*pow(z,6)) +
                                                                                                                                                                                                                                                                                                                                                              7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) + 820*pow(x,2)*pow(z,4) -
                                                                                                                                                                                                                                                                                                                                                                   414*x*pow(z,5) + 84*pow(z,6))) + 28*pow(l,5)*r*pow(z,6)*
                                                                                                                                                                                                                                                                                                                                                                                                        (c*pow(z,2)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                                       105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                                       35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) +
                                                                                                                                                                                                                                                                                                                                                                                                         15*r*(70*pow(s,6) + 42*pow(s,5)*(5*x - 4*z) + 210*pow(s,2)*pow(x,2)*pow(x - z,2) +
                                                                                                                                                                                                                                                                                                                                                                                                               35*pow(s,4)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 70*pow(s,3)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                                                                                                                                                                                                               7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*pow(z,2)))) -
                                                                                                                                                                            36*pow(l,15)*r*(30*r*pow(z,2) + c*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z +
                                                                                                                                                                                                                   10*pow(x,2)*pow(z,2) - 20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                                                                                                                                                                   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3)))) -
                                                                                                                                                                            12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                                                                                                                                                                                      70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                      140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                                                                                                                                                                      42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                                                                                                                                                                                      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) +
                                                                                                                                                                                                        6*r*(385*pow(s,6) + 21*pow(s,5)*(55*x - 56*z) + 105*pow(s,2)*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                                                                                                                                                                             35*pow(s,4)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                                                                                                                                                                             35*pow(s,3)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                                                                                                                                                                             7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)) +
                                                                                                                                                                                                             pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*pow(z,2) - 350*x*pow(z,3) + 245*pow(z,4)))) -
                                                                                                                                                                            12*pow(l,13)*r*(3*r*(45*pow(s,4) + 90*pow(s,3)*x + 9*pow(x,4) - 90*pow(x,2)*pow(z,2) + 180*x*pow(z,3) -
                                                                                                                                                                                                   200*pow(z,4) + 90*pow(s,2)*(pow(x,2) - 3*pow(z,2)) + 45*s*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                                                                                                                                                                              c*(21*pow(s,6) + 3*pow(x,6) + 21*pow(s,5)*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                 100*pow(x,3)*pow(z,3) + 170*pow(x,2)*pow(z,4) - 198*x*pow(z,5) + 142*pow(z,6) +
                                                                                                                                                                                                 15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                                                                                                                                                                                 5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                                                                                                                                                                                 3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) +
                                                                                                                                                                                                 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) -
                                                                                                                                                                                                    396*pow(z,5)))) + 12*pow(l,9)*r*pow(z,2)*
                                                                                                                                                                                                                        (15*r*(42*pow(s,6) + 6*pow(x,6) + 42*pow(s,5)*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                               168*pow(x,3)*pow(z,3) + 182*pow(x,2)*pow(z,4) - 84*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                               210*pow(s,4)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                                                                                                                                                                                               14*pow(s,3)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                                                                                                                                                                                               42*pow(s,2)*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                                                                                                                                                                                         c*pow(z,2)*(350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                       420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                                                                                                                                                                       70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                                                                                                                                                                       70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                                                                                                                                                                       105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                                                                                                                                                                       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                                                                                                                                                                            30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                                                                                                                                                                                                  (126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*pow(z,2) -
                                                                                                                                                                                                                                                                                   495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                                                                                                                                                                                                                                                   30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                                                                                                                                                                                                   10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                                                                                                                                                                                                   6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                                                                                                                                                                                                   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                                                                                                                                                                                        126*pow(z,5)))) + 3*r*(63*pow(s,6) + 189*pow(s,5)*x + 9*pow(x,6) - 135*pow(x,4)*pow(z,2) +
                                                                                                                                                                                                                                                                                                                 450*pow(x,3)*pow(z,3) - 765*pow(x,2)*pow(z,4) + 756*x*pow(z,5) - 595*pow(z,6) +
                                                                                                                                                                                                                                                                                                                 45*pow(s,4)*(7*pow(x,2) - 15*pow(z,2)) + 45*pow(s,3)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                                                                                                                                                                                                 27*pow(s,2)*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                                                                                                                                                                                                 9*s*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                                                                                                                                                                      (2.*pow(l,5)*pow(pow(l,2) - pow(z,2),4)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                                                                                                                                                                       245*pow(z,6)));
        case 'c' :
            result = (x*(36*pow(l,17)*r*(3*pow(s,2) + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*pow(z,2)) +
                         28*pow(l,5)*r*pow(z,8)*(105*pow(s,6) + 105*pow(s,5)*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
                                                     105*pow(s,4)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 105*pow(s,2)*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                     35*pow(s,3)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2)) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*pow(z,2))) -
                         36*pow(l,15)*r*(10*pow(s,4) + 2*pow(x,4) + 20*pow(s,3)*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*pow(z,2) -
                                           20*x*pow(z,3) + 20*pow(z,4) + 10*pow(s,2)*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                           10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3))) -
                         12*pow(l,7)*r*pow(z,6)*(532*pow(s,6) + 84*pow(s,5)*(19*x - 14*z) +
                                                     70*pow(s,4)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                     140*pow(s,3)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                     42*pow(s,2)*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4)) +
                                                     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*pow(z,2) - 175*x*pow(z,3) + 98*pow(z,4))) -
                         12*pow(l,13)*r*(-21*pow(s,6) - 3*pow(x,6) - 21*pow(s,5)*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*pow(z,2) +
                                           100*pow(x,3)*pow(z,3) - 170*pow(x,2)*pow(z,4) + 198*x*pow(z,5) - 142*pow(z,6) -
                                           15*pow(s,4)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) -
                                           5*pow(s,3)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) -
                                           3*pow(s,2)*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)) -
                                           s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) - 400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) - 396*pow(z,5))
                         ) - 8*pow(l,11)*r*pow(z,2)*(126*pow(s,6) + 18*pow(x,6) + 42*pow(s,5)*(9*x - 10*z) - 70*pow(x,5)*z +
                                                         240*pow(x,4)*pow(z,2) - 495*pow(x,3)*pow(z,3) + 540*pow(x,2)*pow(z,4) - 378*x*pow(z,5) + 210*pow(z,6) +
                                                         30*pow(s,4)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                         10*pow(s,3)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                         6*pow(s,2)*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                         6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                              126*pow(z,5))) + 12*pow(l,9)*r*pow(z,4)*
                                                                                 (350*pow(s,6) + 50*pow(x,6) + 42*pow(s,5)*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*pow(z,2) -
                                                                                  420*pow(x,3)*pow(z,3) + 315*pow(x,2)*pow(z,4) - 105*x*pow(z,5) + 49*pow(z,6) +
                                                                                  70*pow(s,4)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                                  70*pow(s,3)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                                  105*pow(s,2)*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                                  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                                       30*pow(z,5)))))/
                     (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                       245*pow(z,6)));
        case 's' :
            result = (x*(36*c*pow(l,17)*r*(6*s + 3*x - 6*z) - 36*c*pow(l,15)*r*
                                                                (40*pow(s,3) + 60*pow(s,2)*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*pow(z,2)) +
                                                                 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*pow(z,2) - 4*pow(z,3))) +
                         105*s*pow(z,9)*(840*pow(s,5) + 2100*pow(s,4)*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
                                           280*pow(s,3)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                           420*pow(s,2)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                         105*pow(z,9)*(140*pow(s,6) + 420*pow(s,5)*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                         210*pow(s,2)*x*pow(x - z,2)*(2*x - z) + 70*pow(s,4)*(10*pow(x,2) - 15*x*z + 6*pow(z,2)) +
                                         pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*pow(z,2) - 35*pow(z,3)) +
                                         140*pow(s,3)*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*pow(z,2) - pow(z,3))) +
                         63*pow(l,12)*(160*pow(s,3) + 30*pow(s,2)*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - pow(z,2)) +
                                         5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*pow(z,2) + 28*pow(z,3))) -
                         14*pow(l,2)*pow(z,7)*(17640*pow(s,6) + 7560*pow(s,5)*(6*x - 7*z) +
                                                   525*pow(s,4)*(120*pow(x,2) - 213*x*z + 116*pow(z,2)) +
                                                   840*pow(s,3)*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*pow(z,2) - 42*pow(z,3)) +
                                                   840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*pow(z,2) - 2*pow(z,3)) +
                                                   105*pow(s,2)*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*pow(z,2) - 396*x*pow(z,3) + 108*pow(z,4)) +
                                                   15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*pow(z,2) - 210*pow(x,2)*pow(z,3) + 84*x*pow(z,4) -
                                                         21*pow(z,5))) - 30*pow(l,10)*(504*pow(s,5) + 105*pow(s,4)*(11*x - 20*z) +
                                                                                           280*pow(s,3)*(5*pow(x,2) - 12*x*z + 6*pow(z,2)) +
                                                                                           105*pow(s,2)*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*pow(z,2) + 20*pow(z,3)) +
                                                                                           168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*pow(z,3) - 20*pow(z,4)) +
                                                                                           7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*pow(z,2) + 140*pow(x,2)*pow(z,3) - 225*x*pow(z,4) + 168*pow(z,5)))\
        - 5*pow(l,8)*z*(12348*pow(s,6) + 1512*pow(s,5)*(21*x - 13*z) + 210*pow(s,4)*(210*pow(x,2) - 183*x*z - 2*pow(z,2)) +
                          1680*pow(s,3)*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*pow(z,2) + 27*pow(z,3)) +
                          42*pow(s,2)*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*pow(z,2) + 1665*x*pow(z,3) - 1980*pow(z,4)) +
                          168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*pow(z,2) + 285*pow(x,2)*pow(z,3) - 450*x*pow(z,4) +
                                 334*pow(z,5)) + 21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*pow(z,2) + 585*pow(x,3)*pow(z,3) -
                                                       1080*pow(x,2)*pow(z,4) + 1068*x*pow(z,5) - 504*pow(z,6))) +
                         60*pow(l,6)*pow(z,3)*(1372*pow(s,6) + 504*pow(s,5)*(7*x + 2*z) +
                                                   35*pow(s,4)*(140*pow(x,2) + 87*x*z - 288*pow(z,2)) +
                                                   280*pow(s,3)*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*pow(z,2) + 56*pow(z,3)) +
                                                   21*pow(s,2)*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*pow(z,2) + 960*x*pow(z,3) - 600*pow(z,4)) +
                                                   56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 180*x*pow(z,4) +
                                                         93*pow(z,5)) + 7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*pow(z,2) + 320*pow(x,3)*pow(z,3) -
                                                                             360*pow(x,2)*pow(z,4) + 261*x*pow(z,5) - 84*pow(z,6))) +
                         15*pow(l,4)*pow(z,5)*(8232*pow(s,6) + 1008*pow(s,5)*(21*x - 41*z) +
                                                   420*pow(s,4)*(70*pow(x,2) - 215*x*z + 186*pow(z,2)) +
                                                   1120*pow(s,3)*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*pow(z,2) - 66*pow(z,3)) +
                                                   42*pow(s,2)*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*pow(z,2) - 2145*x*pow(z,3) + 910*pow(z,4)) +
                                                   168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*pow(z,2) - 275*pow(x,2)*pow(z,3) + 165*x*pow(z,4) -
                                                          59*pow(z,5)) + 7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*pow(z,2) - 1155*pow(x,3)*pow(z,3) +
                                                                              820*pow(x,2)*pow(z,4) - 414*x*pow(z,5) + 84*pow(z,6))) +
                         28*pow(l,5)*r*pow(z,6)*(c*pow(z,2)*(630*pow(s,5) + 525*pow(s,4)*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) +
                                                                   420*pow(s,3)*(5*pow(x,2) - 5*x*z + pow(z,2)) + 210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*pow(z,2)) +
                                                                   105*pow(s,2)*x*(15*pow(x,2) - 20*x*z + 6*pow(z,2))) +
                                                     15*r*(420*pow(s,5) + 210*pow(s,4)*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) +
                                                           140*pow(s,3)*(10*pow(x,2) - 12*x*z + 3*pow(z,2)) + 210*pow(s,2)*x*(5*pow(x,2) - 8*x*z + 3*pow(z,2)) +
                                                           7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*pow(z,2)))) -
                         12*pow(l,7)*r*pow(z,4)*(c*pow(z,2)*(3192*pow(s,5) + 420*pow(s,4)*(19*x - 14*z) +
                                                                   280*pow(s,3)*(38*pow(x,2) - 42*x*z + 15*pow(z,2)) + 14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*pow(z,2)) +
                                                                   420*pow(s,2)*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*pow(z,2) - 5*pow(z,3)) +
                                                                   84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*pow(z,2) - 25*x*pow(z,3) + 7*pow(z,4))) +
                                                     6*r*(2310*pow(s,5) + 105*pow(s,4)*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*pow(z,2)) +
                                                          140*pow(s,3)*(55*pow(x,2) - 84*x*z + 45*pow(z,2)) +
                                                          105*pow(s,2)*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*pow(z,2) - 40*pow(z,3)) +
                                                          7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 105*pow(z,4)))) -
                         12*pow(l,13)*r*(3*r*(180*pow(s,3) + 270*pow(s,2)*x + 180*s*(pow(x,2) - 3*pow(z,2)) +
                                                45*(pow(x,3) - 6*x*pow(z,2) + 8*pow(z,3))) -
                                           c*(126*pow(s,5) + 21*pow(x,5) + 105*pow(s,4)*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*pow(z,2) -
                                              400*pow(x,2)*pow(z,3) + 510*x*pow(z,4) - 396*pow(z,5) + 60*pow(s,3)*(7*pow(x,2) - 7*x*z + 11*pow(z,2)) +
                                              15*pow(s,2)*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*pow(z,2) - 80*pow(z,3)) +
                                              6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*pow(z,2) - 200*x*pow(z,3) + 170*pow(z,4)))) +
                         12*pow(l,9)*r*pow(z,2)*(15*r*(252*pow(s,5) + 210*pow(s,4)*(3*x - 4*z) +
                                                           840*pow(s,3)*(pow(x,2) - 2*x*z + 2*pow(z,2)) + 42*pow(x - z,3)*(pow(x,2) - x*z + 4*pow(z,2)) +
                                                           42*pow(s,2)*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*pow(z,2) - 48*pow(z,3)) +
                                                           84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 13*pow(z,4))) +
                                                     c*pow(z,2)*(2100*pow(s,5) + 210*pow(s,4)*(25*x - 22*z) + 280*pow(s,3)*(25*pow(x,2) - 33*x*z + 21*pow(z,2)) +
                                                                   210*pow(s,2)*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*pow(z,2) - 24*pow(z,3)) +
                                                                   210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*pow(z,2) - 24*x*pow(z,3) + 9*pow(z,4)) +
                                                                   7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*pow(z,2) - 240*pow(x,2)*pow(z,3) + 135*x*pow(z,4) -
                                                                      30*pow(z,5)))) + 8*pow(l,11)*r*(-(c*pow(z,2)*
                                                                                                            (756*pow(s,5) + 210*pow(s,4)*(9*x - 10*z) + 120*pow(s,3)*(21*pow(x,2) - 35*x*z + 40*pow(z,2)) +
                                                                                                             30*pow(s,2)*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*pow(z,2) - 198*pow(z,3)) +
                                                                                                             12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*pow(z,2) - 495*x*pow(z,3) + 270*pow(z,4)) +
                                                                                                             6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*pow(z,2) - 330*pow(x,2)*pow(z,3) + 270*x*pow(z,4) -
                                                                                                                126*pow(z,5)))) + 3*r*(378*pow(s,5) + 945*pow(s,4)*x + 180*pow(s,3)*(7*pow(x,2) - 15*pow(z,2)) +
                                                                                                                                         135*pow(s,2)*(7*pow(x,3) - 30*x*pow(z,2) + 40*pow(z,3)) +
                                                                                                                                         54*s*(7*pow(x,4) - 50*pow(x,2)*pow(z,2) + 100*x*pow(z,3) - 85*pow(z,4)) +
                                                                                                                                         9*(7*pow(x,5) - 75*pow(x,3)*pow(z,2) + 200*pow(x,2)*pow(z,3) - 255*x*pow(z,4) + 168*pow(z,5))))))/
                     (12.*pow(l,5)*pow(pow(l,2) - pow(z,2),3)*(9*pow(l,6) - 63*pow(l,4)*pow(z,2) - 21*pow(l,2)*pow(z,4) -
                                                                       245*pow(z,6)));

    }
    return result;
}
double shape::IntegrateRhoDDerivative(double zeta, char label_i, char label_j) {
    double l,r,z,c,s;
    double result = 0;
    l = _para_l;
    r = _para_r;
    z = _para_z;
    c = _para_c;
    s = _para_s;
    if (label_i=='l'&&label_j=='l'){
        result =
    }
}
double RhoShape(double zeta, void *params) {
    shape* tmp_shape = (shape*) params;
    double result;
    result = tmp_shape->Rho(zeta);
    return result;
}
double shape::CenterOfMassDerivative(const char side, const char label) {
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;
    double result = 0;
    switch (side) {
        case 'L':
            switch (label){
                case 'l' :
                    result = ((pow(l,2) - pow(z,2))*(-1050*pow(l,9) + 39*c*pow(l,12)*r - 210*l*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                         840*pow(l,3)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                                         1260*pow(l,5)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 840*pow(l,7)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                                         14*pow(l,6)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                         35*pow(l,4)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 66*pow(l,10)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                         18*pow(l,8)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                                  7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                  12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                  12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))) -
                             ((pow(l,2) - pow(z,2))*(-1050*pow(l,9)*s + 156*c*pow(l,12)*r*z + 280*pow(l,7)*(15*s - 16*z)*pow(z,2) +
                                                         24*pow(l,3)*(175*s - 184*z)*pow(z,6) - 14*l*(75*s - 16*z)*pow(z,8) - 252*pow(l,5)*pow(z,4)*(25*s + 24*z) +
                                                         44*pow(l,10)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 140*pow(l,4)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                                         108*pow(l,8)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                                         84*pow(l,6)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) +
                                       6*pow(l,4)*(175*s - 184*z)*pow(z,6) - 7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) -
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) -
                                       28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                       12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                       12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)),2)) +
                             (l*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                                 105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                 210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                 210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                 2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                 7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                 2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                              7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                              4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                              12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                              12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)));
                case 'r' :
                    result = ((pow(l,2) - pow(z,2))*(3*c*pow(l,13) + 54*pow(l,11)*r - 54*pow(l,9)*r*(4*s - z)*z - 18*pow(l,7)*r*(8*s - 49*z)*pow(z,3) -
                                                         630*pow(l,5)*r*(4*s - 3*z)*pow(z,5) + 2*pow(l,7)*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                         7*pow(l,5)*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                         2*pow(l,9)*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                                  7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                  12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                  12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))) -
                             ((pow(l,2) - pow(z,2))*(12*c*pow(l,13)*z + 4*pow(l,11)*(54*r*z - 2*c*pow(z,3)) -
                                                         28*pow(l,5)*pow(z,6)*(90*r*z + c*pow(z,3)) - 12*pow(l,9)*pow(z,2)*(6*r*z + 4*c*pow(z,3)) +
                                                         12*pow(l,7)*pow(z,4)*(198*r*z + 6*c*pow(z,3)))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) +
                                       6*pow(l,4)*(175*s - 184*z)*pow(z,6) - 7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) -
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) -
                                       28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                       12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                       12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)),2));
                case 'z' :
                    result = ((pow(l,2) - pow(z,2))*(840*s*(2*s - z)*pow(z,7) - 105*s*pow(z,8) - 105*pow(l,2)*pow(z,6)*(-4*s + 2*z) +
                                                         210*pow(l,4)*pow(z,4)*(-3*s + 4*z) - 210*pow(l,6)*pow(z,2)*(-2*s + 6*z) + 105*pow(l,8)*(-s + 8*z) -
                                                         630*pow(l,2)*pow(z,5)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                         840*pow(l,4)*pow(z,3)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) - 420*pow(l,6)*z*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) +
                                                         2*pow(l,7)*r*pow(z,3)*(441*r + 4*c*(22*s - 7*z)*z - 14*c*pow(z,2)) -
                                                         7*pow(l,5)*r*pow(z,5)*(-270*r + 2*c*(8*s - 5*z)*z - 5*c*pow(z,2)) +
                                                         6*pow(l,7)*r*pow(z,2)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                         35*pow(l,5)*r*pow(z,4)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(6*c*z + 2*c*(-2*s + 3*z)) -
                                                         2*pow(l,9)*r*z*(-27*r + 23*c*pow(z,2) + 2*c*z*(4*s + 23*z)) - 2*pow(l,9)*r*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                                  7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                  12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                  12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))) -
                             ((pow(l,2) - pow(z,2))*(12*c*pow(l,13)*r + 70*pow(l,8)*(15*s - 16*z)*z - 560*pow(l,8)*pow(z,2) -
                                                         1008*pow(l,6)*pow(z,4) + 36*pow(l,4)*(175*s - 184*z)*pow(z,5) - 1104*pow(l,4)*pow(z,6) -
                                                         56*pow(l,2)*(75*s - 16*z)*pow(z,7) + 112*pow(l,2)*pow(z,8) + 1050*s*pow(z,9) -
                                                         168*pow(l,6)*pow(z,3)*(25*s + 24*z) + 4*pow(l,11)*(27*pow(r,2) - 6*c*r*pow(z,2)) -
                                                         28*pow(l,5)*pow(z,6)*(45*pow(r,2) + 3*c*r*pow(z,2)) - 12*pow(l,9)*pow(z,2)*(3*pow(r,2) + 12*c*r*pow(z,2)) +
                                                         12*pow(l,7)*pow(z,4)*(99*pow(r,2) + 18*c*r*pow(z,2)) -
                                                         168*pow(l,5)*pow(z,5)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) - 24*pow(l,9)*z*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                                         48*pow(l,7)*pow(z,3)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) +
                                       6*pow(l,4)*(175*s - 184*z)*pow(z,6) - 7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) -
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) -
                                       28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                       12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                       12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)),2)) -
                             (z*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                                 105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                 210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                 210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                 2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                 7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                 2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                              7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                              4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                              12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                              12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)));
                case 'c' :
                    result = ((pow(l,2) - pow(z,2))*(3*pow(l,13)*r + 4*pow(l,7)*r*(22*s - 7*z)*pow(z,5) - 7*pow(l,5)*r*(8*s - 5*z)*pow(z,7) +
                                                         12*pow(l,11)*r*z*(-2*s + 3*z) - 2*pow(l,9)*r*pow(z,3)*(4*s + 23*z)))/
                             (2.*(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                                  7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                  12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                  12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))) -
                             ((pow(l,2) - pow(z,2))*(12*pow(l,13)*r*z - 8*pow(l,11)*r*pow(z,3) - 48*pow(l,9)*r*pow(z,5) +
                                                         72*pow(l,7)*r*pow(z,7) - 28*pow(l,5)*r*pow(z,9))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) +
                                       6*pow(l,4)*(175*s - 184*z)*pow(z,6) - 7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) -
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) -
                                       28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                       12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                       12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)),2));
                case 's' :
                    result = ((pow(l,2) - pow(z,2))*(105*pow(l,8)*(4*s - z) - 24*c*pow(l,11)*r*z - 210*pow(l,6)*(8*s - 2*z)*pow(z,2) +
                                                         210*pow(l,4)*(12*s - 3*z)*pow(z,4) - 105*pow(l,2)*(16*s - 4*z)*pow(z,6) + 210*s*pow(z,8) + 105*(2*s - z)*pow(z,8) -
                                                         2*pow(l,9)*r*z*(108*r + 4*c*pow(z,2)) - 7*pow(l,5)*r*pow(z,5)*(360*r + 8*c*pow(z,2)) +
                                                         2*pow(l,7)*r*pow(z,3)*(-72*r + 44*c*pow(z,2))))/
                             (2.*(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) + 6*pow(l,4)*(175*s - 184*z)*pow(z,6) -
                                  7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) - 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) - 28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                  12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                  12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))) -
                             ((pow(l,2) - pow(z,2))*(-105*pow(l,10) + 525*pow(l,8)*pow(z,2) - 1050*pow(l,6)*pow(z,4) + 1050*pow(l,4)*pow(z,6) -
                                                         525*pow(l,2)*pow(z,8) + 105*pow(z,10))*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                                                                                                       105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                                                                       210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                                                                                       210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                                                                                       2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                                                                       7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                                                                       2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(-105*pow(l,10)*s + 12*c*pow(l,13)*r*z + 35*pow(l,8)*(15*s - 16*z)*pow(z,2) +
                                       6*pow(l,4)*(175*s - 184*z)*pow(z,6) - 7*pow(l,2)*(75*s - 16*z)*pow(z,8) + 105*s*pow(z,10) -
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 4*pow(l,11)*(18 + 27*pow(r,2)*z - 2*c*r*pow(z,3)) -
                                       28*pow(l,5)*pow(z,6)*(70 + 45*pow(r,2)*z + c*r*pow(z,3)) -
                                       12*pow(l,9)*pow(z,2)*(42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) +
                                       12*pow(l,7)*pow(z,4)*(-14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)),2));
            }
        case 'R':
            switch (label){
                case 'l' :
                    result = -((pow(l,2) - pow(z,2))*(-1050*pow(l,9) + 39*c*pow(l,12)*r - 210*l*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                          840*pow(l,3)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                                          1260*pow(l,5)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 840*pow(l,7)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                                          14*pow(l,6)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                          35*pow(l,4)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 66*pow(l,10)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                          18*pow(l,8)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
                                  105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                  12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                  12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))) +
                             ((pow(l,2) - pow(z,2))*(1050*pow(l,9)*s - 156*c*pow(l,12)*r*z - 24*pow(l,3)*(175*s - 184*z)*pow(z,6) +
                                                         14*l*(75*s - 16*z)*pow(z,8) + 280*pow(l,7)*pow(z,2)*(-15*s + 16*z) + 252*pow(l,5)*pow(z,4)*(25*s + 24*z) +
                                                         140*pow(l,4)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                                         108*pow(l,8)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                                         84*pow(l,6)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + 11*pow(l,10)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) +
                                       7*pow(l,2)*(75*s - 16*z)*pow(z,8) - 105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) +
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)),2))\
    - (l*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
          105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
          210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
          210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
          2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
          7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
          2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
      (105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
       105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
       28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)));
                case 'r' :
                    result = -((pow(l,2) - pow(z,2))*(3*c*pow(l,13) + 54*pow(l,11)*r - 54*pow(l,9)*r*(4*s - z)*z -
                                                          18*pow(l,7)*r*(8*s - 49*z)*pow(z,3) - 630*pow(l,5)*r*(4*s - 3*z)*pow(z,5) +
                                                          2*pow(l,7)*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                          7*pow(l,5)*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                          2*pow(l,9)*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
                                  105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                  12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                  12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))) +
                             ((pow(l,2) - pow(z,2))*(-12*c*pow(l,13)*z + 28*pow(l,5)*pow(z,6)*(90*r*z + c*pow(z,3)) +
                                                         12*pow(l,9)*pow(z,2)*(6*r*z + 4*c*pow(z,3)) - 12*pow(l,7)*pow(z,4)*(198*r*z + 6*c*pow(z,3)) +
                                                         pow(l,11)*(-216*r*z + 8*c*pow(z,3)))*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                                                                                                   105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                                                                   210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                                                                                   210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                                                                                   2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                                                                   7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                                                                   2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) +
                                       7*pow(l,2)*(75*s - 16*z)*pow(z,8) - 105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) +
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)),2));
                case 'z' :
                    result = -((pow(l,2) - pow(z,2))*(840*s*(2*s - z)*pow(z,7) - 105*s*pow(z,8) - 105*pow(l,2)*pow(z,6)*(-4*s + 2*z) +
                                                          210*pow(l,4)*pow(z,4)*(-3*s + 4*z) - 210*pow(l,6)*pow(z,2)*(-2*s + 6*z) + 105*pow(l,8)*(-s + 8*z) -
                                                          630*pow(l,2)*pow(z,5)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                          840*pow(l,4)*pow(z,3)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) - 420*pow(l,6)*z*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) +
                                                          2*pow(l,7)*r*pow(z,3)*(441*r + 4*c*(22*s - 7*z)*z - 14*c*pow(z,2)) -
                                                          7*pow(l,5)*r*pow(z,5)*(-270*r + 2*c*(8*s - 5*z)*z - 5*c*pow(z,2)) +
                                                          6*pow(l,7)*r*pow(z,2)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                          35*pow(l,5)*r*pow(z,4)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(6*c*z + 2*c*(-2*s + 3*z)) -
                                                          2*pow(l,9)*r*z*(-27*r + 23*c*pow(z,2) + 2*c*z*(4*s + 23*z)) - 2*pow(l,9)*r*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
                                  105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                  12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                  12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))) +
                             ((pow(l,2) - pow(z,2))*(-12*c*pow(l,13)*r + 560*pow(l,8)*pow(z,2) + 1008*pow(l,6)*pow(z,4) -
                                                         36*pow(l,4)*(175*s - 184*z)*pow(z,5) + 1104*pow(l,4)*pow(z,6) + 56*pow(l,2)*(75*s - 16*z)*pow(z,7) -
                                                         112*pow(l,2)*pow(z,8) - 1050*s*pow(z,9) + 70*pow(l,8)*z*(-15*s + 16*z) + 168*pow(l,6)*pow(z,3)*(25*s + 24*z) +
                                                         28*pow(l,5)*pow(z,6)*(45*pow(r,2) + 3*c*r*pow(z,2)) + 12*pow(l,9)*pow(z,2)*(3*pow(r,2) + 12*c*r*pow(z,2)) -
                                                         12*pow(l,7)*pow(z,4)*(99*pow(r,2) + 18*c*r*pow(z,2)) + pow(l,11)*(-108*pow(r,2) + 24*c*r*pow(z,2)) +
                                                         168*pow(l,5)*pow(z,5)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) + 24*pow(l,9)*z*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                                         48*pow(l,7)*pow(z,3)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) +
                                       7*pow(l,2)*(75*s - 16*z)*pow(z,8) - 105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) +
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)),2))\
    + (z*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
          105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
          210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
          210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
          2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
          7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
          2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
      (105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
       105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
       28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)));
                case 'c' :
                    result = -((pow(l,2) - pow(z,2))*(3*pow(l,13)*r + 4*pow(l,7)*r*(22*s - 7*z)*pow(z,5) - 7*pow(l,5)*r*(8*s - 5*z)*pow(z,7) +
                                                          12*pow(l,11)*r*z*(-2*s + 3*z) - 2*pow(l,9)*r*pow(z,3)*(4*s + 23*z)))/
                             (2.*(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
                                  105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                  12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                  12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))) +
                             ((pow(l,2) - pow(z,2))*(-12*pow(l,13)*r*z + 8*pow(l,11)*r*pow(z,3) + 48*pow(l,9)*r*pow(z,5) -
                                                         72*pow(l,7)*r*pow(z,7) + 28*pow(l,5)*r*pow(z,9))*
                              (-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                               105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                               210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                               210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                               2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                               7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                               2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) +
                                       7*pow(l,2)*(75*s - 16*z)*pow(z,8) - 105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) +
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)),2));
                case 's' :
                    result = -((pow(l,2) - pow(z,2))*(105*pow(l,8)*(4*s - z) - 24*c*pow(l,11)*r*z - 210*pow(l,6)*(8*s - 2*z)*pow(z,2) +
                                                          210*pow(l,4)*(12*s - 3*z)*pow(z,4) - 105*pow(l,2)*(16*s - 4*z)*pow(z,6) + 210*s*pow(z,8) + 105*(2*s - z)*pow(z,8) -
                                                          2*pow(l,9)*r*z*(108*r + 4*c*pow(z,2)) - 7*pow(l,5)*r*pow(z,5)*(360*r + 8*c*pow(z,2)) +
                                                          2*pow(l,7)*r*pow(z,3)*(-72*r + 44*c*pow(z,2))))/
                             (2.*(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) + 7*pow(l,2)*(75*s - 16*z)*pow(z,8) -
                                  105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) + 42*pow(l,6)*pow(z,4)*(25*s + 24*z) +
                                  28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                  12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                  12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)))) +
                             ((pow(l,2) - pow(z,2))*(105*pow(l,10) - 525*pow(l,8)*pow(z,2) + 1050*pow(l,6)*pow(z,4) - 1050*pow(l,4)*pow(z,6) +
                                                         525*pow(l,2)*pow(z,8) - 105*pow(z,10))*(-105*pow(l,10) + 3*c*pow(l,13)*r + 105*s*(2*s - z)*pow(z,8) -
                                                                                                       105*pow(l,2)*pow(z,6)*(8*pow(s,2) - 4*s*z + pow(z,2)) +
                                                                                                       210*pow(l,4)*pow(z,4)*(6*pow(s,2) - 3*s*z + 2*pow(z,2)) -
                                                                                                       210*pow(l,6)*pow(z,2)*(4*pow(s,2) - 2*s*z + 3*pow(z,2)) + 105*pow(l,8)*(2*pow(s,2) - s*z + 4*pow(z,2)) +
                                                                                                       2*pow(l,7)*r*pow(z,3)*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*pow(z,2)) -
                                                                                                       7*pow(l,5)*r*pow(z,5)*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*pow(z,2)) + 6*pow(l,11)*r*(9*r + 2*c*z*(-2*s + 3*z)) -
                                                                                                       2*pow(l,9)*r*z*(27*r*(4*s - z) + c*pow(z,2)*(4*s + 23*z))))/
                             (2.*pow(105*pow(l,10)*s - 12*c*pow(l,13)*r*z - 6*pow(l,4)*(175*s - 184*z)*pow(z,6) +
                                       7*pow(l,2)*(75*s - 16*z)*pow(z,8) - 105*s*pow(z,10) + 35*pow(l,8)*pow(z,2)*(-15*s + 16*z) +
                                       42*pow(l,6)*pow(z,4)*(25*s + 24*z) + 28*pow(l,5)*pow(z,6)*(-70 + 45*pow(r,2)*z + c*r*pow(z,3)) +
                                       12*pow(l,9)*pow(z,2)*(-42 + 3*pow(r,2)*z + 4*c*r*pow(z,3)) -
                                       12*pow(l,7)*pow(z,4)*(14 + 99*pow(r,2)*z + 6*c*r*pow(z,3)) + pow(l,11)*(72 - 108*pow(r,2)*z + 8*c*r*pow(z,3)),2));
            }
    }
    return result;
}
