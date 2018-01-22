#include "../include/shape.h"
#include <iostream>
#include <math>

double* store(string file_path);

int* shape::grid(double * starting_point, double * step_length){
    double shape_para[5]={_para_l,_para_r,_para_z,_para_c,_para_s};
    int lower_limit[5];
    for(i=0;i<5;i++){
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
    double term0,term1,term2;
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
double shape::energy(double * starting_point, double * step_length, double* storation, int* step) {
    int* location;
    location = shape::grid(starting_point, step_length);
    return storation[location[0]+step[0],location[1]+step[1],location[2]+step[2],location[3]+step[3],location[4]+step[4]];
}
double shape::energy_cubic_spline(){
    double weight ;
    return weight*
}