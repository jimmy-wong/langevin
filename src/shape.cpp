#include "../include/shape.h"
#include <boost/assert.hpp>
#include <algorithm>

using namespace std;

void shape::set_gs(int *gs){
	for(int i=0; i<5; i++){
		_gs[i] = gs[i];
	}

}
void shape::set_steps(int *steps) {
    for(int i=0; i<5; i++){
        _steps[i] = steps[i];
    }
}
void shape::grid(double *starting_point, double *step_length, int* lower_limit){
    double shape_para[5]={_para_l,_para_r,_para_z,_para_c,_para_s};
    for(int i=0;i<5;i++){
        lower_limit[i] = floor((shape_para[i]-starting_point[i])/step_length[i]);
    }
}
void shape::coefficiency() {
    // shape_para unit: Rcn, not fm
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;
    // auxiliary variable
    double term0;
    double l2,l3,l4,l5,l6,l7,l8,z2,z3,z4,z5,z6,z7,z8,r2;

    l2=l*l;  l3=l2*l;  l4=l2*l2;  l5=l3*l2;  l6=l3*l3;  l7=l3*l4;  l8=l4*l4;
    z2=z*z;  z3=z2*z;  z4=z2*z2;  z5=z3*z2;  z6=z3*z3;  z7=z3*z4;  z8=z4*z4;
    r2=r*r;

    term0 = l2-z2;
    _a0 = r*r/term0;
    _a1 = 2.*_a0*z/term0;
    _a2 = (c*r+_a0+2.*_a1*z)/term0;
    _a3 = (-7*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z -
	       165*l8*s*z2 + 40*l8*z3 - 8*c*l5*l6*r*z3 +
	       410*l6*s*z4 + 240*l6*z5 + 68*c*l4*l5*r*z5 -
	       312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
	       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
	       140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
	(l5*pow((l2 - z2),3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245 *z6));
    _a4 = (-7*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z -
	       60*l8*z2 - 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 +
	       140*l6*s*z3 + 270*l6*z4 + 50*c*l4*l5*r*z4 -
	       330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
	       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 +
	       105*l2*z8 + 35*c*l5*r*z8 + 175*s*z4*z5))/
	(l5*pow((l2 - z2),3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
}
double shape::grid_energy(double* storation, int* step) {
    int index;
    index = step[0]*_steps[1]*_steps[2]*_steps[3]*_steps[4]+
	step[1]*_steps[2]*_steps[3]*_steps[4]+
	step[2]*_steps[3]*_steps[4]+
	step[3]*_steps[4]+
	step[4];
    if (equal(step,step+5,_gs)) {
		return storation[index];
    }
    else{return storation[index]+
				50*(exp(1/pow(step[0]-_steps[0],4))+
					exp(1/pow(step[1]-_steps[1],4))+
					exp(1/pow(step[2]-_steps[2],4))+
					exp(1/pow(step[3]-_steps[3],4))+
					exp(1/pow(step[4]-_steps[4],4))+
					exp(1/pow(step[0],4))+
					exp(1/pow(step[1],4))+
					exp(1/pow(step[2],4))+
					exp(1/pow(step[3],4))+
					exp(1/pow(step[4],4))-10);}
}
double shape::AH(gsl_vector* generalized_coordinates){
    _para_l = gsl_vector_get(generalized_coordinates,0);
    _para_r = gsl_vector_get(generalized_coordinates,1);
    _para_z = gsl_vector_get(generalized_coordinates,2);
    _para_c = gsl_vector_get(generalized_coordinates,3);
    _para_s = gsl_vector_get(generalized_coordinates,4);
    coefficiency();

    double term12,term22;
    double l = _para_l, z = _para_z, s=_para_s;
    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4;
    double z2=z*2, z3=z*z2, z4=z2*z2;
    double s2=s*2, s3=s*s2, s4=s2*s2, s5=s2*s3, s6=s3*s3, s7=s3*s4;

    term12 = (2*_a0*l3)/3. - (_a1*l4)/4. + (2*_a2*l5)/15. - (_a3*l6)/12. + (2*_a4*l7)/35. + _a0*l2*s +
             (_a1*l2*s2)/2. - (_a0*s3)/3. + (_a2*l2*s3)/3. - (_a1*s4)/4. +
             (_a3*l2*s4)/4. - (_a2*s5)/5. + (_a4*l2*s5)/5. - (_a3*s6)/6. - (_a4*s7)/7. -
             (2*_a1*l3*z)/3. + (_a2*l4*z)/2. - (2*_a3*l5*z)/5. + (_a4*l6*z)/3. - _a1*l2*s*z -
             _a2*l2*s2*z + (_a1*s3*z)/3. - _a3*l2*s3*z + (_a2*s4*z)/2. - _a4*l2*s4*z +
             (3*_a3*s5*z)/5. + (2*_a4*s6*z)/3. + (2*_a2*l3*z2)/3. - (3*_a3*l4*z2)/4. +
             (4*_a4*l5*z2)/5. + _a2*l2*s*z2 + (3*_a3*l2*s2*z2)/2. -
             (_a2*s3*z2)/3. + 2*_a4*l2*s3*z2 - (3*_a3*s4*z2)/4. -
             (6*_a4*s5*z2)/5. - (2*_a3*l3*z3)/3. + _a4*l4*z3 - _a3*l2*s*z3 -
             2*_a4*l2*s2*z3 + (_a3*s3*z3)/3. + _a4*s4*z3 + (2*_a4*l3*z4)/3. +
             _a4*l2*s*z4 - (_a4*s3*z4)/3.;
    term22 = (2*_a0*l3)/3. + (_a1*l4)/4. + (2*_a2*l5)/15. + (_a3*l6)/12. + (2*_a4*l7)/35. - _a0*l2*s -
             (_a1*l2*s2)/2. + (_a0*s3)/3. - (_a2*l2*s3)/3. + (_a1*s4)/4. -
             (_a3*l2*s4)/4. + (_a2*s5)/5. - (_a4*l2*s5)/5. + (_a3*s6)/6. + (_a4*s7)/7. -
             (2*_a1*l3*z)/3. - (_a2*l4*z)/2. - (2*_a3*l5*z)/5. - (_a4*l6*z)/3. + _a1*l2*s*z +
             _a2*l2*s2*z - (_a1*s3*z)/3. + _a3*l2*s3*z - (_a2*s4*z)/2. + _a4*l2*s4*z -
             (3*_a3*s5*z)/5. - (2*_a4*s6*z)/3. + (2*_a2*l3*z2)/3. + (3*_a3*l4*z2)/4. +
             (4*_a4*l5*z2)/5. - _a2*l2*s*z2 - (3*_a3*l2*s2*z2)/2. +
             (_a2*s3*z2)/3. - 2*_a4*l2*s3*z2 + (3*_a3*s4*z2)/4. +
             (6*_a4*s5*z2)/5. - (2*_a3*l3*z3)/3. - _a4*l4*z3 + _a3*l2*s*z3 +
             2*_a4*l2*s2*z3 - (_a3*s3*z3)/3. - _a4*s4*z3 + (2*_a4*l3*z4)/3. -
             _a4*l2*s*z4 + (_a4*s3*z4)/3.;
    return _Acn*(term12/(term12+term22));
}
double shape::Rho(double x){
    double result;
    double l = _para_l;
    double z = _para_z;
    double s = _para_s;
    result = pow(_Rcn,2)*(pow(l,2) - pow(s + x,2))*(_a0 + _a1*(s + x - z) + _a2*pow(s + x - z,2) + _a3*pow(s + x - z,3) + _a4*pow(s + x - z,4));
    return result;
}
double shape::RhoDerivative(double x, const char label) {
//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    double result = 0;
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;

    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4, l8=l4*l4;
    double r2=r*2;
    double z2=z*2, z3=z*z2, z4=z2*z2, z5=z2*z3, z6=z3*z3, z7=z3*z4, z8=z4*z4;

    if(label=='x'){
	result = (l2 - pow(s + x, 2)) *
	    (_a1 + 2 * _a2 * (s + x - z) + 3 * _a3 * pow(s + x - z, 2) + 4 * _a4 * pow(s + x - z, 3)) -
	    2 * (s + x) * (_a0 + _a1 * (s + x - z) + _a2 * pow(s + x - z, 2) + _a3 * pow(s + x - z, 3) +
			   _a4 * pow(s + x - z, 4));
    }
    else if (label=='l') {
	result = (l2 - pow(s + x, 2)) * ((-8 * l * r2 * (s + x - z) * z) / pow(l2 - z2, 3) -
						(2 * l * r2) / pow(l2 - z2, 2) -
						(7 * pow(s + x - z, 4) * (-150 * l4*l5 + 39 * c * l6*l6 * r +
									  198 * l5*l5 * r2 -
									  840 * pow(l, 7) * s * z -
									  480 * pow(l, 7) * z2 -
									  132 * c * l5*l5 * r * z2 +
									  810 * l8 * r2 * z2 +
									  840 * l5 * s * z3 +
									  1620 * l5 * z4 +
									  450 * c * l8 * r * z4 -
									  2310 * pow(l, 6) * r2 * z4 +
									  840 * l3 * s * z5 -
									  1200 * l3 * z6 -
									  532 * c * pow(l, 6) * r * z6 +
									  1750 * l4 * r2 * z6 -
									  840 * l * s * z7 +
									  210 * l * z8 +
									  175 * c * l4 * r * z8)) /
						(l5 * pow(l2 - z2, 3) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) - (7 * pow(s + x - z, 3) *
								       (-150 * l4*l5 * s - 600 * l4*l5 * z +
									78 * c * l6*l6 * r * z +
									792 * l5*l5 * r2 * z -
									1320 * pow(l, 7) * s * z2 +
									320 * pow(l, 7) * z3 -
									88 * c * l5*l5 * r * z3 +
									2460 * l5 * s * z4 +
									1440 * l5 * z5 +
									612 * c * l8 * r * z5 -
									2184 * pow(l, 6) * r2 * z5 -
									360 * l3 * s * z6 -
									1440 * l3 * z7 -
									952 * c * pow(l, 6) * r * z7 +
									2800 * l4 * r2 * z7 -
									630 * l * s * z8 + 280 * l * z4*z5 +
									350 * c * l4 * r * z4*z5)) /
						(l5 * pow(l2 - z2, 3) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) + (7 * pow(s + x - z, 4) *
								       (54 * l5 -
									252 * l3 *
									z2 -
									42 * l * z4) *
								       (-15 * l5*l5 +
									3 * c * pow(l, 13) * r +
									18 * l5*l6 *
									r2 -
									105 * l8 * s * z -
									60 * l8 *
									z2 -
									12 * c * l5*l6 * r *
									z2 +
									90 * l4*l5 *
									r2 * z2 +
									140 * pow(l, 6) * s *
									z3 +
									270 * pow(l, 6) *
									z4 +
									50 * c * l4*l5 * r *
									z4 -
									330 * pow(l, 7) *
									r2 * z4 +
									210 * l4 * s *
									z5 -
									300 * l4 *
									z6 -
									76 * c * pow(l, 7) * r *
									z6 +
									350 * l5 *
									r2 * z6 -
									420 * l2 * s *
									z7 +
									105 * l2 *
									z8 +
									35 * c * l5 * r *
									z8 +
									175 * s * z4*z5)) /
						(l5 *
						 pow(l2 - z2,
						     3) *
						 pow(9 * pow(l, 6) -
						     63 * l4 *
						     z2 -
						     21 * l2 *
						     z4 -
						     245 * z6, 2)) +
						(42 * pow(s + x - z, 4) *
						 (-15 * l5*l5 + 3 * c * pow(l, 13) * r +
						  18 * l5*l6 * r2 - 105 * l8 * s * z -
						  60 * l8 * z2 -
						  12 * c * l5*l6 * r * z2 +
						  90 * l4*l5 * r2 * z2 +
						  140 * pow(l, 6) * s * z3 +
						  270 * pow(l, 6) * z4 + 50 * c * l4*l5 * r * z4 -
						  330 * pow(l, 7) * r2 * z4 +
						  210 * l4 * s * z5 - 300 * l4 * z6 -
						  76 * c * pow(l, 7) * r * z6 +
						  350 * l5 * r2 * z6 -
						  420 * l2 * s * z7 + 105 * l2 * z8 +
						  35 * c * l5 * r * z8 + 175 * s * z4*z5)) /
						(l4 * pow(l2 - z2, 4) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) + (35 * pow(s + x - z, 4) *
								       (-15 * l5*l5 + 3 * c * pow(l, 13) * r +
									18 * l5*l6 * r2 -
									105 * l8 * s * z -
									60 * l8 * z2 -
									12 * c * l5*l6 * r * z2 +
									90 * l4*l5 * r2 * z2 +
									140 * pow(l, 6) * s * z3 +
									270 * pow(l, 6) * z4 +
									50 * c * l4*l5 * r * z4 -
									330 * pow(l, 7) * r2 * z4 +
									210 * l4 * s * z5 -
									300 * l4 * z6 -
									76 * c * pow(l, 7) * r * z6 +
									350 * l5 * r2 * z6 -
									420 * l2 * s * z7 +
									105 * l2 * z8 +
									35 * c * l5 * r * z8 +
									175 * s * z4*z5)) /
						(pow(l, 6) * pow(l2 - z2, 3) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) + (7 * pow(s + x - z, 3) *
								       (54 * l5 -
									252 * l3 *
									z2 -
									42 * l * z4) *
								       (-15 * l5*l5 * s -
									60 * l5*l5 * z +
									6 * c * pow(l, 13) * r *
									z + 72 * l5*l6 *
									r2 * z -
									165 * l8 * s *
									z2 +
									40 * l8 *
									z3 -
									8 * c * l5*l6 * r *
									z3 +
									410 * pow(l, 6) * s *
									z4 +
									240 * pow(l, 6) *
									z5 +
									68 * c * l4*l5 * r *
									z5 -
									312 * pow(l, 7) *
									r2 * z5 -
									90 * l4 * s *
									z6 -
									360 * l4 *
									z7 -
									136 * c * pow(l, 7) * r *
									z7 +
									560 * l5 *
									r2 * z7 -
									315 * l2 * s *
									z8 +
									140 * l2 *
									z4*z5 +
									70 * c * l5 * r *
									z4*z5 +
									175 * s * z5*z5)) /
						(l5 *
						 pow(l2 - z2,
						     3) *
						 pow(9 * pow(l, 6) -
						     63 * l4 *
						     z2 -
						     21 * l2 *
						     z4 -
						     245 * z6, 2)) +
						(42 * pow(s + x - z, 3) *
						 (-15 * l5*l5 * s - 60 * l5*l5 * z +
						  6 * c * pow(l, 13) * r * z + 72 * l5*l6 * r2 * z -
						  165 * l8 * s * z2 +
						  40 * l8 * z3 - 8 * c * l5*l6 * r * z3 +
						  410 * pow(l, 6) * s * z4 + 240 * pow(l, 6) * z5 +
						  68 * c * l4*l5 * r * z5 -
						  312 * pow(l, 7) * r2 * z5 -
						  90 * l4 * s * z6 -
						  360 * l4 * z7 - 136 * c * pow(l, 7) * r * z7 +
						  560 * l5 * r2 * z7 -
						  315 * l2 * s * z8 + 140 * l2 * z4*z5 +
						  70 * c * l5 * r * z4*z5 + 175 * s * z5*z5)) /
						(l4 * pow(l2 - z2, 4) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) + (35 * pow(s + x - z, 3) *
								       (-15 * l5*l5 * s - 60 * l5*l5 * z +
									6 * c * pow(l, 13) * r * z +
									72 * l5*l6 * r2 * z -
									165 * l8 * s * z2 +
									40 * l8 * z3 -
									8 * c * l5*l6 * r * z3 +
									410 * pow(l, 6) * s * z4 +
									240 * pow(l, 6) * z5 +
									68 * c * l4*l5 * r * z5 -
									312 * pow(l, 7) * r2 * z5 -
									90 * l4 * s * z6 -
									360 * l4 * z7 -
									136 * c * pow(l, 7) * r * z7 +
									560 * l5 * r2 * z7 -
									315 * l2 * s * z8 +
									140 * l2 * z4*z5 +
									70 * c * l5 * r * z4*z5 +
									175 * s * z5*z5)) /
						(pow(l, 6) * pow(l2 - z2, 3) *
						 (9 * pow(l, 6) - 63 * l4 * z2 -
						  21 * l2 * z4 -
						  245 * z6)) + (pow(s + x - z, 2) *
								       ((-16 * l * r2 *
									 z2) /
									pow(l2 - z2,
									    3) -
									(2 * l * r2) /
									pow(l2 - z2,
									    2))) /
						(l2 - z2) -
						(2 * l * pow(s + x - z, 2) *
						 (c * r + (4 * r2 * z2) / pow(l2 - z2, 2) +
						  r2 / (l2 - z2))) /
						pow(l2 - z2, 2)) +
	    2 * l * ((2 * r2 * (s + x - z) * z) / pow(l2 - z2, 2) +
		     r2 / (l2 - z2) - (7 * pow(s + x - z, 4) *
							    (-15 * l5*l5 + 3 * c * pow(l, 13) * r +
							     18 * l5*l6 * r2 -
							     105 * l8 * s * z - 60 * l8 * z2 -
							     12 * c * l5*l6 * r * z2 +
							     90 * l4*l5 * r2 * z2 +
							     140 * pow(l, 6) * s * z3 +
							     270 * pow(l, 6) * z4 +
							     50 * c * l4*l5 * r * z4 -
							     330 * pow(l, 7) * r2 * z4 +
							     210 * l4 * s * z5 -
							     300 * l4 * z6 -
							     76 * c * pow(l, 7) * r * z6 +
							     350 * l5 * r2 * z6 -
							     420 * l2 * s * z7 +
							     105 * l2 * z8 +
							     35 * c * l5 * r * z8 +
							     175 * s * z4*z5)) /
		     (l5 * pow(l2 - z2, 3) *
		      (9 * pow(l, 6) - 63 * l4 * z2 -
		       21 * l2 * z4 -
		       245 * z6)) - (7 * pow(s + x - z, 3) *
					    (-15 * l5*l5 * s -
					     60 * l5*l5 * z +
					     6 * c * pow(l, 13) * r * z +
					     72 * l5*l6 * r2 *
					     z - 165 * l8 * s *
					     z2 +
					     40 * l8 * z3 -
					     8 * c * l5*l6 * r *
					     z3 +
					     410 * pow(l, 6) * s *
					     z4 +
					     240 * pow(l, 6) * z5 +
					     68 * c * l4*l5 * r *
					     z5 -
					     312 * pow(l, 7) * r2 *
					     z5 -
					     90 * l4 * s * z6 -
					     360 * l4 * z7 -
					     136 * c * pow(l, 7) * r *
					     z7 +
					     560 * l5 * r2 *
					     z7 -
					     315 * l2 * s *
					     z8 +
					     140 * l2 * z4*z5 +
					     70 * c * l5 * r *
					     z4*z5 +
					     175 * s * z5*z5)) /
		     (l5 *
		      pow(l2 - z2, 3) *
		      (9 * pow(l, 6) -
		       63 * l4 * z2 -
		       21 * l2 * z4 -
		       245 * z6)) +
		     (pow(s + x - z, 2) * (c * r + (4 * r2 * z2) / pow(l2 - z2, 2) +
					   r2 / (l2 - z2))) / (l2 - z2));
    }
    else if(label=='r'){
	result = (l2 - pow(s + x,2))*((4*r*(s + x - z)*z)/pow(l2 - z2,2) + (2*r)/(l2 - z2) -
					    (7*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								 50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
								 35*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    - (7*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								   624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (pow(s + x - z,2)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/
					    (l2 - z2));
    }
    else if(label=='z'){
	result = (l2 - pow(s + x,2))*((8*r2*(s + x - z)*z2)/pow(l2 - z2,3) +
					    (2*r2*(s + x - z))/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z +
								 420*l6*s*z2 + 1080*l6*z3 + 200*c*l4*l5*r*z3 -
								 1320*l7*r2*z3 + 1050*l4*s*z4 - 1800*l4*z5 -
								 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								 840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					       (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
						12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 +
						270*l6*z4 + 50*c*l4*l5*r*z4 - 330*l7*r2*z4 +
						210*l4*s*z5 - 300*l4*z6 - 76*c*l7*r*z6 +
						350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
						35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 -
										     245*z6,2)) - (42*pow(s + x - z,4)*z*
													 (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
													  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 +
													  270*l6*z4 + 50*c*l4*l5*r*z4 - 330*l7*r2*z4 +
													  210*l4*s*z5 - 300*l4*z6 - 76*c*l7*r*z6 +
													  350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
													  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z -
								    60*l8*z2 - 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 +
								    270*l6*z4 + 50*c*l4*l5*r*z4 - 330*l7*r2*z4 +
								    210*l4*s*z5 - 300*l4*z6 - 76*c*l7*r*z6 +
								    350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								    35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    - (7*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z +
								   120*l8*z2 - 24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 +
								   340*c*l4*l5*r*z4 - 1560*l7*r2*z4 - 540*l4*s*z5 -
								   2520*l4*z6 - 952*c*l7*r*z6 + 3920*l5*r2*z6 -
								   2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 + 1750*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					       (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
						40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
						68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
						136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
						140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 -
										     245*z6,2)) - (42*pow(s + x - z,3)*z*
													 (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
													  40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
													  68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
													  136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
													  140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z -
								    165*l8*s*z2 + 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 +
								    240*l6*z5 + 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 -
								    360*l4*z7 - 136*c*l7*r*z7 + 560*l5*r2*z7 -
								    315*l2*s*z8 + 140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    + (pow(s + x - z,2)*((16*r2*z3)/pow(l2 - z2,3) +
								 (10*r2*z)/pow(l2 - z2,2)))/(l2 - z2) +
					    (2*pow(s + x - z,2)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    pow(l2 - z2,2) - (2*(s + x - z)*
									  (c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    (l2 - z2));
    }
    else if (label=='c'){
	result = (l2 - pow(s + x,2))*((r*pow(s + x - z,2))/(l2 - z2) -
					    (7*pow(s + x - z,4)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								 35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))
					    - (7*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 -
								   136*l7*r*z7 + 70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if(label=='s'){
	result = (l2 - pow(s + x,2))*((2*r2*z)/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(-105*l8*z + 140*l6*z3 +
								 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) - (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 -
																						    105*l8*s*z - 60*l8*z2 - 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 +
																						    140*l6*s*z3 + 270*l6*z4 + 50*c*l4*l5*r*z4 -
																						    330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
																						    76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 +
																						    105*l2*z8 + 35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
										  245*z6)) - (7*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 -
															90*l4*z6 - 315*l2*z8 + 175*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
										  245*z6)) - (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z +
															 72*l5*l6*r2*z - 165*l8*s*z2 + 40*l8*z3 - 8*c*l5*l6*r*z3 +
															 410*l6*s*z4 + 240*l6*z5 + 68*c*l4*l5*r*z5 -
															 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
															 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
															 140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
										  245*z6)) + (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
														   r2/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((2*r2*(s + x - z)*z)/pow(l2 - z2,2) + r2/(l2 - z2) -
		       (7*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z -
					    60*l8*z2 - 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 +
					    140*l6*s*z3 + 270*l6*z4 + 50*c*l4*l5*r*z4 -
					    330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					    76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 +
					    105*l2*z8 + 35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
							     245*z6)) - (7*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z +
												   72*l5*l6*r2*z - 165*l8*s*z2 + 40*l8*z3 - 8*c*l5*l6*r*z3 +
												   410*l6*s*z4 + 240*l6*z5 + 68*c*l4*l5*r*z5 -
												   312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
												   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
												   140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
							     245*z6)) + (pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
												 r2/(l2 - z2)))/(l2 - z2));
    }
    else{
        BOOST_ASSERT_MSG(false, "RhoDerivative!");
    }
//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    std::cout << label << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;
    result = result*_Rcn;
    return result;
}
double shape::RhoDDerivative(double x, char label_i, char label_j) {
//	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	double result;
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;

    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4, l8=l4*l4;
    double r2=r*2;
    double z2=z*2, z3=z*z2, z4=z2*z2, z5=z2*z3, z6=z3*z3, z7=z3*z4, z8=z4*z4;

    if (label_i=='x'&&label_j=='l'){
	result = (l2 - pow(s + x,2))*((-8*l*r2*z)/pow(l2 - z2,3) -
					    (28*pow(s + x - z,3)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								  132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								  450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								  532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
								  320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
								  612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
								  952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (168*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (140*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (126*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (105*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
					    (l2 - z2) - (4*l*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
										      r2/(l2 - z2)))/pow(l2 - z2,2)) -
	    2*(s + x)*((-8*l*r2*(s + x - z)*z)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2) -
		       (7*pow(s + x - z,4)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
					    132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
					    450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
					    532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
		       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
					    320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
					    612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
					    952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
		       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
			(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
			 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
			 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
			 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
			 35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
		       (42*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (35*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
			(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
			 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
			 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
			 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
			 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
		       (42*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (35*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
		       (l2 - z2) - (2*l*pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
								      r2/(l2 - z2)))/pow(l2 - z2,2)) +
	    2*l*((2*r2*z)/pow(l2 - z2,2) - (28*pow(s + x - z,3)*
							      (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
							       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
							       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
							       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
							       35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='x' && label_j=='r'){
	result = (l2 - pow(s + x,2))*((4*r*z)/pow(l2 - z2,2) -
					    (28*pow(s + x - z,3)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								  50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
								  35*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								  624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((4*r*(s + x - z)*z)/pow(l2 - z2,2) + (2*r)/(l2 - z2) -
		       (7*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
					    50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
					    35*c*l5*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
					    624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='x'&& label_j=='z'){
	result = (l2 - pow(s + x,2))*((8*r2*z2)/pow(l2 - z2,3) + (2*r2)/pow(l2 - z2,2) -
					    (28*pow(s + x - z,3)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								  1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								  1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								  840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (168*pow(s + x - z,3)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								     35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (84*pow(s + x - z,2)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								  24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								  1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								  3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								  1750*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (126*pow(s + x - z,2)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (42*(s + x - z)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
							     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
							     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
							     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
							     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
					    (l2 - z2) + (4*(s + x - z)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
										      r2/(l2 - z2)))/pow(l2 - z2,2) -
					    (2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((8*r2*(s + x - z)*z2)/pow(l2 - z2,3) + (2*r2*(s + x - z))/pow(l2 - z2,2) -
		       (7*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
					    1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
					    1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
					    840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
			(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
			 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
			 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
			 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
			 35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		       (42*pow(s + x - z,4)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					       35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
					    24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
					    1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
					    3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
					    1750*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
			(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
			 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
			 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
			 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
			 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		       (42*pow(s + x - z,3)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
		       (l2 - z2) + (2*pow(s + x - z,2)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
								      r2/(l2 - z2)))/pow(l2 - z2,2) -
		       (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='x'&&label_j=='c'){
	result = (l2 - pow(s + x,2))*((2*r*(s + x - z))/(l2 - z2) -
					    (28*pow(s + x - z,3)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								  35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								  70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))) -
	    2*(s + x)*((r*pow(s + x - z,2))/(l2 - z2) -
		       (7*pow(s + x - z,4)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
					    35*l5*r*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
					    70*l5*r*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if (label_i=='x'&& label_j=='s'){
	result = (l2 - pow(s + x,2))*((-28*pow(s + x - z,3)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 -
								   420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (84*pow(s + x - z,2)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								  315*l2*z8 + 175*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*(s + x - z)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
							     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
							     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
							     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
							     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((2*r2*z)/pow(l2 - z2,2) -
		       (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2))
	    - 2*(s + x)*((2*r2*z)/pow(l2 - z2,2) -
			 (7*pow(s + x - z,4)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
			 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
			 (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					       35*c*l5*r*z8 + 175*s*z4*z5))/
			 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
			 (7*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
					      315*l2*z8 + 175*z5*z5))/
			 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
			 (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
			 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
			 (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2))
	    - 2*((2*r2*(s + x - z)*z)/pow(l2 - z2,2) + r2/(l2 - z2) -
		 (7*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
				      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
				      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
				      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
				      35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (7*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
		 (l2 - z2));
    }
    else if (label_i=='l'&&label_j=='l'){
	result = (l2 - pow(s + x,2))*((48*l2*r2*(s + x - z)*z)/pow(l2 - z2,4) +
					    (8*l2*r2)/pow(l2 - z2,3) - (8*r2*(s + x - z)*z)/pow(l2 - z2,3) -
					    (2*r2)/pow(l2 - z2,2) - (7*pow(s + x - z,4)*
										       (-1350*l8 + 468*c*l5*l6*r + 1980*l4*l5*r2 - 5880*l6*s*z - 3360*l6*z2 -
											1320*c*l4*l5*r*z2 + 6480*l7*r2*z2 + 4200*l4*s*z3 + 8100*l4*z4 +
											3600*c*l7*r*z4 - 13860*l5*r2*z4 + 2520*l2*s*z5 - 3600*l2*z6 -
											3192*c*l5*r*z6 + 7000*l3*r2*z6 - 840*s*z7 + 210*z8 + 700*c*l3*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (14*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
					      132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
					      450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
					      532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (84*pow(s + x - z,4)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								  132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								  450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								  532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (70*pow(s + x - z,4)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								  132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								  450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								  532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(-1350*l8*s - 5400*l8*z + 936*c*l5*l6*r*z + 7920*l4*l5*r2*z -
								 9240*l6*s*z2 + 2240*l6*z3 - 880*c*l4*l5*r*z3 + 12300*l4*s*z4 +
								 7200*l4*z5 + 4896*c*l7*r*z5 - 13104*l5*r2*z5 - 1080*l2*s*z6 -
								 4320*l2*z7 - 5712*c*l5*r*z7 + 11200*l3*r2*z7 - 630*s*z8 + 280*z4*z5 +
								 1400*c*l3*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (14*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
					      320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
					      612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
					      952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (84*pow(s + x - z,3)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
								  320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
								  612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
								  952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (70*pow(s + x - z,3)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
								  320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
								  612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
								  952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,4)*pow(54*l5 - 252*l3*z2 - 42*l*z4,2)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,4)*(270*l4 - 756*l2*z2 - 42*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (84*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (70*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (336*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l3*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (378*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (210*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l7*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,3)*pow(54*l5 - 252*l3*z2 - 42*l*z4,2)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,3)*(270*l4 - 756*l2*z2 - 42*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (84*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (70*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (336*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l3*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (378*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (210*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l7*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((96*l2*r2*z2)/pow(l2 - z2,4) +
							       (8*l2*r2)/pow(l2 - z2,3) - (16*r2*z2)/pow(l2 - z2,3) -
							       (2*r2)/pow(l2 - z2,2)))/(l2 - z2) -
					    (4*l*pow(s + x - z,2)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
					    pow(l2 - z2,2) + (8*l2*pow(s + x - z,2)*
									  (c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/pow(l2 - z2,3) -
					    (2*pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    pow(l2 - z2,2)) + 4*l*((-8*l*r2*(s + x - z)*z)/pow(l2 - z2,3) -
									       (2*l*r2)/pow(l2 - z2,2) - (7*pow(s + x - z,4)*
															    (-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
															     132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
															     450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
															     532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
									       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
									       (7*pow(s + x - z,3)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
												    320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
												    612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
												    952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
									       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
									       (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
										(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
										 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
										 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
										 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
										 35*c*l5*r*z8 + 175*s*z4*z5))/
									       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
									       (42*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
												     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
												     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
												     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
												     35*c*l5*r*z8 + 175*s*z4*z5))/
									       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
									       (35*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
												     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
												     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
												     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
												     35*c*l5*r*z8 + 175*s*z4*z5))/
									       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
									       (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
										(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
										 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
										 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
										 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
										 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
									       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
									       (42*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
												     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
												     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
												     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
												     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
									       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
									       (35*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
												     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
												     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
												     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
												     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
									       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
									       (pow(s + x - z,2)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
									       (l2 - z2) - (2*l*pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
															      r2/(l2 - z2)))/pow(l2 - z2,2)) +
	    2*((2*r2*(s + x - z)*z)/pow(l2 - z2,2) + r2/(l2 - z2) -
	       (7*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
				    12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
				    50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
				    76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
				    35*c*l5*r*z8 + 175*s*z4*z5))/
	       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	       (7*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				    40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				    68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				    136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				    70*c*l5*r*z4*z5 + 175*s*z5*z5))/
	       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	       (pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
	       (l2 - z2));
    }
    else if (label_i=='l'&&label_j=='r'){
	result = (l2 - pow(s + x,2))*((-16*l*r*(s + x - z)*z)/pow(l2 - z2,3) - (4*l*r)/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(39*c*l6*l6 + 396*l5*l5*r - 132*c*l5*l5*z2 + 1620*l8*r*z2 +
								 450*c*l8*z4 - 4620*l6*r*z4 - 532*c*l6*z6 + 3500*l4*r*z6 +
								 175*c*l4*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 + 50*c*l4*l5*z4 -
					      660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 + 35*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								  50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
								  35*c*l5*z8))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								  50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
								  35*c*l5*z8))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(78*c*l6*l6*z + 1584*l5*l5*r*z - 88*c*l5*l5*z3 + 612*c*l8*z5 -
								 4368*l6*r*z5 - 952*c*l6*z7 + 5600*l4*r*z7 + 350*c*l4*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 - 624*l7*r*z5 -
					      136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								  624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								  624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((-32*l*r*z2)/pow(l2 - z2,3) - (4*l*r)/pow(l2 - z2,2)))/(l2 - z2) -
					    (2*l*pow(s + x - z,2)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/pow(l2 - z2,2))
	    + 2*l*((4*r*(s + x - z)*z)/pow(l2 - z2,2) + (2*r)/(l2 - z2) -
		   (7*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
					50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
					35*c*l5*z8))/
		   (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		   (7*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
					624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
		   (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		   (pow(s + x - z,2)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/(l2 - z2));
    }
    else if(label_i=='l'&&label_j=='z'){
	result = (l2 - pow(s + x,2))*((-48*l*r2*(s + x - z)*z2)/pow(l2 - z2,4) -
					    (8*l*r2*(s + x - z))/pow(l2 - z2,3) -
					    (7*pow(s + x - z,4)*(-840*l7*s - 960*l7*z - 264*c*l5*l5*r*z + 1620*l8*r2*z + 2520*l5*s*z2 +
								 6480*l5*z3 + 1800*c*l8*r*z3 - 9240*l6*r2*z3 + 4200*l3*s*z4 -
								 7200*l3*z5 - 3192*c*l6*r*z5 + 10500*l4*r2*z5 - 5880*l*s*z6 +
								 1680*l*z7 + 1400*c*l4*r*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
					      132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
					      450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
					      532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,4)*z*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								    132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								    450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								    532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								  132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								  450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								  532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(-600*l4*l5 + 78*c*l6*l6*r + 792*l5*l5*r2 - 2640*l7*s*z + 960*l7*z2 -
								 264*c*l5*l5*r*z2 + 9840*l5*s*z3 + 7200*l5*z4 + 3060*c*l8*r*z4 -
								 10920*l6*r2*z4 - 2160*l3*s*z5 - 10080*l3*z6 - 6664*c*l6*r*z6 +
								 19600*l4*r2*z6 - 5040*l*s*z7 + 2520*l*z8 + 3150*c*l4*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
					      1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
					      1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
					      840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								  1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								  1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								  840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								  1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								  1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								  840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
					      320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
					      612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
					      952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,3)*z*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z -
								    1320*l7*s*z2 + 320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 +
								    1440*l5*z5 + 612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 -
								    1440*l3*z7 - 952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 +
								    350*c*l4*r*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
								  320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
								  612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
								  952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,4)*(-504*l3*z - 168*l*z3)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,4)*z*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (28*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (35*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (336*pow(s + x - z,4)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								     35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (168*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (210*pow(s + x - z,4)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								     35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (140*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
					      24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
					      1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
					      3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
					      1750*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								  24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								  1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								  3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								  1750*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								  24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								  1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								  3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								  1750*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,3)*(-504*l3*z - 168*l*z3)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,3)*z*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (21*pow(s + x - z,2)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (35*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (336*pow(s + x - z,3)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (126*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (210*pow(s + x - z,3)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (105*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((-96*l*r2*z3)/pow(l2 - z2,4) - (40*l*r2*z)/pow(l2 - z2,3)))/
					    (l2 - z2) + (2*pow(s + x - z,2)*z*((-16*l*r2*z2)/pow(l2 - z2,3) -
											   (2*l*r2)/pow(l2 - z2,2)))/pow(l2 - z2,2) -
					    (2*(s + x - z)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
					    (l2 - z2) - (2*l*pow(s + x - z,2)*((16*r2*z3)/pow(l2 - z2,3) +
											   (10*r2*z)/pow(l2 - z2,2)))/pow(l2 - z2,2) -
					    (8*l*pow(s + x - z,2)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    pow(l2 - z2,3) + (4*l*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
											   r2/(l2 - z2)))/pow(l2 - z2,2)) +
	    2*l*((8*r2*(s + x - z)*z2)/pow(l2 - z2,3) + (2*r2*(s + x - z))/pow(l2 - z2,2) -
		 (7*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
				      1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
				      1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
				      840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
		  (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
		   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
		   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
		   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
		   35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		 (42*pow(s + x - z,4)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					 35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
				       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
				       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
				       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
				       35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (7*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
				      24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
				      1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
				      3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
				      1750*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
		  (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
		   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
		   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
		   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
		   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		 (42*pow(s + x - z,3)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (pow(s + x - z,2)*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
		 (l2 - z2) + (2*pow(s + x - z,2)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
								r2/(l2 - z2)))/pow(l2 - z2,2) -
		 (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='l'&&label_j=='c'){
	result = (l2 - pow(s + x,2))*((-2*l*r*pow(s + x - z,2))/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(39*l6*l6*r - 132*l5*l5*r*z2 + 450*l8*r*z4 - 532*l6*r*z6 +
								 175*l4*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 + 35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,4)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								  35*l5*r*z8))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,4)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								  35*l5*r*z8))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(78*l6*l6*r*z - 88*l5*l5*r*z3 + 612*l8*r*z5 - 952*l6*r*z7 +
								 350*l4*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 + 70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								  70*l5*r*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								  70*l5*r*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))) +
	    2*l*((r*pow(s + x - z,2))/(l2 - z2) - (7*pow(s + x - z,4)*
							       (3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 + 35*l5*r*z8))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (7*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
				      70*l5*r*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if (label_i=='l'&&label_j=='s'){
	result = (l2 - pow(s + x,2))*((-8*l*r2*z)/pow(l2 - z2,3) -
					    (7*pow(s + x - z,4)*(-840*l7*z + 840*l5*z3 + 840*l3*z5 - 840*l*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(-150*l4*l5 - 1320*l7*z2 + 2460*l5*z4 - 360*l3*z6 - 630*l*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (28*pow(s + x - z,3)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
								  132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
								  450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
								  532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,4)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,4)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
								  320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
								  612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
								  952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
					    /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (168*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (140*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								   12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								   50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								   76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								   35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 - 315*l2*z8 +
					      175*z5*z5))/(l5*pow(l2 - z2,3)*
							       pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (42*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								  315*l2*z8 + 175*z5*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (35*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								  315*l2*z8 + 175*z5*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(54*l5 - 252*l3*z2 - 42*l*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (126*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (105*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								   40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								   68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								   136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								   70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
					    (l2 - z2) - (4*l*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
										      r2/(l2 - z2)))/pow(l2 - z2,2)) -
	    2*(s + x)*((-8*l*r2*(s + x - z)*z)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2) -
		       (7*pow(s + x - z,4)*(-150*l4*l5 + 39*c*l6*l6*r + 198*l5*l5*r2 - 840*l7*s*z - 480*l7*z2 -
					    132*c*l5*l5*r*z2 + 810*l8*r2*z2 + 840*l5*s*z3 + 1620*l5*z4 +
					    450*c*l8*r*z4 - 2310*l6*r2*z4 + 840*l3*s*z5 - 1200*l3*z6 -
					    532*c*l6*r*z6 + 1750*l4*r2*z6 - 840*l*s*z7 + 210*l*z8 + 175*c*l4*r*z8))
		       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(-150*l4*l5*s - 600*l4*l5*z + 78*c*l6*l6*r*z + 792*l5*l5*r2*z - 1320*l7*s*z2 +
					    320*l7*z3 - 88*c*l5*l5*r*z3 + 2460*l5*s*z4 + 1440*l5*z5 +
					    612*c*l8*r*z5 - 2184*l6*r2*z5 - 360*l3*s*z6 - 1440*l3*z7 -
					    952*c*l6*r*z7 + 2800*l4*r2*z7 - 630*l*s*z8 + 280*l*z4*z5 + 350*c*l4*r*z4*z5))
		       /(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,4)*(54*l5 - 252*l3*z2 - 42*l*z4)*
			(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
			 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
			 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
			 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
			 35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
		       (42*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (35*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,3)*(54*l5 - 252*l3*z2 - 42*l*z4)*
			(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
			 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
			 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
			 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
			 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
		       (42*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (35*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*((-16*l*r2*z2)/pow(l2 - z2,3) - (2*l*r2)/pow(l2 - z2,2)))/
		       (l2 - z2) - (2*l*pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
								      r2/(l2 - z2)))/pow(l2 - z2,2)) +
	    2*l*((2*r2*z)/pow(l2 - z2,2) - (7*pow(s + x - z,4)*
							      (-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
				       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
				       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
				       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
				       35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (7*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
				      315*l2*z8 + 175*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='r'&&label_j=='r'){
	result = (l2 - pow(s + x,2))*((4*(s + x - z)*z)/pow(l2 - z2,2) + 2/(l2 - z2) -
					    (7*pow(s + x - z,4)*(36*l5*l6 + 180*l4*l5*z2 - 660*l7*z4 + 700*l5*z6))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(144*l5*l6*z - 624*l7*z5 + 1120*l5*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((8*z2)/pow(l2 - z2,2) + 2/(l2 - z2)))/(l2 - z2));
    }
    else if(label_i=='r'&&label_j=='z'){
	result = (l2 - pow(s + x,2))*((16*r*(s + x - z)*z2)/pow(l2 - z2,3) + (4*r*(s + x - z))/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(-24*c*l5*l6*z + 360*l4*l5*r*z + 200*c*l4*l5*z3 - 2640*l7*r*z3 -
								 456*c*l7*z5 + 4200*l5*r*z5 + 280*c*l5*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 + 50*c*l4*l5*z4 -
					      660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 + 35*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,4)*z*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								    50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 + 35*c*l5*z8
						))/(l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								  50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 + 35*c*l5*z8
						))/(l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(6*c*l6*l7 + 144*l5*l6*r - 24*c*l5*l6*z2 + 340*c*l4*l5*z4 -
								 3120*l7*r*z4 - 952*c*l7*z6 + 7840*l5*r*z6 + 630*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 - 624*l7*r*z5 -
					      136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,3)*z*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								    624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								  624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((32*r*z3)/pow(l2 - z2,3) + (20*r*z)/pow(l2 - z2,2)))/(l2 - z2) +
					    (2*pow(s + x - z,2)*z*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/
					    pow(l2 - z2,2) - (2*(s + x - z)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/
					    (l2 - z2));
    }
    else if (label_i=='r'&&label_j=='c'){
	result = (l2 - pow(s + x,2))*(pow(s + x - z,2)/(l2 - z2) -
					    (7*pow(s + x - z,4)*(3*l6*l7 - 12*l5*l6*z2 + 50*l4*l5*z4 - 76*l7*z6 + 35*l5*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(6*l6*l7*z - 8*l5*l6*z3 + 68*l4*l5*z5 - 136*l7*z7 +
								 70*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if (label_i=='r'&&label_j=='s'){
	result = (l2 - pow(s + x,2))*((4*r*z)/pow(l2 - z2,2) -
					    (28*pow(s + x - z,3)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
								  50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
								  35*c*l5*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
								  624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((4*r*(s + x - z)*z)/pow(l2 - z2,2) + (2*r)/(l2 - z2) -
		       (7*pow(s + x - z,4)*(3*c*l6*l7 + 36*l5*l6*r - 12*c*l5*l6*z2 + 180*l4*l5*r*z2 +
					    50*c*l4*l5*z4 - 660*l7*r*z4 - 76*c*l7*z6 + 700*l5*r*z6 +
					    35*c*l5*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(6*c*l6*l7*z + 144*l5*l6*r*z - 8*c*l5*l6*z3 + 68*c*l4*l5*z5 -
					    624*l7*r*z5 - 136*c*l7*z7 + 1120*l5*r*z7 + 70*c*l5*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*(c + (8*r*z2)/pow(l2 - z2,2) + (2*r)/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='z'&&label_j=='z'){
	result = (l2 - pow(s + x,2))*((48*r2*(s + x - z)*z3)/pow(l2 - z2,4) +
					    (24*r2*(s + x - z)*z)/pow(l2 - z2,3) - (8*r2*z2)/pow(l2 - z2,3) -
					    (2*r2)/pow(l2 - z2,2) - (7*pow(s + x - z,4)*
										       (-120*l8 - 24*c*l5*l6*r + 180*l4*l5*r2 + 840*l6*s*z + 3240*l6*z2 +
											600*c*l4*l5*r*z2 - 3960*l7*r2*z2 + 4200*l4*s*z3 - 9000*l4*z4 -
											2280*c*l7*r*z4 + 10500*l5*r2*z4 - 17640*l2*s*z5 + 5880*l2*z6 +
											1960*c*l5*r*z6 + 12600*s*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (14*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
					      1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
					      1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
					      840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (84*pow(s + x - z,4)*z*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								    1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								    1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								    840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (56*pow(s + x - z,3)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								  1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								  1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								  840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(-330*l8*s + 240*l8*z - 48*c*l5*l6*r*z + 4920*l6*s*z2 + 4800*l6*z3 +
								 1360*c*l4*l5*r*z3 - 6240*l7*r2*z3 - 2700*l4*s*z4 - 15120*l4*z5 -
								 5712*c*l7*r*z5 + 23520*l5*r2*z5 - 17640*l2*s*z6 + 10080*l2*z7 +
								 5040*c*l5*r*z7 + 15750*s*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,4)*pow(-126*l4*z - 84*l2*z3 - 1470*z5,2)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,4)*(-126*l4 - 252*l2*z2 - 7350*z4)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (84*pow(s + x - z,4)*z*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (56*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (336*pow(s + x - z,4)*z2*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z -
									    60*l8*z2 - 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 +
									    270*l6*z4 + 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 -
									    300*l4*z6 - 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 +
									    105*l2*z8 + 35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (336*pow(s + x - z,3)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								     35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (84*pow(s + x - z,2)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (14*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
					      24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
					      1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
					      3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
					      1750*s*z4*z5))/(l5*pow(l2 - z2,3)*
								 pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (84*pow(s + x - z,3)*z*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								    24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								    1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								    3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								    1750*s*z4*z5))/(l5*pow(l2 - z2,4)*
										       (9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (42*pow(s + x - z,2)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								  24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								  1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								  3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								  1750*s*z4*z5))/(l5*pow(l2 - z2,3)*
										     (9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (14*pow(s + x - z,3)*pow(-126*l4*z - 84*l2*z3 - 1470*z5,2)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) +
					    (7*pow(s + x - z,3)*(-126*l4 - 252*l2*z2 - 7350*z4)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
					    (84*pow(s + x - z,3)*z*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,2)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (336*pow(s + x - z,3)*z2*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z -
									    165*l8*s*z2 + 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 +
									    240*l6*z5 + 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 -
									    360*l4*z7 - 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 +
									    140*l2*z4*z5 + 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								  40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								  68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								  136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								  70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (252*pow(s + x - z,2)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*(s + x - z)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
							     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
							     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
							     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
							     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (pow(s + x - z,2)*((96*r2*z4)/pow(l2 - z2,4) + (88*r2*z2)/pow(l2 - z2,3) +
							       (10*r2)/pow(l2 - z2,2)))/(l2 - z2) +
					    (4*pow(s + x - z,2)*z*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
					    pow(l2 - z2,2) - (4*(s + x - z)*((16*r2*z3)/pow(l2 - z2,3) +
											 (10*r2*z)/pow(l2 - z2,2)))/(l2 - z2) +
					    (8*pow(s + x - z,2)*z2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    pow(l2 - z2,3) + (2*pow(s + x - z,2)*
									  (c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/pow(l2 - z2,2) -
					    (8*(s + x - z)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    pow(l2 - z2,2) + (2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
					    (l2 - z2));
    }
    else if (label_i=='z'&&label_j=='c'){
	result = (l2 - pow(s + x,2))*((2*r*pow(s + x - z,2)*z)/pow(l2 - z2,2) - (2*r*(s + x - z))/(l2 - z2) -
					    (7*pow(s + x - z,4)*(-24*l5*l6*r*z + 200*l4*l5*r*z3 - 456*l7*r*z5 + 280*l5*r*z7))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 + 35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,4)*z*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								    35*l5*r*z8))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								  35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(6*l6*l7*r - 24*l5*l6*r*z2 + 340*l4*l5*r*z4 - 952*l7*r*z6 +
								 630*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 + 70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,3)*z*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								    70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								  70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if (label_i=='z'&&label_j=='s'){
	result = (l2 - pow(s + x,2))*((8*r2*z2)/pow(l2 - z2,3) + (2*r2)/pow(l2 - z2,2) -
					    (7*pow(s + x - z,4)*(-105*l8 + 420*l6*z2 + 1050*l4*z4 - 2940*l2*z6 + 1575*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (28*pow(s + x - z,3)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
								  1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
								  1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
								  840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,4)*z*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (7*pow(s + x - z,3)*(-330*l8*z + 1640*l6*z3 - 540*l4*z5 - 2520*l2*z7 + 1750*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (28*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					      35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (168*pow(s + x - z,3)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								     35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (84*pow(s + x - z,2)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
								  24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
								  1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
								  3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
								  1750*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 - 315*l2*z8 +
					      175*z5*z5))/(l5*pow(l2 - z2,3)*
							       pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (42*pow(s + x - z,3)*z*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								    315*l2*z8 + 175*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								  315*l2*z8 + 175*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (21*pow(s + x - z,2)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
					     (-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
					    (126*pow(s + x - z,2)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
								     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
								     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
								     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
								     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (42*(s + x - z)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
							     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
							     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
							     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
							     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(s + x - z)*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
					    (l2 - z2) + (4*(s + x - z)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
										      r2/(l2 - z2)))/pow(l2 - z2,2) -
					    (2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2)) -
	    2*(s + x)*((8*r2*(s + x - z)*z2)/pow(l2 - z2,3) + (2*r2*(s + x - z))/pow(l2 - z2,2) -
		       (7*pow(s + x - z,4)*(-105*l8*s - 120*l8*z - 24*c*l5*l6*r*z + 180*l4*l5*r2*z + 420*l6*s*z2 +
					    1080*l6*z3 + 200*c*l4*l5*r*z3 - 1320*l7*r2*z3 + 1050*l4*s*z4 -
					    1800*l4*z5 - 456*c*l7*r*z5 + 2100*l5*r2*z5 - 2940*l2*s*z6 +
					    840*l2*z7 + 280*c*l5*r*z7 + 1575*s*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
			(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
			 12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
			 50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
			 76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
			 35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		       (42*pow(s + x - z,4)*z*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					       12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					       50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					       76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					       35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(-60*l5*l5 + 6*c*l6*l7*r + 72*l5*l6*r2 - 330*l8*s*z + 120*l8*z2 -
					    24*c*l5*l6*r*z2 + 1640*l6*s*z3 + 1200*l6*z4 + 340*c*l4*l5*r*z4 -
					    1560*l7*r2*z4 - 540*l4*s*z5 - 2520*l4*z6 - 952*c*l7*r*z6 +
					    3920*l5*r2*z6 - 2520*l2*s*z7 + 1260*l2*z8 + 630*c*l5*r*z8 +
					    1750*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (7*pow(s + x - z,3)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
			(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
			 40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
			 68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
			 136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
			 70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
		       (42*pow(s + x - z,3)*z*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					       40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					       68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					       136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					       70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (pow(s + x - z,2)*((16*r2*z3)/pow(l2 - z2,3) + (10*r2*z)/pow(l2 - z2,2)))/
		       (l2 - z2) + (2*pow(s + x - z,2)*z*(c*r + (4*r2*z2)/pow(l2 - z2,2) +
								      r2/(l2 - z2)))/pow(l2 - z2,2) -
		       (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2));
    }
    else if (label_i=='c'&&label_j=='c'){
	result = 0.;
    }
    else if (label_i=='c'&&label_j=='s'){
	result = (l2 - pow(s + x,2))*((2*r*(s + x - z))/(l2 - z2) -
					    (28*pow(s + x - z,3)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
								  35*l5*r*z8))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (21*pow(s + x - z,2)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
								  70*l5*r*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6))) -
	    2*(s + x)*((r*pow(s + x - z,2))/(l2 - z2) -
		       (7*pow(s + x - z,4)*(3*l6*l7*r - 12*l5*l6*r*z2 + 50*l4*l5*r*z4 - 76*l7*r*z6 +
					    35*l5*r*z8))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(6*l6*l7*r*z - 8*l5*l6*r*z3 + 68*l4*l5*r*z5 - 136*l7*r*z7 +
					    70*l5*r*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)));
    }
    else if(label_i=='s'&&label_j=='s'){
	result = (l2 - pow(s + x,2))*((-56*pow(s + x - z,3)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 -
								   420*l2*z7 + 175*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (84*pow(s + x - z,2)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
								  12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
								  50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
								  76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
								  35*c*l5*r*z8 + 175*s*z4*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*pow(s + x - z,2)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
								  315*l2*z8 + 175*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
					    (42*(s + x - z)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
							     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
							     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
							     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
							     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
					    (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
					    (2*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2)) -
	    4*(s + x)*((2*r2*z)/pow(l2 - z2,2) -
		       (7*pow(s + x - z,4)*(-105*l8*z + 140*l6*z3 + 210*l4*z5 - 420*l2*z7 + 175*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (28*pow(s + x - z,3)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
					     12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
					     50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
					     76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
					     35*c*l5*r*z8 + 175*s*z4*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (7*pow(s + x - z,3)*(-15*l5*l5 - 165*l8*z2 + 410*l6*z4 - 90*l4*z6 -
					    315*l2*z8 + 175*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		       (21*pow(s + x - z,2)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
					     40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
					     68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
					     136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
					     70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		       (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		       (2*(s + x - z)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/(l2 - z2))
	    - 2*((2*r2*(s + x - z)*z)/pow(l2 - z2,2) + r2/(l2 - z2) -
		 (7*pow(s + x - z,4)*(-15*l5*l5 + 3*c*l6*l7*r + 18*l5*l6*r2 - 105*l8*s*z - 60*l8*z2 -
				      12*c*l5*l6*r*z2 + 90*l4*l5*r2*z2 + 140*l6*s*z3 + 270*l6*z4 +
				      50*c*l4*l5*r*z4 - 330*l7*r2*z4 + 210*l4*s*z5 - 300*l4*z6 -
				      76*c*l7*r*z6 + 350*l5*r2*z6 - 420*l2*s*z7 + 105*l2*z8 +
				      35*c*l5*r*z8 + 175*s*z4*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
		 (7*pow(s + x - z,3)*(-15*l5*l5*s - 60*l5*l5*z + 6*c*l6*l7*r*z + 72*l5*l6*r2*z - 165*l8*s*z2 +
				      40*l8*z3 - 8*c*l5*l6*r*z3 + 410*l6*s*z4 + 240*l6*z5 +
				      68*c*l4*l5*r*z5 - 312*l7*r2*z5 - 90*l4*s*z6 - 360*l4*z7 -
				      136*c*l7*r*z7 + 560*l5*r2*z7 - 315*l2*s*z8 + 140*l2*z4*z5 +
				      70*c*l5*r*z4*z5 + 175*s*z5*z5))/
		 (l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
		 (pow(s + x - z,2)*(c*r + (4*r2*z2)/pow(l2 - z2,2) + r2/(l2 - z2)))/
		 (l2 - z2));
    }
    else{
        BOOST_ASSERT_MSG(false,"RhoDDerivative!");
    }
//	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
	return result;
}
double shape::IntegrateRhoDerivative(double x, char label){
    double l,r,z,c,s;
    double result = 0;
    l = _para_l;
    r = _para_r;
    z = _para_z;
    c = _para_c;
    s = _para_s;

    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4, l8=l4*l4;
    double z2=z*2, z3=z*z2, z4=z2*z2, z5=z2*z3, z6=z3*z3, z7=z3*z4, z8=z4*z4;
    double s2=s*2, s3=s*s2, s4=s2*s2, s5=s2*s3, s6=s3*s3, s7=s3*s4;

    if(label=='l'){
	result = (x*(612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		     756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				    4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		     28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
				    252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
					    21*z5)) - 300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
									  280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
									  70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
									  35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
									  84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
									  7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
									       168*z5)) + 360*l5*z3*
		     (196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
		      70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
		      7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
		      28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
				   93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
						       360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
		      2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
			   42*x*z5 + 42*z6)) - 40*l7*z*
		     (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
		      420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
		      14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
		      84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
				   334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
							 1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
		      8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
			   378*x*z5 + 63*z6)) + 60*l3*z5*
		     (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		      280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		      14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		      84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
				   59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
						       70*pow(x,2)*z4 + 49*z6) +
		      7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
			   414*x*z5 + 84*z6)) + 140*l4*r*z6*
		     (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
				  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
				  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			    35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			    7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
							 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
							 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		     84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
							 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							 42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							 pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					     6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
						  35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						  7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						  pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		     156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
					   200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
					 100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
					 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
					 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
					 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
					 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
					    396*z5))) + 108*l8*r*z2*
		     (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
			    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
			    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
			    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
			    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
		      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
				  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
				  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
				  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
				  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
				  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
				       30*z5))) + 88*l5*l5*r*(-(c*z2*
									  (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
									   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
									   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
									   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
									   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
									   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
										126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
												       450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
												       45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
												       27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
												       9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
						      245*z6)) - (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
									(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
									 105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
											 210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
											 pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
											 140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
									 63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
										       5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
										       4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
									 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
											       252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
											       210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
											       420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
											       35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
											       15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
												       21*z5)) - 30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
																     280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
																     70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
																     35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
																     84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
																     7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
																	  168*z5)) + 60*l6*z3*
									 (196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
									  70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
										       93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
													   360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
									  2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
									       42*x*z5 + 42*z6)) - 5*l8*z*
									 (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
									  420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
									  14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
									  84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
										       334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
													     1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
									  8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
									       378*x*z5 + 63*z6)) + 15*l4*z5*
									 (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
									  280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
									  14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
									  84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
										       59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
													   70*pow(x,2)*z4 + 49*z6) +
									  7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
									       414*x*z5 + 84*z6)) + 28*l5*r*z6*
									 (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
										      105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
										      35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
									  15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
										35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
										7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
									 36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
													    10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
													    10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
									 12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
													     70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
													     140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
													     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
													     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
												 6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
												      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
												      35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
												      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
												      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
									 12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
											      200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
											 c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
											    100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
											    15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
											    5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
											    3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
											    s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
											       396*z5))) + 12*l4*l5*r*z2*
									 (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
										168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
										210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
										14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
										42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
									  c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
										      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
										      70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
										      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
										      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
										      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
											   30*z5))) + 8*l5*l6*r*(-(c*z2*
															     (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
															      495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
															      30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
															      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
															      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
															      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
																   126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
																			  450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
																			  45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
																			  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
																			  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 -
							 245*z6,2)) - (x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
										105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
												210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
												pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
												140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
										63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
											      5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
											      4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
										14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
												      252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
												      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
												      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
												      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
												      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
													      21*z5)) - 30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
																	    280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
																	    70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
																	    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
																	    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
																	    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
																		 168*z5)) + 60*l6*z3*
										(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
										 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										 7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										 28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
											      93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
														  360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
										 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
										      42*x*z5 + 42*z6)) - 5*l8*z*
										(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
										 420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
										 14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
										 84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
											      334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
														    1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
										 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
										      378*x*z5 + 63*z6)) + 15*l4*z5*
										(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
										 280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
										 14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
										 84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
											      59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
														  70*pow(x,2)*z4 + 49*z6) +
										 7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
										      414*x*z5 + 84*z6)) + 28*l5*r*z6*
										(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
											     105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
											     35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
										 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
										       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
										       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
										36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
														   10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
														   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
										12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
														    70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
														    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
														    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
														    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
													6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
													     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
													     35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
													     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
													     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
										12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
												     200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
												c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
												   100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
												   15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
												   5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
												   3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
												   s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
												      396*z5))) + 12*l4*l5*r*z2*
										(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
										       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
										       210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
										       14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
										       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
										 c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
											     420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
											     70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
											     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
											     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
											     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
												  30*z5))) + 8*l5*l6*r*(-(c*z2*
																    (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
																     495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
																     30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
																     10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
																     6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
																     6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
																	  126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
																				 450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
																				 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
																				 27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
																				 9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
						     245*z6)) - (5*x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
									    105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
											    210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
											    pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
											    140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
									    63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
											  5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
											  4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
									    14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
												  252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
												  210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
												  420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
												  35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
												  15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
													  21*z5)) - 30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
																	280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
																	70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
																	35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
																	84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
																	7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
																	     168*z5)) + 60*l6*z3*
									    (196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
									     70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									     7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									     28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
											  93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
													      360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
									     2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
										  42*x*z5 + 42*z6)) - 5*l8*z*
									    (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
									     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
									     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
									     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
											  334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
														1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
									     8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
										  378*x*z5 + 63*z6)) + 15*l4*z5*
									    (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
									     280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
									     14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
									     84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
											  59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
													      70*pow(x,2)*z4 + 49*z6) +
									     7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
										  414*x*z5 + 84*z6)) + 28*l5*r*z6*
									    (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
											 105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
											 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
									     15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
										   35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
										   7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
									    36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
													       10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
													       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
									    12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
														70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
														140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
														42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
														pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
												    6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
													 35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
													 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
													 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
													 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
									    12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
												 200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
											    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
											       100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
											       15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
											       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
											       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
											       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
												  396*z5))) + 12*l4*l5*r*z2*
									    (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
										   168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
										   210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
										   14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
										   42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
									     c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
											 420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
											 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
											 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
											 105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
											 7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
											      30*z5))) + 8*l5*l6*r*(-(c*z2*
																(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
																 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
																 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
																 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
																 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
																 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
																      126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
																			     450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
																			     45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
																			     27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
																			     9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
						      245*z6));
//	cout<<x<<' '<<label<<' '<<result<<endl;
    }
    else if (label=='r'){
	result = (x*(108*l*l8*l8*r - 1080*l7*l8*r*z2 +
		     36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		     420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					      35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					      7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		     36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				     90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		     180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z +
					      84*pow(x,4)*z2 - 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					      210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					      14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		     72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) +
					     105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					     35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		     24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 -
				     765*pow(x,2)*z4 + 756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				     45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				     27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				     9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		     28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
						 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
						 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
						      10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						      10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		     12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
						       70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					   6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
						35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		     12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
					200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				   c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
				      100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
				      15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
					 396*z5))) + 12*l4*l5*z2*
		     (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
			    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
			    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
			    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
			    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
		      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
				  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
				  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
				  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
				  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
				  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
				       30*z5))) + 8*l5*l6*(-(c*z2*
								       (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
									495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
									30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
									10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
									6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
									6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
									     126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
												    450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
												    45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
												    27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
												    9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
						      245*z6));
//        cout<<x<<' '<<label<<' '<<result<<endl;
    }
    else if(label=='z') {
	result = (x*(36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*
		     (-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*s2*x*(x - z)*(2*x - z) +
		      70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) +
		      140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		     63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
				   4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
		     945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
				     210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
				     pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				     140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		     14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 +
					   350*pow(x,2)*z4 + 105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
					   420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
					   840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					   35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
					   15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		     30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 +
				   2520*x*z4 - 1680*z5 + 70*s4*(-12*x + 12*z) +
				   35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) + 84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
				   7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		     98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
					   252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					   210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					   420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					   35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					   15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
						   21*z5)) + 60*l6*z3*
		     (168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
		      7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
		      28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
		      7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 -
			   504*z5) - 2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 +
						210*x*z4 + 252*z5) - 2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 +
										    210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) +
		     180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
					    7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
							 93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
									     360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
					    2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
						 42*x*z5 + 42*z6)) + 15*l4*z5*
		     (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		      14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		      84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		      4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		      7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 +
			   504*z5) + 4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
					      70*pow(x,2)*z4 + 49*z6)) -
		     5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
				   14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				   84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				   21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 -
					 3024*z5) + 8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 -
							       1890*x*z4 + 378*z5) + 8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 -
												    630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) -
		     5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				 420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				 14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				 84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
					      334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
								    1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
				 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
				      378*x*z5 + 63*z6)) + 75*l4*z4*
		     (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		      280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		      14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		      84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
				   59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
						       70*pow(x,2)*z4 + 49*z6) +
		      7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
			   414*x*z5 + 84*z6)) - 36*l7*l8*r*
		     (60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 + 10*s2*(-3*x + 6*z) +
				  10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
		     28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
							 105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					     15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
						   7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					     2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						    105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						    35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		     168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
							  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
						    35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
						    7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 +
					  45*s*(-12*x*z + 24*z2)) - c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 +
									     680*pow(x,2)*z3 - 990*x*z4 + 852*z5 + 15*s4*(-7*x + 22*z) +
									     5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
									     3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
									     s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		     12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) +
							 14*s*x*pow(x - z,2)*(-8*x + 42*z) + 140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) -
							 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
							 pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					     6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
						  35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) -
						  210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
						  pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					     2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		     48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
							 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							 42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							 pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					     6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
						  35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						  7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						  pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		     8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 +
						  2160*pow(x,2)*z3 - 1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) +
						  10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
						  6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
						  6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				    3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 -
					 3570*z5 + 45*s3*(-60*x*z + 120*z2) +
					 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
					 9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				    2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
						126*z5))) + 12*l4*l5*r*z2*
		     (15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
			    420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
			    14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
			    42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
		      c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 -
				  525*x*z4 + 294*z5 + 70*s4*(-33*x + 42*z) +
				  70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
				  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
				  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
		      2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
			     420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
			     70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
			     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
			     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
			     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
				  30*z5))) + 24*l4*l5*r*z*(15*r*
								    (42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
								     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
								     210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
								     14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
								     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
								    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
										420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
										70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
										70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
										105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
										7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
										     30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
						      245*z6)) - (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
									(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
									 105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
											 210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
											 pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
											 140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
									 63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
										       5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
										       4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
									 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
											       252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
											       210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
											       420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
											       35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
											       15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
												       21*z5)) - 30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
																     280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
																     70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
																     35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
																     84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
																     7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
																	  168*z5)) + 60*l6*z3*
									 (196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
									  70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
										       93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
													   360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
									  2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
									       42*x*z5 + 42*z6)) - 5*l8*z*
									 (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
									  420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
									  14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
									  84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
										       334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
													     1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
									  8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
									       378*x*z5 + 63*z6)) + 15*l4*z5*
									 (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
									  280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
									  14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
									  84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
										       59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
													   70*pow(x,2)*z4 + 49*z6) +
									  7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
									       414*x*z5 + 84*z6)) + 28*l5*r*z6*
									 (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
										      105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
										      35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
									  15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
										35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
										7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
									 36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
													    10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
													    10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
									 12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
													     70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
													     140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
													     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
													     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
												 6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
												      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
												      35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
												      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
												      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
									 12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
											      200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
											 c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
											    100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
											    15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
											    5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
											    3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
											    s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
											       396*z5))) + 12*l4*l5*r*z2*
									 (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
										168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
										210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
										14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
										42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
									  c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
										      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
										      70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
										      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
										      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
										      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
											   30*z5))) + 8*l5*l6*r*(-(c*z2*
															     (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
															      495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
															      30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
															      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
															      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
															      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
																   126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
																			  450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
																			  45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
																			  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
																			  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 -
							 245*z6,2)) + (x*z*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
										  105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
												  210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
												  pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
												  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
										  63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
												5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
												4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
										  14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
													252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
													210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
													420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
													35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
													15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
														21*z5)) - 30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 +
																	      280*pow(x,3)*z3 - 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 +
																	      70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
																	      35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
																	      84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
																	      7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 +
																		   168*z5)) + 60*l6*z3*
										  (196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
										   70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										   7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										   28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
												93*z5) + 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
														    360*pow(x,2)*z4 + 261*x*z5 - 84*z6) -
										   2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
											42*x*z5 + 42*z6)) - 5*l8*z*
										  (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
										   420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
										   14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
										   84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
												334*z5) + 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
														      1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6) +
										   8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
											378*x*z5 + 63*z6)) + 15*l4*z5*
										  (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
										   280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
										   14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
										   84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
												59*z5) + 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 -
														    70*pow(x,2)*z4 + 49*z6) +
										   7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
											414*x*z5 + 84*z6)) + 28*l5*r*z6*
										  (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
											       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
											       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
										   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
											 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
											 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
										  36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z +
														     10*pow(x,2)*z2 - 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
														     10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
										  12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
														      70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
														      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
														      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
														      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
													  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
													       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
													       35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
													       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
													       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
										  12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 -
												       200*z4 + 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
												  c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 -
												     100*pow(x,3)*z3 + 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 +
												     15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
												     5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
												     3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
												     s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 -
													396*z5))) + 12*l4*l5*r*z2*
										  (15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
											 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
											 210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
											 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
											 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
										   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
											       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
											       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
											       70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
											       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
											       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
												    30*z5))) + 8*l5*l6*r*(-(c*z2*
																      (126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
																       495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
																       30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
																       10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
																       6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
																       6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
																	    126*z5))) + 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 +
																				   450*pow(x,3)*z3 - 765*pow(x,2)*z4 + 756*x*z5 - 595*z6 +
																				   45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
																				   27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
																				   9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
//        cout<<x<<' '<<label<<' '<<result<<endl;
    }
    else if(label=='c') {
        result = (x * (36 * l*l8*l8 * r *
                       (3 * s2 + 3 * s * x + pow(x, 2) - 6 * s * z - 3 * x * z + 3 * z2) +
                       28 * l5 * r * z8 *
                       (105 * s6 + 105 * s5 * (3 * x - 2 * z) + 105 * s * pow(x, 3) * pow(x - z, 2) +
                        105 * s4 * (5 * pow(x, 2) - 5 * x * z + z2) +
                        105 * s2 * pow(x, 2) * (3 * pow(x, 2) - 5 * x * z + 2 * z2) +
                        35 * s3 * x * (15 * pow(x, 2) - 20 * x * z + 6 * z2) +
                        pow(x, 4) * (15 * pow(x, 2) - 35 * x * z + 21 * z2)) -
                       36 * l7*l8 * r *
                       (10 * s4 + 2 * pow(x, 4) + 20 * s3 * (x - z) - 5 * pow(x, 3) * z +
                        10 * pow(x, 2) * z2 -
                        20 * x * z3 + 20 * z4 +
                        10 * s2 * (2 * pow(x, 2) - 3 * x * z + 3 * z2) +
                        10 * s * (pow(x, 3) - 2 * pow(x, 2) * z + 3 * x * z2 - 4 * z3)) -
                       12 * pow(l, 7) * r * z6 * (532 * s6 + 84 * s5 * (19 * x - 14 * z) +
                                                         70 * s4 *
                                                         (38 * pow(x, 2) - 42 * x * z + 15 * z2) +
                                                         14 * s * x * pow(x - z, 2) *
                                                         (38 * pow(x, 2) - 8 * x * z + 21 * z2) +
                                                         140 * s3 *
                                                         (19 * pow(x, 3) - 28 * pow(x, 2) * z + 15 * x * z2 -
                                                          5 * z3) +
                                                         42 * s2 * (38 * pow(x, 4) - 70 * pow(x, 3) * z +
                                                                           50 * pow(x, 2) * z2 -
                                                                           25 * x * z3 + 7 * z4) +
                                                         pow(x, 2) * (76 * pow(x, 4) - 196 * pow(x, 3) * z +
                                                                      210 * pow(x, 2) * z2 -
                                                                      175 * x * z3 + 98 * z4)) -
                       12 * pow(l, 13) * r *
                       (-21 * s6 - 3 * pow(x, 6) - 21 * s5 * (3 * x - 2 * z) + 7 * pow(x, 5) * z -
                        33 * pow(x, 4) * z2 +
                        100 * pow(x, 3) * z3 - 170 * pow(x, 2) * z4 + 198 * x * z5 -
                        142 * z6 -
                        15 * s4 * (7 * pow(x, 2) - 7 * x * z + 11 * z2) -
                        5 * s3 * (21 * pow(x, 3) - 28 * pow(x, 2) * z + 66 * x * z2 - 80 * z3) -
                        3 * s2 *
                        (21 * pow(x, 4) - 35 * pow(x, 3) * z + 110 * pow(x, 2) * z2 - 200 * x * z3 +
                         170 * z4) -
                        s * (21 * pow(x, 5) - 42 * pow(x, 4) * z + 165 * pow(x, 3) * z2 -
                             400 * pow(x, 2) * z3 + 510 * x * z4 - 396 * z5)
                       ) - 8 * l5*l6 * r * z2 *
                           (126 * s6 + 18 * pow(x, 6) + 42 * s5 * (9 * x - 10 * z) - 70 * pow(x, 5) * z +
                            240 * pow(x, 4) * z2 - 495 * pow(x, 3) * z3 + 540 * pow(x, 2) * z4 -
                            378 * x * z5 + 210 * z6 +
                            30 * s4 * (21 * pow(x, 2) - 35 * x * z + 40 * z2) +
                            10 * s3 *
                            (63 * pow(x, 3) - 140 * pow(x, 2) * z + 240 * x * z2 - 198 * z3) +
                            6 * s2 *
                            (63 * pow(x, 4) - 175 * pow(x, 3) * z + 400 * pow(x, 2) * z2 - 495 * x * z3 +
                             270 * z4) +
                            6 * s * (21 * pow(x, 5) - 70 * pow(x, 4) * z + 200 * pow(x, 3) * z2 -
                                     330 * pow(x, 2) * z3 + 270 * x * z4 -
                                     126 * z5)) + 12 * l4*l5 * r * z4 *
                                                         (350 * s6 + 50 * pow(x, 6) +
                                                          42 * s5 * (25 * x - 22 * z) - 154 * pow(x, 5) * z +
                                                          294 * pow(x, 4) * z2 -
                                                          420 * pow(x, 3) * z3 + 315 * pow(x, 2) * z4 -
                                                          105 * x * z5 + 49 * z6 +
                                                          70 * s4 *
                                                          (25 * pow(x, 2) - 33 * x * z + 21 * z2) +
                                                          70 * s3 *
                                                          (25 * pow(x, 3) - 44 * pow(x, 2) * z + 42 * x * z2 -
                                                           24 * z3) +
                                                          105 * s2 * (10 * pow(x, 4) - 22 * pow(x, 3) * z +
                                                                             28 * pow(x, 2) * z2 -
                                                                             24 * x * z3 + 9 * z4) +
                                                          7 * s * (50 * pow(x, 5) - 132 * pow(x, 4) * z +
                                                                   210 * pow(x, 3) * z2 -
                                                                   240 * pow(x, 2) * z3 + 135 * x * z4 -
                                                                   30 * z5)))) /
                 (12. * l5 * pow(l2 - z2, 3) *
                  (9 * pow(l, 6) - 63 * l4 * z2 - 21 * l2 * z4 -
                   245 * z6));
//        cout<<x<<' '<<label<<' '<<result<<endl;
    }
    else if(label=='s'){
            result = (x*(36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
                                                                (40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) +
                                                                 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) +
                         105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
                                           280*s3*(10*pow(x,2) - 15*x*z + 6*z2) +
                                           420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
                         105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) +
                                         210*s2*x*pow(x - z,2)*(2*x - z) + 70*s4*(10*pow(x,2) - 15*x*z + 6*z2) +
                                         pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
                                         140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
                         63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
                                         5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
                         14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) +
                                                   525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
                                                   840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
                                                   840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
                                                   105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
                                                   15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 -
                                                         21*z5)) - 30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) +
                                                                                           280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
                                                                                           105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
                                                                                           168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
                                                                                           7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5))
        - 5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
                          1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
                          42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
                          168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 +
                                 334*z5) + 21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 -
                                                       1080*pow(x,2)*z4 + 1068*x*z5 - 504*z6)) +
                         60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
                                                   35*s4*(140*pow(x,2) + 87*x*z - 288*z2) +
                                                   280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
                                                   21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
                                                   56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 +
                                                         93*z5) + 7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 -
                                                                             360*pow(x,2)*z4 + 261*x*z5 - 84*z6)) +
                         15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
                                                   420*s4*(70*pow(x,2) - 215*x*z + 186*z2) +
                                                   1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
                                                   42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
                                                   168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 -
                                                          59*z5) + 7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 +
                                                                              820*pow(x,2)*z4 - 414*x*z5 + 84*z6)) +
                         28*l5*r*z6*(c*z2*(630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) +
                                                                   420*s3*(5*pow(x,2) - 5*x*z + z2) + 210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
                                                                   105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
                                                     15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) +
                                                           140*s3*(10*pow(x,2) - 12*x*z + 3*z2) + 210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) +
                                                           7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
                         12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) +
                                                                   280*s3*(38*pow(x,2) - 42*x*z + 15*z2) + 14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
                                                                   420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
                                                                   84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
                                                     6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
                                                          140*s3*(55*pow(x,2) - 84*x*z + 45*z2) +
                                                          105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
                                                          7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
                         12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) +
                                                45*(pow(x,3) - 6*x*z2 + 8*z3)) -
                                           c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 -
                                              400*pow(x,2)*z3 + 510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
                                              15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
                                              6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
                         12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) +
                                                           840*s3*(pow(x,2) - 2*x*z + 2*z2) + 42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
                                                           42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
                                                           84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
                                                     c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
                                                                   210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
                                                                   210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
                                                                   7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 -
                                                                      30*z5))) + 8*l5*l6*r*(-(c*z2*
                                                                                                            (756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
                                                                                                             30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
                                                                                                             12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
                                                                                                             6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 -
                                                                                                                126*z5))) + 3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
                                                                                                                                         135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
                                                                                                                                         54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
                                                                                                                                         9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
                     (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 -
                                                                       245*z6));
//        cout<<x<<' '<<label<<' '<<result<<endl;
    }
    else{
        BOOST_ASSERT_MSG(false, "IntegrateRhoDerivative!");
    }
    result = result*pow(_Rcn,2);
    return result;
}
double shape::IntegrateRhoDDerivative(double x, char label_i, char label_j) {
    double l,r,z,c,s;
    double result = 0;
    l = _para_l;
    r = _para_r;
    z = _para_z;
    c = _para_c;
    s = _para_s;

    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4, l8=l4*l4;
    double z2=z*2, z3=z*z2, z4=z2*z2, z5=z2*z3, z6=z3*z3, z7=z3*z4, z8=z4*z4;
    double s2=s*2, s3=s*s2, s4=s2*s2, s5=s2*s3, s6=s3*s3, s7=s3*s4;

    if (label_i=='l'&&label_j=='l'){
        result = (x*(9792*l7*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
                     8316*l5*l5*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				     5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				     4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
                     28*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				  70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				  210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				  420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				  35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				  15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
                     2700*l8*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
                     1800*l4*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					     70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					     7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					     28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					     7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
						  84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								      175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
                     280*l6*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				     21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
					   1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
										  700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
                     180*l2*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
					    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
					    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
					    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
						 414*x*z5 + 84*z6)) + 560*l3*r*z6*
		     (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
				  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
				  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			    35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			    7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
                     7560*l6*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
							  20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
							  10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
                     504*l5*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) +
							  70*s4*(38*pow(x,2) - 42*x*z + 15*z2) + 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
						   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
                     1872*l5*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					    90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				       c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
					  170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
					  5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
					  3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
					  s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
                     864*l7*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
						    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
						    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
							  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
							  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
							  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
                     880*l4*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
						   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
					  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
					  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			     5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			     4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
	      28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
			     70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
			     210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
			     420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
			     35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
			     15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      360*l5*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				     70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				     7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				     28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				     7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					  84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
							      175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
	      40*l7*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			     21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				   1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									  700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
	      60*l3*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					 414*x*z5 + 84*z6)) + 140*l4*r*z6*
	      (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
	       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		     35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		     7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						  20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						  10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				    90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			       c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				  170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				  5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				  3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				  s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      108*l8*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					     210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					     14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				       c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						   70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						   70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      88*l5*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					    495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					    30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			      3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				   756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				   45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				   27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				   9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*(612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			       5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			       4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
			       70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
			       210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
			       420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
			       35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
			       15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			      490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			      35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			      84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			      7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		360*l5*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				       70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				       7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				       28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				       7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					    84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
		40*l7*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				     1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									    700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
		60*l3*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				      280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				      14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				      84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				      4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				      7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					   414*x*z5 + 84*z6)) + 140*l4*r*z6*
		(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			     105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			     35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						    20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						    10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					     35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				      90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				 c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				    170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				    5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				    3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				    s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		108*l8*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					       210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					       14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					 c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						     420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						     70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		88*l5*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					      495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					      30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				     756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				     45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				     27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				     9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*(612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				 5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		  28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				 210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				 420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				 35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				 15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		  360*l5*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					 7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					 28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					      84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								  175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
		  40*l7*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				 420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				 14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				 84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				       1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									      700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
		  60*l3*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
					280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
					4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
					7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					     414*x*z5 + 84*z6)) + 140*l4*r*z6*
		  (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						      20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						      10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					       35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				   c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				      170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  108*l8*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
						 210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
						 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						       70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  88*l5*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
						10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				  3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				       756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				       45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				       27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				       9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (x*pow(54*l5 - 252*l3*z2 - 42*l*z4,2)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			    4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
	      14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
				    252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
							     175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
	      5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				  1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									 700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
	      15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					 414*x*z5 + 84*z6)) + 28*l5*r*z6*
	      (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
	       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		     35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		     7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) -
	    (x*(270*l4 - 756*l2*z2 - 42*z4)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			    4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
	      14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
				    252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
							     175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
	      5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				  1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									 700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
	      15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					 414*x*z5 + 84*z6)) + 28*l5*r*z6*
	      (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
	       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		     35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		     7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			    4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
	      14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
				    252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
							     175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
	      5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				  1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									 700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
	      15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					 414*x*z5 + 84*z6)) + 28*l5*r*z6*
	      (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
	       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		     35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		     7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (5*x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
			    4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
	      14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
				    252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
							     175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
	      5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				  1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									 700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
	      15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
				    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
				    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
				    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
				    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
				    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					 414*x*z5 + 84*z6)) + 28*l5*r*z6*
	      (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
	       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
		     35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
		     7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					   35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
					    210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
					    14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (4*x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				  70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		  14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
					252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		  60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					     84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								 175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
		  5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				      1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									     700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
		  15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
					280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
					4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
					7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					     414*x*z5 + 84*z6)) + 28*l5*r*z6*
		  (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						     20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						     10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					       35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				       90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				  c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				     170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				     5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				     3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				     s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
						210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
						14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					  c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						      70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					       495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					       30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					       10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					       6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					       6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				      45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l3*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (9*x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				  70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		  14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
					252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		  60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					     84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								 175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
		  5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				      1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									     700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
		  15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
					280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
					4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
					7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					     414*x*z5 + 84*z6)) + 28*l5*r*z6*
		  (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						     20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						     10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					       35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				       90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				  c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				     170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				     5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				     3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				     s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
						210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
						14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					  c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						      70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					       495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					       30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					       10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					       6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					       6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				      45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (5*x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				  70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) +
				4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4)) -
		  14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 -
					252*pow(x,4)*z3 + 70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		  60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					     84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 -
								 175*pow(x,2)*z4 + 42*x*z5 + 42*z6)) -
		  5*l8*z*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
				420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 +
				      1068*x*z5 - 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 +
									     700*pow(x,2)*z4 - 378*x*z5 + 63*z6)) +
		  15*l4*z5*(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
					280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
					4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
					7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 -
					     414*x*z5 + 84*z6)) + 28*l5*r*z6*
		  (c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
			       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
			       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
		   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
			 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
			 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						     20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						     10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) +
					       35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				       90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				  c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				     170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				     5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				     3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				     s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 +
						210*s4*(pow(x,2) - 2*x*z + 2*z2) + 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) +
						14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					  c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						      70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					       495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					       30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					       10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					       6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					       6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				 3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				      45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l7*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if(label_i=='l'&&label_j=='r'){
	result = (x*(1836*l8*l8*r - 16200*l7*l7*r*z2 + 612*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		     2100*l4*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		     468*l6*l6*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				      90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		     1620*l8*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					       42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		     504*l6*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		     264*l5*l5*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		     140*l4*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
							105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     540*l7*l7*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						       20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		     84*l6*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						       14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					   6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		     156*l6*l6*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				       170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		     108*l8*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
						  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
							420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
							70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     88*l5*l5*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				   3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
					756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
					45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (108*l*l8*l8*r - 1080*l7*l8*r*z2 + 36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
				       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
				       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
	      36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
			      90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
	      180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
				       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				       42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
	      72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
				      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
				      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
				      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
	      24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
			      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
			      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
			      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
	      28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
				    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
					       20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
					       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				    6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					 35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
			       170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
			       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
			       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
			       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			   3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*(108*l*l8*l8*r - 1080*l7*l8*r*z2 + 36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
				      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					    70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					    pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		12*l4*l5*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		8*l5*l6*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*(108*l*l8*l8*r - 1080*l7*l8*r*z2 + 36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					   35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					   7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		  36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				  90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		  180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					   168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					   42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					   42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		  72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					  35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					  7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					  pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		  24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		  28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						    105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						    35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					      70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					      pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						   20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				     90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				   170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				   5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				   3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				   s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  12*l4*l5*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					      168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					      42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						    420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						    70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						    105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						    7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					     495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					     30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					     6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					     6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			       3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				    756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				    45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				    27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				    9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='l'&&label_j=='z'){
	result = (x*(612*c*l8*l8*r*(-6*s - 3*x + 6*z) + 756*l5*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
									   4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) -
		     28*l*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
				    105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
				    420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
				    840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
				    15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		     300*l4*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
				   70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
				   84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
				   7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		     196*l*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				     70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				     210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				     420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				     35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				     15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		     360*l5*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
					    7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
					    28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
					    7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
					    2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
					    2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
					       42*z6)) + 1080*l5*z2*(196*s7 + 84*s6*(7*x + 2*z) +
										       7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										       7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										       28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
										       7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
											    84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
														42*x*z5 + 42*z6)) + 60*l3*z5*
		     (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		      14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		      84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		      4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		      7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
		      4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
		     40*l7*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
				    14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				    84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				    21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
				    8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
				    8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
				       63*z6)) - 40*l7*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
								    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
								    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
								    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
								    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
									  504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
											       378*x*z5 + 63*z6)) + 300*l3*z4*
		     (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		      280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		      14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		      84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		      4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		      7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			   84*z6)) - 540*l7*l7*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
									10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
		     140*l4*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
							  105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					      15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
						    7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					      2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
						     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
						     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		     840*l4*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
							  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						    70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						    pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     156*l6*l6*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
				      c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
					 15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
					 3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
					 s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		     84*l6*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
							 140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
							 pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					     6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
						  35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
						  pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					     2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		     336*l6*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
							  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		     88*l5*l5*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
						   1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
						   6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
						   6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				     3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
					  45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
					  9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				     2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
					    540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		     108*l8*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
						    420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
						    14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
						    42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					      c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
							  294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
							  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
							  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					      2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
						     315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     216*l8*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					     42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				       c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						   70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*(-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) -
								    420*s2*x*(x - z)*(2*x - z) + 70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) +
								    140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
	      63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
			    4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
	      945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
	      14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
				    105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
				    420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
				    840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
				    15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
	      30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
			    70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
			    84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
			    7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
	      98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				    70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
	      60*l6*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
				    7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
				    28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
				    7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
				    2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
				    2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
				       42*z6)) + 180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
									      7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									      7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									      28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									      7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										   84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
												       42*x*z5 + 42*z6)) + 15*l4*z5*
	      (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
	       14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
	       84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
	       4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
	       7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
	       4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
	      5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
			    14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
			    84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
			    21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
			    8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
			    8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
			       63*z6)) - 5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
							   420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
							   14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
							   84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
							   21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
								 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
										      378*x*z5 + 63*z6)) + 75*l4*z4*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) - 36*l7*l8*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
	      28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						  105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
				      15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
					    7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
				      2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
					     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
					     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
	      168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
				       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					     70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					     pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
			      c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				 15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				 3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				 s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
	      12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						  140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						  pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
				      6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					   35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					   pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
				      2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					     14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
					     140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
	      48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
	      8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					   1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					   6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					   6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
			     3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				  45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				  9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
			     2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
				    540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
				    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
				    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
				    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
	      12*l4*l5*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
					    420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
					    14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
					    42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
				      c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						  294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
				      2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
					     315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      24*l4*l5*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
				     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				     42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
			       c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
					   70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*(36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*(-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) -
								      420*s2*x*(x - z)*(2*x - z) + 70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) +
								      140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
			      4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
		945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
				      105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
				      420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
				      840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
				      15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
			      70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
			      84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
			      7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		60*l6*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
				      7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
				      28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
				      7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
				      2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
				      2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
					 42*z6)) + 180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
										7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
										7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										     84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
													 42*x*z5 + 42*z6)) + 15*l4*z5*
		(-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		 14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		 84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		 4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		 7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
		 4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
		5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
			      14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
			      84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
			      21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
			      8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
			      8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
				 63*z6)) - 5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
							     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
							     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
							     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
							     21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
								   504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
											378*x*z5 + 63*z6)) + 75*l4*z4*
		(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		 280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		 14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		 84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		 7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		      84*z6)) - 36*l7*l8*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								  10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
		28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						    105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
					      7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
					       105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
					       pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						     105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						     35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					       70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					       pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
				c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				   15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				   3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				   s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						    140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						    pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					     35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					     pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					       14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
					       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					     1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					     6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					     6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
			       3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				    45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				    9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
			       2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
				      540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
				      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
				      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
				      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		12*l4*l5*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
					      420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
					      14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
					      42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						    294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						    105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						    7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
					       315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					       70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		24*l4*l5*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
				       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				       42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				 c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					     420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
					     70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (2.*l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*(36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*
		  (-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*s2*x*(x - z)*(2*x - z) +
		   70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) + 140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		  63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
				4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
		  945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				  70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		  14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
					105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
					420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
					840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
					15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		  30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
				70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
				84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
				7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		  98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		  60*l6*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
					7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
					28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
					7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
					2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
					2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
					   42*z6)) + 180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
										  7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
										  7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										       84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
													   42*x*z5 + 42*z6)) + 15*l4*z5*
		  (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		   14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		   84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		   4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		   7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
		   4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
		  5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
				14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
				8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
				8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
				   63*z6)) - 5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
							       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
							       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
							       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
							       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
								     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
											  378*x*z5 + 63*z6)) + 75*l4*z4*
		  (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		   280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		   14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		   84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		   4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		   7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			84*z6)) - 36*l7*l8*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								    10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
		  28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						      105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					  15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
						7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					  2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
						 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
						 pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		  168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						 pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
				  c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				     15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				     3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				     s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		  12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						      140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						      pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					  6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					       35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					       pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					  2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						 42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						 pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		  48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		  8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					       1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					       6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					       6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				 3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				      45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				      9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				 2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
					540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		  12*l4*l5*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
						420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
						14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
						42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					  c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						      294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						      105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						      7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					  2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
						 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						 105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						 7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  24*l4*l5*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
					       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			     5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
			       70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
			       210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
			       420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
			       35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
			       15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      360*l5*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				     70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				     7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				     28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				     7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					  84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							      42*x*z5 + 42*z6)) - 40*l7*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 60*l3*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 140*l4*r*z6*(c*z2*
							     (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							      105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							      pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							     15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								   70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								   pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						  20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						  10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				    90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			       c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				  170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				  5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				  3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				  s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      108*l8*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					     42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				       c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						   70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      88*l5*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					    495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					    30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			      3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				   756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				   45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				   27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				   9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*z*(612*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  756*l5*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				 5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
		  - 28*l*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				   70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				   210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				   420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				   35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				   15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  300*l4*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		  360*l5*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					 7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					 28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					      84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
								  42*x*z5 + 42*z6)) - 40*l7*z*
		  (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
		   420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
		   14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
		   84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
		   21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
			 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					      378*x*z5 + 63*z6)) + 60*l3*z5*
		  (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		   280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		   14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		   84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		   4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		   7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			84*z6)) + 140*l4*r*z6*(c*z2*
								 (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
								  105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
								  pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
								 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								       70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								       pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  540*l7*l7*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						      20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						      10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  84*l6*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  156*l6*l6*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				   c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				      170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  108*l8*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
						 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  88*l5*l5*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				  3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				       756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				       45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				       27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				       9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) -
	    (x*(-504*l3*z - 168*l*z3)*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
						   105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
								   70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
								   140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
						   63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
								 5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
						   - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
									   70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
									   210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
									   420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
									   35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
									   15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
						   30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
								 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
								 35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
								 84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
								 7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
						   60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
									 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									 7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									 28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									      84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
												  42*x*z5 + 42*z6)) - 5*l8*z*
						   (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
						    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
						    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
						    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
						    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
							  504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
									       378*x*z5 + 63*z6)) + 15*l4*z5*
						   (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
						    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
						    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
						    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
						    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
						    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
							 84*z6)) + 28*l5*r*z6*(c*z2*
												 (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
												  105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
												  pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
												 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
												       70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
												       pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
						   36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
										      20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
										      10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
						   12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
										       14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
										       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
										       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
										       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
									   6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
										35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
										7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
										pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
						   12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
									90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
								   c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
								      170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
								      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
								      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
								      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
						   12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
										 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
										 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
										 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
									   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
										       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
										       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
										       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
										       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
						   8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
										495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
										30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
										6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
										6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
								  3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
								       756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
								       45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
								       27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
								       9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*z*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l4*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (5*x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l6*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (4*x*z*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		    105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				    70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				    140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		    63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				  5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
		    - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					    70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		    30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				  490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				  35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				  84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				  7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		    60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					  70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					  7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					       84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
								   42*x*z5 + 42*z6)) - 5*l8*z*
		    (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
		     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
		     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
		     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
		     21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
			   504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
						378*x*z5 + 63*z6)) + 15*l4*z5*
		    (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		     280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		     14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		     84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		     4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		     7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			  84*z6)) + 28*l5*r*z6*(c*z2*
								  (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
								   105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
								   pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
								  15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
									70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
									pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		    36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						       20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		    12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
							14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					    6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						 35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		    12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				       170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		    12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
						  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
							420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
							70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		    8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				   3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
					756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
					45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l4*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*z*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		    105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				    70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				    140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		    63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
				  5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
		    - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					    70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		    30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
				  490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
				  35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				  84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				  7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		    60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
					  70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
					  7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					       84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
								   42*x*z5 + 42*z6)) - 5*l8*z*
		    (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
		     420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
		     14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
		     84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
		     21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
			   504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
						378*x*z5 + 63*z6)) + 15*l4*z5*
		    (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		     280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		     14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		     84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		     4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		     7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			  84*z6)) + 28*l5*r*z6*(c*z2*
								  (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
								   105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
								   pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
								  15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
									70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
									pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		    36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						       20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		    12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
							14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					    6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						 35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		    12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
					 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				       170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		    12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
						  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
						  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
							420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
							70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		    8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
						 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
						 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				   3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
					756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
					45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l6*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='l'&&label_j=='c'){
	result = (x*(612*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		     140*l4*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					      105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					      35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		     540*l7*l7*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				      10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		     84*l6*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					     14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		     156*l6*l6*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				      170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		     88*l5*l5*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					      495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		     108*l8*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*l*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
	      28*l5*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
				      105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
				      35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
	      36*l7*l8*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
			      10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
	      12*l7*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
				      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
				      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
				      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
	      12*l6*l7*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
			      170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
			      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
			      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
			      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
	      8*l5*l6*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
				      495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
				      10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
				      6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
				      6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
	      12*l4*l5*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
				      420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
				      70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
				      105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
				      7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*(36*l*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		28*l5*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		36*l7*l8*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		12*l7*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		12*l6*l7*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		8*l5*l6*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		12*l4*l5*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (2.*l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*(36*l*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		  28*l5*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		  36*l7*l8*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				  10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		  12*l7*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		  12*l6*l7*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				  170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				  5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				  3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				  s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		  8*l5*l6*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					  495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					  10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					  6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					  6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		  12*l4*l5*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='l'&&label_j=='s'){
	result = (x*(612*c*l8*l8*r*(6*s + 3*x - 6*z) - 540*c*l7*l7*r*(40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) +
									      10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) +
		     756*l5*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
				    5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
		     28*l*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
				    840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		     300*l4*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
				   105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				   168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				   7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
		     40*l7*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
				    1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				    42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				    168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				    21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
					504*z6)) + 360*l5*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
										35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
										7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										   84*z6)) + 60*l3*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
															 420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
															 42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
															 168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
															 7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
															    84*z6)) + 140*l4*r*z6*(c*z2*
																				     (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																				      210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																				     15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																					   210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		     84*l6*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
							 14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							 84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					     6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						  7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
		     156*l6*l6*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
				      c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
					 510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
					 15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
					 6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		     108*l8*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
						    42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						    84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					      c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
							  210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							  210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							  7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     88*l5*l5*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
						   30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						   12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						   6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				     3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
					  135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					  54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					  9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(54*l5 - 252*l3*z2 - 42*l*z4)*
	     (36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
	      (40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) + 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))
	      + 105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			    70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			    140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
			    5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
	      14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
				    840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
			    105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
	      5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
			    1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				504*z6)) + 60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
								       35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
								       21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
								       56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
								       7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									  84*z6)) + 15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
														420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
														   84*z6)) + 28*l5*r*z6*(c*z2*
																			   (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																			    210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																			   15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																				 210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
	      12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
				      6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
	      12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				 510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				 15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
	      12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
					    42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						  210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					   30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				  135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*(36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
		(40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) + 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))
		+ 105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				  280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
			      5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
		14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
				      840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
			      105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			      168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			      7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
		5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
			      1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			      42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			      168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			      21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				  504*z6)) + 60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
									 35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									 21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									 56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									 7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									    84*z6)) + 15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
														  420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														  42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														  168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														  7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
														     84*z6)) + 28*l5*r*z6*(c*z2*
																			     (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																			      210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																			     15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																				   210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
		12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
				c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				   510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				   15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				   6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
					      42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						    210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						    210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						    7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		8*l5*l6*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					     30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					     12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					     6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			       3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				    135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				    54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				    9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l4*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (5*x*(36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
		  (40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) + 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))
		  + 105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				    280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
				5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
		  14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
					840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
				105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
		  5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
				1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				    504*z6)) + 60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
									   35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									   21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									   56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									   7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									      84*z6)) + 15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
														    420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														    42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														    168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														    7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
														       84*z6)) + 28*l5*r*z6*(c*z2*
																			       (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																				210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																			       15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																				     210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		  12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					  6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
		  12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
				  c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				     510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				     15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				     6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		  12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
						42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					  c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						      210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						      210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						      7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					       30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					       12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					       6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				 3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				      135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l6*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='r'&&label_j=='r'){
	result = (x*(216*l*l8*l8 - 2160*l7*l8*z2 + 840*l5*z6*
		     (70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
		      70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
		      pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		     72*l6*l7*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		     360*l4*l5*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		     144*l7*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					    35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					    7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					    pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		     48*l5*l6*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				   756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				   27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				   9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='r'&&label_j=='z'){
	result = (x*(-2160*l7*l8*r*z + 36*c*l*l8*l8*(-6*s - 3*x + 6*z) +
		     420*l5*r*z6*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
					      7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
		     2520*l5*r*z5*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		     36*l6*l7*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) +
		     180*l4*l5*r*z2*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
					      420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
					      14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
					      42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) -
		     72*l7*r*z4*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					     35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					     pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
		     360*l4*l5*r*z*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 - 168*pow(x,3)*z3 +
				       182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				       42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		     288*l7*r*z3*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		     24*l5*l6*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 -
				     3570*z5 + 45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				     9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
		     36*l7*l8*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 + 10*s2*(-3*x + 6*z) +
					       10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) + 28*l5*z6*
		     (c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) + 105*s2*pow(x,2)*(-5*x + 4*z) +
				  35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
		      15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
			    7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
		      2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
			     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
			     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		     168*l5*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
							105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     12*l6*l7*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
				   c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				      15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				      3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				      s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		     12*l7*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						       140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						       42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						       pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					   6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
						35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
						7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
						pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					   2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		     48*l7*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						       14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					   6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		     8*l5*l6*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
						1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
						6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
						6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				  3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				       45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				       9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				  2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
					 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		     12*l4*l5*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
						 420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
						 14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
						 42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					   c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						       294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						       105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						       7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					   2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
						  315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     24*l4*l5*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (108*l*l8*l8*r - 1080*l7*l8*r*z2 + 36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
				       35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
				       7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
	      36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
			      90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
	      180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
				       168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				       42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				       42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
	      72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
				      35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
				      7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
				      pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
	      24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
			      756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
			      27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
			      9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
	      28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
				    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
					       20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
					       10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				    6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					 35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					 7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					 pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				 90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			    c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
			       170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
			       5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
			       3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
			       s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					 495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			   3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*z*(108*l*l8*l8*r - 1080*l7*l8*r*z2 + 36*l*l8*l8*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		  420*l5*r*z6*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) +
					   35*s4*(10*pow(x,2) - 12*x*z + 3*z2) + 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) +
					   7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) + pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2)) -
		  36*l6*l7*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				  90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) +
		  180*l4*l5*r*z2*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					   168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					   42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					   42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		  72*l7*r*z4*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					  35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					  7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					  pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4)) +
		  24*l5*l6*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) + 45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		  28*l5*z6*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						    105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						    35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					      70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					      pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  36*l7*l8*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						   20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		  12*l7*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		  12*l6*l7*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				     90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				   170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				   5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				   3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				   s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		  12*l4*l5*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					      168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					      42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						    420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						    70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						    105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						    7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					     495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					     30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					     6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					     6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			       3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				    756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				    45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				    27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				    9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='r'&&label_j=='c'){
	result = (x*(36*l*l8*l8*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		     28*l5*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		     36*l7*l8*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				   10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		     12*l7*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					   14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					   42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					   pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		     12*l6*l7*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				   170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				   5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				   3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				   s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		     8*l5*l6*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					   10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		     12*l4*l5*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					   70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='r'&&label_j=='s'){
	result = (x*(36*c*l*l8*l8*(6*s + 3*x - 6*z) + 420*l5*r*z6*
		     (420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
		      210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2)) -
		     36*c*l7*l8*(40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) +
				     10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		     36*l6*l7*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) +
		     180*l4*l5*r*z2*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
					      42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) -
		     72*l7*r*z4*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4)) +
		     24*l5*l6*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				     135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) + 54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				     9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)) +
		     28*l5*z6*(c*z2*(630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) +
						       420*s3*(5*pow(x,2) - 5*x*z + z2) + 210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						       105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
					   15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
						 210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		     12*l7*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						       14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						       84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					   6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
		     12*l6*l7*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
				   c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				      510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				      15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				      6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		     12*l4*l5*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
						 42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						 84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					   c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						       210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						       210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						       7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     8*l5*l6*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
						30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
						6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				  3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				       135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				       54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				       9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='z'&&label_j=='z'){
	result = (x*(216*c*l*l8*l8*r + 105*s*(840*s4 + pow(x,3)*(168*x - 210*z) + 140*s3*(12*x - 6*z) + 840*s2*x*(x - z) +
						840*s*pow(x,2)*(x - z) + 420*s2*x*(2*x - z))*z4*z5 +
		     63*l6*l6*(-120*s2 + 5*s*(-36*x + 168*z) + 4*(-20*pow(x,2) + 120*x*z - 180*z2)) +
		     1890*s*z8*(-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*s2*x*(x - z)*(2*x - z) +
				      70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) + 140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) -
		     36*l7*l8*r*(60*r + c*(60*s2 + 20*pow(x,2) + 10*s*(6*x - 24*z) - 120*x*z + 240*z2)) -
		     14*l2*z7*(24360*s5 + 560*pow(x,5) + 210*s4*(240*x - 252*z) + 420*s2*(10*x - 12*z)*pow(x - z,2) -
					   1512*pow(x,4)*z + 1400*pow(x,2)*z3 - 1680*s2*(x - z)*(-12*pow(x,2) + 10*x*z - 6*z2) +
					   35*s3*(1520*pow(x,2) - 2376*x*z + 1296*z2) + 15*s*x*(504*pow(x,3) - 1260*pow(x,2)*z + 1008*x*z2 - 420*z3) +
					   840*s2*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3)) +
		     7560*s*z7*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		     30*l5*l5*(840*s4 - 168*pow(x,4) + 1680*pow(x,3)*z - 5880*pow(x,2)*z2 + 10080*x*z3 - 8400*z4 +
				   35*s3*(24*x + 120*z) + 84*s2*(90*x*z - 240*z2) + 7*s*(-60*pow(x,3) + 840*pow(x,2)*z - 2700*x*z2 + 3360*z3)
			 ) - 196*l2*z6*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
						    105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
						    420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
						    840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
						    35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
						    15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		     588*l2*z5*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					    70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		     60*l6*z3*(-4032*s5 + 70*s4*(-132*x + 336*z) + 7*s3*(-1600*pow(x,2) + 5760*x*z - 7200*z2) +
					   28*s2*(-270*pow(x,3) + 1200*pow(x,2)*z - 2160*x*z2 + 1860*z3) +
					   7*s*(-384*pow(x,4) + 1920*pow(x,3)*z - 4320*pow(x,2)*z2 + 5220*x*z3 - 2520*z4) -
					   2*z*(-336*pow(x,4) + 1260*pow(x,3)*z - 2100*pow(x,2)*z2 + 840*x*z3 + 1260*z4) -
					   4*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5)) +
		     15*l4*z5*(31248*s5 + 280*s4*(240*x - 396*z) + 14*s3*(5360*pow(x,2) - 12870*x*z + 10920*z2) +
					   84*s2*(540*pow(x,3) - 1650*pow(x,2)*z + 1980*x*z2 - 1180*z3) +
					   4*z*(-462*pow(x,4) + 840*pow(x,3)*z - 840*pow(x,2)*z2 + 1470*z4) +
					   7*s*(1968*pow(x,4) - 6930*pow(x,3)*z + 9840*pow(x,2)*z2 - 8280*x*z3 + 2520*z4) +
					   8*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5)) -
		     5*l8*z*(-168*s5 + 420*s4*(-12*x + 162*z) + 14*s3*(-920*pow(x,2) + 9990*x*z - 23760*z2) +
				   84*s2*(-170*pow(x,3) + 1710*pow(x,2)*z - 5400*x*z2 + 6680*z3) +
				   21*s*(-360*pow(x,4) + 3510*pow(x,3)*z - 12960*pow(x,2)*z2 + 21360*x*z3 - 15120*z4) +
				   8*z*(630*pow(x,4) - 3780*pow(x,3)*z + 8400*pow(x,2)*z2 - 7560*x*z3 + 1890*z4) +
				   16*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5)) +
		     360*l6*z2*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
					    7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
					    28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
					    7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
					    2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
					    2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
					       42*z6)) + 360*l6*z*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
									       70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									       7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									       28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									       7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										    84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
													42*x*z5 + 42*z6)) + 150*l4*z4*
		     (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		      14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		      84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		      4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		      7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
		      4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
		     10*l8*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
				  14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				  84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				  21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
				  8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
				  8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
				     63*z6)) + 300*l4*z3*(1176*s7 + 168*s6*(21*x - 41*z) +
									    84*s5*(70*pow(x,2) - 215*x*z + 186*z2) + 280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
									    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
									    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
									    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
									    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
										 84*z6)) + 28*l5*r*z6*(15*r*
															 (210*s4 + 420*s3*x + 420*s2*pow(x,2) + 210*s*pow(x,3) + 42*pow(x,4)) +
															 c*(210*s4 + 420*s3*x + 420*s2*pow(x,2) + 210*s*pow(x,3) + 42*pow(x,4))*z2 +
															 4*c*z*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) + 105*s2*pow(x,2)*(-5*x + 4*z) +
																35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
															 2*c*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
															      105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
															      pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		     336*l5*r*z5*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
							  105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					      15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
						    7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					      2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
						     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
						     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		     840*l5*r*z4*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
							  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						    70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						    pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		     12*l6*l7*r*(3*r*(-540*s2 - 180*pow(x,2) + 1080*x*z - 2400*z2 + 45*s*(-12*x + 48*z)) -
				     c*(330*s4 + 66*pow(x,4) + 5*s3*(132*x - 480*z) - 600*pow(x,3)*z + 2040*pow(x,2)*z2 - 3960*x*z3 +
					4260*z4 + 3*s2*(220*pow(x,2) - 1200*x*z + 2040*z2) +
					s*(330*pow(x,3) - 2400*pow(x,2)*z + 6120*x*z2 - 7920*z3))) -
		     12*l7*r*z4*(c*z2*(2100*s4 + 140*s3*(30*x - 30*z) + 588*s*x*pow(x - z,2) - 56*s*x*(x - z)*(-8*x + 42*z) +
							 28*s*x*(38*pow(x,2) - 8*x*z + 21*z2) + 42*s2*(100*pow(x,2) - 150*x*z + 84*z2) +
							 pow(x,2)*(420*pow(x,2) - 1050*x*z + 1176*z2)) +
					     6*r*(3150*s4 + 35*s3*(180*x - 240*z) + 1470*s2*pow(x - z,2) - 420*s2*(x - z)*(-6*x + 14*z) +
						  210*s2*(11*pow(x,2) - 6*x*z + 7*z2) + 7*s*x*(450*pow(x,2) - 1200*x*z + 1260*z2) +
						  pow(x,2)*(630*pow(x,2) - 2100*x*z + 2940*z2)) +
					     4*c*z*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						    140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						    pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					     2*c*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		     96*l7*r*z3*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
							 140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
							 pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					     6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
						  35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
						  pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					     2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		     144*l7*r*z2*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
							  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
						   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		     8*l5*l6*r*(-(c*z2*(2400*s4 + 480*pow(x,4) + 10*s3*(480*x - 1188*z) - 2970*pow(x,3)*z + 6480*pow(x,2)*z2 -
						  7560*x*z3 + 6300*z4 + 6*s2*(800*pow(x,2) - 2970*x*z + 3240*z2) +
						  6*s*(400*pow(x,3) - 1980*pow(x,2)*z + 3240*x*z2 - 2520*z3))) +
				    3*r*(-1350*s4 - 270*pow(x,4) + 2700*pow(x,3)*z - 9180*pow(x,2)*z2 + 15120*x*z3 - 17850*z4 +
					 45*s3*(-60*x + 240*z) + 27*s2*(-100*pow(x,2) + 600*x*z - 1020*z2) +
					 9*s*(-150*pow(x,3) + 1200*pow(x,2)*z - 3060*x*z2 + 3360*z3)) -
				    4*c*z*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 - 1890*x*z4 +
					   1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					   6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					   6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4)) -
				    2*c*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
					 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					 6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					 6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		     12*l4*l5*r*z2*(15*r*(840*s4 + 168*pow(x,4) + 14*s3*(120*x - 288*z) + 336*s*pow(x - z,3) - 1008*pow(x,3)*z +
						   2184*pow(x,2)*z2 - 1680*x*z3 + 1470*z4 - 252*s*pow(x - z,2)*(-x + 8*z) +
						   252*s*(x - z)*(pow(x,2) - x*z + 4*z2) + 42*s2*(40*pow(x,2) - 144*x*z + 156*z2)) +
					     c*z2*(2940*s4 + 588*pow(x,4) + 70*s3*(84*x - 144*z) - 2520*pow(x,3)*z + 3780*pow(x,2)*z2 - 2100*x*z3 +
							 1470*z4 + 105*s2*(56*pow(x,2) - 144*x*z + 108*z2) +
							 7*s*(420*pow(x,3) - 1440*pow(x,2)*z + 1620*x*z2 - 600*z3)) +
					     4*c*z*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						    294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						    105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						    7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					     2*c*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
						  315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     48*l4*l5*r*z*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 - 420*x*z4 +
					    294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) + 14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) -
					    126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) + 42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
				      c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						  294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
				      2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
					     315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     24*l4*l5*r*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					  168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					  42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					  42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				    c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*(-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) -
								    420*s2*x*(x - z)*(2*x - z) + 70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) +
								    140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
	      63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
			    4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
	      945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
	      14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
				    105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
				    420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
				    840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
				    15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
	      30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
			    70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
			    84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
			    7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
	      98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				    70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				    210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
	      60*l6*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
				    7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
				    28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
				    7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
				    2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
				    2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
				       42*z6)) + 180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
									      7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									      7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									      28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									      7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										   84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
												       42*x*z5 + 42*z6)) + 15*l4*z5*
	      (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
	       14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
	       84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
	       4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
	       7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
	       4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
	      5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
			    14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
			    84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
			    21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
			    8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
			    8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
			       63*z6)) - 5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
							   420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
							   14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
							   84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
							   21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
								 504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
										      378*x*z5 + 63*z6)) + 75*l4*z4*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) - 36*l7*l8*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
	      28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						  105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
				      15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
					    7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
				      2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
					     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
					     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
	      168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						   105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						   35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
				       15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
					     70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
					     pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
			      c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				 15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				 3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				 s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
	      12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						  140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						  pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
				      6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					   35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					   pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
				      2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					     14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
					     140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
	      48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
	      8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					   1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					   6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					   6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
			     3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				  45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				  9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
			     2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
				    540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
				    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
				    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
				    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
	      12*l4*l5*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
					    420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
					    14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
					    42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
				      c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						  294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
				      2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
					     315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      24*l4*l5*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
				     168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
				     42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
				     42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
			       c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					   420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
					   70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					   105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					   7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (6.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*z*(36*c*l*l8*l8*r*(-6*s - 3*x + 6*z) + 105*s*z4*z5*
		  (-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*s2*x*(x - z)*(2*x - z) +
		   70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) + 140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		  63*l6*l6*(-60*s3 + 60*s2*(-x - 2*z) + 5*s*(-4*pow(x,2) - 36*x*z + 84*z2) +
				4*(-20*pow(x,2)*z + 60*x*z2 - 60*z3)) +
		  945*s*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				  70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				  140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		  14*l2*z7*(-8820*s6 - 90*pow(x,6) + 560*pow(x,5)*z - 756*pow(x,4)*z2 + 350*pow(x,2)*z4 +
					105*s5*(-213*x + 232*z) + 210*s4*(-145*pow(x,2) + 240*x*z - 126*z2) +
					420*s2*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) -
					840*s2*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
					15*s*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		  30*l5*l5*(-420*s5 - 168*pow(x,4)*z + 840*pow(x,3)*z2 - 1960*pow(x,2)*z3 + 2520*x*z4 - 1680*z5 +
				70*s4*(-12*x + 12*z) + 35*s3*(-24*pow(x,2) + 24*x*z + 60*z2) +
				84*s2*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
				7*s*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		  98*l2*z6*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		  60*l6*z3*(168*s6 + 7*s5*(87*x - 576*z) + 70*s4*(17*pow(x,2) - 132*x*z + 168*z2) +
					7*s3*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
					28*s2*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
					7*s*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5) -
					2*z*(98*pow(x,5) - 336*pow(x,4)*z + 630*pow(x,3)*z2 - 700*pow(x,2)*z3 + 210*x*z4 + 252*z5) -
					2*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 + 42*x*z5 +
					   42*z6)) + 180*l6*z2*(196*s7 + 84*s6*(7*x + 2*z) +
										  7*s5*(140*pow(x,2) + 87*x*z - 288*z2) + 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
										  7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
										  28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
										  7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										       84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
													   42*x*z5 + 42*z6)) + 15*l4*z5*
		  (-6888*s6 + 84*s5*(-215*x + 372*z) + 280*s4*(-92*pow(x,2) + 240*x*z - 198*z2) +
		   14*s3*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
		   84*s2*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
		   4*z*(196*pow(x,5) - 462*pow(x,4)*z + 420*pow(x,3)*z2 - 280*pow(x,2)*z3 + 294*z5) +
		   7*s*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5) +
		   4*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6)) -
		  5*l8*z*(-3276*s6 + 42*s5*(-183*x - 4*z) + 420*s4*(-22*pow(x,2) - 12*x*z + 81*z2) +
				14*s3*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				84*s2*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				21*s*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5) +
				8*z*(-98*pow(x,5) + 630*pow(x,4)*z - 1890*pow(x,3)*z2 + 2800*pow(x,2)*z3 - 1890*x*z4 + 378*z5) +
				8*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 - 378*x*z5 +
				   63*z6)) - 5*l8*(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
							       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
							       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
							       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
							       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
								     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
											  378*x*z5 + 63*z6)) + 75*l4*z4*
		  (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		   280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		   14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		   84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		   4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		   7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
			84*z6)) - 36*l7*l8*r*(60*r*z + c*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								    10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2))) +
		  28*l5*r*z6*(c*z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						      105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					  15*r*(-168*s5 - 420*s2*pow(x,2)*(x - z) + 35*s4*(-12*x + 6*z) + 70*s3*x*(-8*x + 6*z) +
						7*s*pow(x,3)*(-24*x + 30*z) + pow(x,4)*(-28*x + 42*z)) +
					  2*c*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
						 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
						 pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) +
		  168*l5*r*z5*(c*z2*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
						       105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
						       35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
					   15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
						 70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
						 pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		  12*l6*l7*r*(3*r*(-540*s2*z - 180*pow(x,2)*z + 540*x*z2 - 800*z3 + 45*s*(-12*x*z + 24*z2)) -
				  c*(-42*s5 - 7*pow(x,5) + 66*pow(x,4)*z - 300*pow(x,3)*z2 + 680*pow(x,2)*z3 - 990*x*z4 + 852*z5 +
				     15*s4*(-7*x + 22*z) + 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) +
				     3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) +
				     s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4))) -
		  12*l7*r*z4*(c*z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						      140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						      pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					  6*r*(-1176*s5 + 105*s2*pow(x - z,2)*(-6*x + 14*z) + 35*s4*(-84*x + 90*z) +
					       35*s3*(-112*pow(x,2) + 180*x*z - 120*z2) - 210*s2*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       7*s*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3) +
					       pow(x,2)*(-196*pow(x,3) + 630*pow(x,2)*z - 1050*x*z2 + 980*z3)) +
					  2*c*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						 14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						 42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						 pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) -
		  48*l7*r*z3*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						      140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						      pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					  6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					       pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) +
		  8*l5*l6*r*(-(c*z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					       1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					       6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					       6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				 3*r*(-1350*s4*z - 270*pow(x,4)*z + 1350*pow(x,3)*z2 - 3060*pow(x,2)*z3 + 3780*x*z4 - 3570*z5 +
				      45*s3*(-60*x*z + 120*z2) + 27*s2*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
				      9*s*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				 2*c*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
					540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		  12*l4*l5*r*z2*(15*r*(-168*s5 - 28*pow(x,5) + 168*pow(x,4)*z - 504*pow(x,3)*z2 + 728*pow(x,2)*z3 -
						420*x*z4 + 294*z5 + 210*s4*(-2*x + 4*z) + 42*s*pow(x - z,3)*(-x + 8*z) +
						14*s3*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*s*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
						42*s2*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					  c*z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 - 525*x*z4 +
						      294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						      105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						      7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					  2*c*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
						 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
						 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						 105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						 7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  24*l4*l5*r*z*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
					       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (x*pow(-126*l4*z - 84*l2*z3 - 1470*z5,2)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (6.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,3)) -
	    (x*(-126*l4 - 252*l2*z2 - 7350*z4)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) -
	    (x*z*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
	      105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			      70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			      140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			    5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
	      - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
				      70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
				      210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				      420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				      35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				      15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			    490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			    35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
	      60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				    70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				    7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				    28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				    7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					 84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							     42*x*z5 + 42*z6)) - 5*l8*z*
	      (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
	       420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
	       14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
	       84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
	       21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		     504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					  378*x*z5 + 63*z6)) + 15*l4*z5*
	      (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
	       280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
	       14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
	       84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
	       4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
	       7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		    84*z6)) + 28*l5*r*z6*(c*z2*
							    (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							    15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								  70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								  pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
	      36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						 20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
	      12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
				      6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					   pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
	      12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				   90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				 170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				 5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				 s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
	      12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					    168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					    42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						  70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					   495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					   30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				  756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				  45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l5*pow(l2 - z2,4)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (4*x*z2*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
			   105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
					   70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
					   140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
			   63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
					 5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
			   - 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
						   70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
						   210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
						   420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
						   35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
						   15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
			   30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
					 490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
					 35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
					 84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
					 7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
			   60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
						 70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
						 7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
						 28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
						 7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
						      84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
									  42*x*z5 + 42*z6)) - 5*l8*z*
			   (1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
			    420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				  504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
						       378*x*z5 + 63*z6)) + 15*l4*z5*
			   (1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
			    280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
			    14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
			    84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
			    4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
			    7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
				 84*z6)) + 28*l5*r*z6*(c*z2*
									 (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
									  105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
									  pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
									 15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
									       70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
									       pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
			   36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
							      20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
							      10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
			   12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
							       14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
							       140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							       42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
							       pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
						   6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
							35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
							7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
							pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
			   12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
						90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
					   c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
					      170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
					      5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
					      3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
					      s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
			   12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
							 168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
							 42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
							 42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
						   c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
							       420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
							       70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							       105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
							       7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
			   8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
							495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
							30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
							6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
							6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
					  3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
					       756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
					       45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
					       27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
					       9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (l5*pow(l2 - z2,5)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (x*(36*l*l8*l8*r*(3*r + c*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2)) +
		105*s*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		63*l6*l6*(40*s4 + 10*s3*(7*x - 6*z) + 60*s2*(pow(x,2) - x*z - z2) +
			      5*s*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3) + 4*(pow(x,4) - 10*pow(x,2)*z2 + 20*x*z3 - 15*z4))
		- 14*l2*z7*(2520*s7 + 1260*s6*(6*x - 7*z) - 90*pow(x,6)*z + 280*pow(x,5)*z2 - 252*pow(x,4)*z3 +
					70*pow(x,2)*z5 + 105*s5*(120*pow(x,2) - 213*x*z + 116*z2) +
					210*s4*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					420*s2*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					35*s3*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*s*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		30*l5*l5*(84*s6 + 6*pow(x,6) + 21*s5*(11*x - 20*z) - 84*pow(x,4)*z2 + 280*pow(x,3)*z3 -
			      490*pow(x,2)*z4 + 504*x*z5 - 280*z6 + 70*s4*(5*pow(x,2) - 12*x*z + 6*z2) +
			      35*s3*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			      84*s2*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			      7*s*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) +
		60*l6*z3*(196*s7 + 84*s6*(7*x + 2*z) + 7*s5*(140*pow(x,2) + 87*x*z - 288*z2) +
				      70*s4*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
				      7*s3*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
				      28*s2*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
				      7*s*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
					   84*z6) - 2*z*(-27*pow(x,6) + 98*pow(x,5)*z - 168*pow(x,4)*z2 + 210*pow(x,3)*z3 - 175*pow(x,2)*z4 +
							       42*x*z5 + 42*z6)) - 5*l8*z*
		(1764*s7 + 252*s6*(21*x - 13*z) + 42*s5*(210*pow(x,2) - 183*x*z - 2*z2) +
		 420*s4*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
		 14*s3*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
		 84*s2*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
		 21*s*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
		       504*z6) + 8*z*(18*pow(x,6) - 98*pow(x,5)*z + 315*pow(x,4)*z2 - 630*pow(x,3)*z3 + 700*pow(x,2)*z4 -
					    378*x*z5 + 63*z6)) + 15*l4*z5*
		(1176*s7 + 168*s6*(21*x - 41*z) + 84*s5*(70*pow(x,2) - 215*x*z + 186*z2) +
		 280*s4*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
		 14*s3*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
		 84*s2*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
		 4*z*(-60*pow(x,6) + 196*pow(x,5)*z - 231*pow(x,4)*z2 + 140*pow(x,3)*z3 - 70*pow(x,2)*z4 + 49*z6) +
		 7*s*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
		      84*z6)) + 28*l5*r*z6*(c*z2*
							      (105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
							       105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
							       pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) +
							      15*r*(70*s6 + 42*s5*(5*x - 4*z) + 210*s2*pow(x,2)*pow(x - z,2) + 35*s4*(10*pow(x,2) - 12*x*z + 3*z2) +
								    70*s3*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*s*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2) +
								    pow(x,4)*(10*pow(x,2) - 28*x*z + 21*z2))) -
		36*l7*l8*r*(30*r*z2 + c*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 -
						   20*x*z3 + 20*z4 + 10*s2*(2*pow(x,2) - 3*x*z + 3*z2) +
						   10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))) -
		12*l7*r*z4*(c*z2*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
						    140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
						    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
					6*r*(385*s6 + 21*s5*(55*x - 56*z) + 105*s2*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					     35*s4*(55*pow(x,2) - 84*x*z + 45*z2) + 35*s3*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					     7*s*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4) +
					     pow(x,2)*(55*pow(x,4) - 196*pow(x,3)*z + 315*pow(x,2)*z2 - 350*x*z3 + 245*z4))) -
		12*l6*l7*r*(3*r*(45*s4 + 90*s3*x + 9*pow(x,4) - 90*pow(x,2)*z2 + 180*x*z3 - 200*z4 +
				     90*s2*(pow(x,2) - 3*z2) + 45*s*(pow(x,3) - 6*x*z2 + 8*z3)) -
				c*(21*s6 + 3*pow(x,6) + 21*s5*(3*x - 2*z) - 7*pow(x,5)*z + 33*pow(x,4)*z2 - 100*pow(x,3)*z3 +
				   170*pow(x,2)*z4 - 198*x*z5 + 142*z6 + 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) +
				   5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				   3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) +
				   s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5))) +
		12*l4*l5*r*z2*(15*r*(42*s6 + 6*pow(x,6) + 42*s5*(3*x - 4*z) - 28*pow(x,5)*z + 84*pow(x,4)*z2 -
					      168*pow(x,3)*z3 + 182*pow(x,2)*z4 - 84*x*z5 + 49*z6 + 210*s4*(pow(x,2) - 2*x*z + 2*z2) +
					      42*s*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 14*s3*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					      42*s2*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					c*z2*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
						    420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 +
						    70*s4*(25*pow(x,2) - 33*x*z + 21*z2) + 70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						    105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						    7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		8*l5*l6*r*(-(c*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					     495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 +
					     30*s4*(21*pow(x,2) - 35*x*z + 40*z2) + 10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					     6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					     6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			       3*r*(63*s6 + 189*s5*x + 9*pow(x,6) - 135*pow(x,4)*z2 + 450*pow(x,3)*z3 - 765*pow(x,2)*z4 +
				    756*x*z5 - 595*z6 + 45*s4*(7*pow(x,2) - 15*z2) +
				    45*s3*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				    27*s2*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				    9*s*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='z'&&label_j=='c'){
	result = -(x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
		   (36*l*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		    28*l5*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					    105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					    35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		    36*l7*l8*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				    10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		    12*l7*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					    14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					    42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					    pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		    12*l6*l7*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				    170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				    5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				    3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				    s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		    8*l5*l6*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					    495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		    12*l4*l5*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					    420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					    70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					    105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					    7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*z*(36*l*l8*l8*r*(3*s2 + 3*s*x + pow(x,2) - 6*s*z - 3*x*z + 3*z2) +
		  28*l5*r*z8*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					  105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					  35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		  36*l7*l8*r*(10*s4 + 2*pow(x,4) + 20*s3*(x - z) - 5*pow(x,3)*z + 10*pow(x,2)*z2 - 20*x*z3 + 20*z4 +
				  10*s2*(2*pow(x,2) - 3*x*z + 3*z2) + 10*s*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		  12*l7*r*z6*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					  14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					  42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					  pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) -
		  12*l6*l7*r*(-21*s6 - 3*pow(x,6) - 21*s5*(3*x - 2*z) + 7*pow(x,5)*z - 33*pow(x,4)*z2 + 100*pow(x,3)*z3 -
				  170*pow(x,2)*z4 + 198*x*z5 - 142*z6 - 15*s4*(7*pow(x,2) - 7*x*z + 11*z2) -
				  5*s3*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				  3*s2*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4) -
				  s*(21*pow(x,5) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 + 510*x*z4 - 396*z5)) -
		  8*l5*l6*r*z2*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 -
					  495*pow(x,3)*z3 + 540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
					  10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					  6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					  6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		  12*l4*l5*r*z4*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					  420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					  70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					  105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					  7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) +
	    (x*(36*l*l8*l8*r*(-6*s - 3*x + 6*z) - 36*l7*l8*r*(-20*s3 - 5*pow(x,3) + 20*pow(x,2)*z - 60*x*z2 + 80*z3 +
								    10*s2*(-3*x + 6*z) + 10*s*(-2*pow(x,2) + 6*x*z - 12*z2)) +
		168*l5*r*z7*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) +
					 105*s4*(5*pow(x,2) - 5*x*z + z2) + 105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
					 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) + pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2)) -
		12*l6*l7*r*(42*s5 + 7*pow(x,5) - 66*pow(x,4)*z + 300*pow(x,3)*z2 - 680*pow(x,2)*z3 + 990*x*z4 -
				852*z5 - 15*s4*(-7*x + 22*z) - 5*s3*(-28*pow(x,2) + 132*x*z - 240*z2) -
				3*s2*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3) -
				s*(-42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4)) -
		48*l7*r*z5*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4)) +
		24*l4*l5*r*z3*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 -
					420*pow(x,3)*z3 + 315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)) +
		28*l5*r*z6*(z2*(-210*s5 - 210*s*pow(x,3)*(x - z) + 105*s4*(-5*x + 2*z) +
						  105*s2*pow(x,2)*(-5*x + 4*z) + 35*s3*x*(-20*x + 12*z) + pow(x,4)*(-35*x + 42*z)) +
					2*z*(105*s6 + 105*s5*(3*x - 2*z) + 105*s*pow(x,3)*pow(x - z,2) + 105*s4*(5*pow(x,2) - 5*x*z + z2) +
					     105*s2*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 35*s3*x*(15*pow(x,2) - 20*x*z + 6*z2) +
					     pow(x,4)*(15*pow(x,2) - 35*x*z + 21*z2))) -
		12*l7*r*z4*(z2*(-1176*s5 + 70*s4*(-42*x + 30*z) + 14*s*x*pow(x - z,2)*(-8*x + 42*z) +
						  140*s3*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*s*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
						  42*s2*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3) +
						  pow(x,2)*(-196*pow(x,3) + 420*pow(x,2)*z - 525*x*z2 + 392*z3)) +
					2*z*(532*s6 + 84*s5*(19*x - 14*z) + 70*s4*(38*pow(x,2) - 42*x*z + 15*z2) +
					     14*s*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) +
					     140*s3*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					     42*s2*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4) +
					     pow(x,2)*(76*pow(x,4) - 196*pow(x,3)*z + 210*pow(x,2)*z2 - 175*x*z3 + 98*z4))) +
		8*l5*l6*r*(-(z2*(-420*s5 - 70*pow(x,5) + 480*pow(x,4)*z - 1485*pow(x,3)*z2 + 2160*pow(x,2)*z3 -
					   1890*x*z4 + 1260*z5 + 30*s4*(-35*x + 80*z) + 10*s3*(-140*pow(x,2) + 480*x*z - 594*z2) +
					   6*s2*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
					   6*s*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) -
			       2*z*(126*s6 + 18*pow(x,6) + 42*s5*(9*x - 10*z) - 70*pow(x,5)*z + 240*pow(x,4)*z2 - 495*pow(x,3)*z3 +
				    540*pow(x,2)*z4 - 378*x*z5 + 210*z6 + 30*s4*(21*pow(x,2) - 35*x*z + 40*z2) +
				    10*s3*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
				    6*s2*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
				    6*s*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		12*l4*l5*r*z2*(z2*(-924*s5 - 154*pow(x,5) + 588*pow(x,4)*z - 1260*pow(x,3)*z2 + 1260*pow(x,2)*z3 -
						  525*x*z4 + 294*z5 + 70*s4*(-33*x + 42*z) + 70*s3*(-44*pow(x,2) + 84*x*z - 72*z2) +
						  105*s2*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
						  7*s*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					2*z*(350*s6 + 50*pow(x,6) + 42*s5*(25*x - 22*z) - 154*pow(x,5)*z + 294*pow(x,4)*z2 - 420*pow(x,3)*z3 +
					     315*pow(x,2)*z4 - 105*x*z5 + 49*z6 + 70*s4*(25*pow(x,2) - 33*x*z + 21*z2) +
					     70*s3*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     105*s2*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*s*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='z'&&label_j=='s'){
	result = (x*(-216*c*l*l8*l8*r - 36*c*l7*l8*r*(-60*s2 + 20*s*(-3*x + 6*z) + 10*(-2*pow(x,2) + 6*x*z - 12*z2)) +
		     105*s*z4*z5*(-2100*s4 - 420*s*x*pow(x - z,2) - 420*pow(x,2)*pow(x - z,2) - 840*s*x*(x - z)*(2*x - z) +
				     280*s3*(-15*x + 12*z) + 420*s2*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		     105*z4*z5*(-420*s5 - 210*s2*x*pow(x - z,2) - 420*s*pow(x,2)*pow(x - z,2) - 420*s2*x*(x - z)*(2*x - z) +
				   70*s4*(-15*x + 12*z) + pow(x,3)*(-70*pow(x,2) + 168*x*z - 105*z2) + 140*s3*(-10*pow(x,2) + 12*x*z - 3*z2)) +
		     63*l6*l6*(-180*s2 + 120*s*(-x - 2*z) + 5*(-4*pow(x,2) - 36*x*z + 84*z2)) +
		     945*s*z8*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				     280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		     945*z8*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				   70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				   140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		     14*l2*z7*(-52920*s5 + 525*s4*(-213*x + 232*z) + 840*s3*(-145*pow(x,2) + 240*x*z - 126*z2) +
					   840*s*pow(x - z,2)*(-12*pow(x,2) + 10*x*z - 6*z2) - 1680*s*(x - z)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					   105*s2*(-675*pow(x,3) + 1520*pow(x,2)*z - 1188*x*z2 + 432*z3) +
					   15*x*(-133*pow(x,4) + 504*pow(x,3)*z - 630*pow(x,2)*z2 + 336*x*z3 - 105*z4)) -
		     30*l5*l5*(-2100*s4 + 280*s3*(-12*x + 12*z) + 105*s2*(-24*pow(x,2) + 24*x*z + 60*z2) +
				   168*s*(-5*pow(x,3) + 45*x*z2 - 80*z3) +
				   7*(-12*pow(x,4) - 60*pow(x,3)*z + 420*pow(x,2)*z2 - 900*x*z3 + 840*z4)) -
		     5*l8*z*(-19656*s5 + 210*s4*(-183*x - 4*z) + 1680*s3*(-22*pow(x,2) - 12*x*z + 81*z2) +
				   42*s2*(-405*pow(x,3) - 920*pow(x,2)*z + 4995*x*z2 - 7920*z3) +
				   168*s*(-15*pow(x,4) - 170*pow(x,3)*z + 855*pow(x,2)*z2 - 1800*x*z3 + 1670*z4) +
				   21*(14*pow(x,5) - 360*pow(x,4)*z + 1755*pow(x,3)*z2 - 4320*pow(x,2)*z3 + 5340*x*z4 - 3024*z5)) +
		     60*l6*z3*(1008*s5 + 35*s4*(87*x - 576*z) + 280*s3*(17*pow(x,2) - 132*x*z + 168*z2) +
					   21*s2*(195*pow(x,3) - 1600*pow(x,2)*z + 2880*x*z2 - 2400*z3) +
					   56*s*(33*pow(x,4) - 270*pow(x,3)*z + 600*pow(x,2)*z2 - 720*x*z3 + 465*z4) +
					   7*(49*pow(x,5) - 384*pow(x,4)*z + 960*pow(x,3)*z2 - 1440*pow(x,2)*z3 + 1305*x*z4 - 504*z5)) -
		     98*l2*z6*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
					   840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					   840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					   105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					   15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) +
		     15*l4*z5*(-41328*s5 + 420*s4*(-215*x + 372*z) + 1120*s3*(-92*pow(x,2) + 240*x*z - 198*z2) +
					   42*s2*(-1530*pow(x,3) + 5360*pow(x,2)*z - 6435*x*z2 + 3640*z3) +
					   168*s*(-122*pow(x,4) + 540*pow(x,3)*z - 825*pow(x,2)*z2 + 660*x*z3 - 295*z4) +
					   7*(-364*pow(x,5) + 1968*pow(x,4)*z - 3465*pow(x,3)*z2 + 3280*pow(x,2)*z3 - 2070*x*z4 + 504*z5)) -
		     5*l8*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
				 1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				 42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				 168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				 21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				     504*z6)) + 180*l6*z2*(1372*s6 + 504*s5*(7*x + 2*z) +
									     35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									     21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									     56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									     7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
										84*z6)) + 75*l4*z4*(8232*s6 + 1008*s5*(21*x - 41*z) +
														      420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														      42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														      168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														      7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
															 84*z6)) + 28*l5*r*z6*(c*z2*
																				 (-1050*s4 - 210*pow(x,3)*(x - z) + 420*s3*(-5*x + 2*z) + 210*s*pow(x,2)*(-5*x + 4*z) + 105*s2*x*(-20*x + 12*z)) +
																				 15*r*(-840*s4 - 840*s*pow(x,2)*(x - z) + 140*s3*(-12*x + 6*z) + 210*s2*x*(-8*x + 6*z) + 7*pow(x,3)*(-24*x + 30*z)) +
																				 2*c*z*(630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																					210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2))) +
		     168*l5*r*z5*(c*z2*(630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) +
							  420*s3*(5*pow(x,2) - 5*x*z + z2) + 210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) +
							  105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
					      15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
						    210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		     12*l6*l7*r*(3*r*(-1080*s*z + 45*(-12*x*z + 24*z2)) -
				     c*(-210*s4 - 42*pow(x,4) + 330*pow(x,3)*z - 1200*pow(x,2)*z2 + 2040*x*z3 - 1980*z4 +
					60*s3*(-7*x + 22*z) + 15*s2*(-28*pow(x,2) + 132*x*z - 240*z2) +
					6*s*(-35*pow(x,3) + 220*pow(x,2)*z - 600*x*z2 + 680*z3))) -
		     12*l7*r*z4*(c*z2*(-5880*s4 + 280*s3*(-42*x + 30*z) + 14*x*pow(x - z,2)*(-8*x + 42*z) +
							 420*s2*(-28*pow(x,2) + 30*x*z - 15*z2) - 28*x*(x - z)*(38*pow(x,2) - 8*x*z + 21*z2) +
							 84*s*(-70*pow(x,3) + 100*pow(x,2)*z - 75*x*z2 + 28*z3)) +
					     6*r*(-5880*s4 + 210*s*pow(x - z,2)*(-6*x + 14*z) + 140*s3*(-84*x + 90*z) +
						  105*s2*(-112*pow(x,2) + 180*x*z - 120*z2) - 420*s*(x - z)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  7*x*(-168*pow(x,3) + 450*pow(x,2)*z - 600*x*z2 + 420*z3)) +
					     2*c*z*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						    14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						    84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4))) -
		     48*l7*r*z3*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
							 14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							 84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					     6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
						  7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) +
		     8*l5*l6*r*(-(c*z2*(-2100*s4 + 120*s3*(-35*x + 80*z) + 30*s2*(-140*pow(x,2) + 480*x*z - 594*z2) +
						  12*s*(-175*pow(x,3) + 800*pow(x,2)*z - 1485*x*z2 + 1080*z3) +
						  6*(-70*pow(x,4) + 400*pow(x,3)*z - 990*pow(x,2)*z2 + 1080*x*z3 - 630*z4))) +
				    3*r*(-5400*s3*z + 135*s2*(-60*x*z + 120*z2) + 54*s*(-100*pow(x,2)*z + 300*x*z2 - 340*z3) +
					 9*(-150*pow(x,3)*z + 600*pow(x,2)*z2 - 1020*x*z3 + 840*z4)) -
				    2*c*z*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					   30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
		     12*l4*l5*r*z2*(15*r*(-840*s4 + 840*s3*(-2*x + 4*z) + 42*pow(x - z,3)*(-x + 8*z) +
						   42*s2*(-40*pow(x,2) + 120*x*z - 144*z2) - 126*pow(x - z,2)*(pow(x,2) - x*z + 4*z2) +
						   84*s*(-10*pow(x,3) + 40*pow(x,2)*z - 72*x*z2 + 52*z3)) +
					     c*z2*(-4620*s4 + 280*s3*(-33*x + 42*z) + 210*s2*(-44*pow(x,2) + 84*x*z - 72*z2) +
							 210*s*(-22*pow(x,3) + 56*pow(x,2)*z - 72*x*z2 + 36*z3) +
							 7*(-132*pow(x,4) + 420*pow(x,3)*z - 720*pow(x,2)*z2 + 540*x*z3 - 150*z4)) +
					     2*c*z*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						    210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						    210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						    7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		     24*l4*l5*r*z*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
					    42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						  210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6)) -
	    (x*(-126*l4*z - 84*l2*z3 - 1470*z5)*
	     (36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
	      (40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) + 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))
	      + 105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
			    70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
			    140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
	      63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
			    5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
	      14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
				    840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
				    840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
				    105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
				    15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
	      30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
			    105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
			    168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
			    7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
	      5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
			    1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
			    42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
			    168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
			    21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				504*z6)) + 60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
								       35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
								       21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
								       56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
								       7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									  84*z6)) + 15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
														420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
														   84*z6)) + 28*l5*r*z6*(c*z2*
																			   (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																			    210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																			   15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																				 210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
	      12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						  14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						  84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
				      6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					   140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					   7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
	      12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
			      c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				 510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				 15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				 6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
	      12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
					    42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
					    84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
				      c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						  210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						  210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						  7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
	      8*l5*l6*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					   30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					   12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					   6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
			     3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				  135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				  54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				  9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (12.*l5*pow(l2 - z2,3)*pow(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6,2)) +
	    (x*z*(36*c*l*l8*l8*r*(6*s + 3*x - 6*z) - 36*c*l7*l8*r*
		  (40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) + 10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3))
		  + 105*s*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				    280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  105*z4*z5*(140*s6 + 420*s5*(x - z) + 140*s*pow(x,2)*pow(x - z,3) + 210*s2*x*pow(x - z,2)*(2*x - z) +
				70*s4*(10*pow(x,2) - 15*x*z + 6*z2) + pow(x,3)*(20*pow(x,3) - 70*pow(x,2)*z + 84*x*z2 - 35*z3) +
				140*s3*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		  63*l6*l6*(160*s3 + 30*s2*(7*x - 6*z) + 120*s*(pow(x,2) - x*z - z2) +
				5*(5*pow(x,3) - 4*pow(x,2)*z - 18*x*z2 + 28*z3)) -
		  14*l2*z7*(17640*s6 + 7560*s5*(6*x - 7*z) + 525*s4*(120*pow(x,2) - 213*x*z + 116*z2) +
					840*s3*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					840*s*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					105*s2*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4) +
					15*x*(24*pow(x,5) - 133*pow(x,4)*z + 252*pow(x,3)*z2 - 210*pow(x,2)*z3 + 84*x*z4 - 21*z5)) -
		  30*l5*l5*(504*s5 + 105*s4*(11*x - 20*z) + 280*s3*(5*pow(x,2) - 12*x*z + 6*z2) +
				105*s2*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) +
				168*s*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4) +
				7*(7*pow(x,5) - 12*pow(x,4)*z - 30*pow(x,3)*z2 + 140*pow(x,2)*z3 - 225*x*z4 + 168*z5)) -
		  5*l8*z*(12348*s6 + 1512*s5*(21*x - 13*z) + 210*s4*(210*pow(x,2) - 183*x*z - 2*z2) +
				1680*s3*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				42*s2*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				168*s*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5) +
				21*(12*pow(x,6) + 14*pow(x,5)*z - 180*pow(x,4)*z2 + 585*pow(x,3)*z3 - 1080*pow(x,2)*z4 + 1068*x*z5 -
				    504*z6)) + 60*l6*z3*(1372*s6 + 504*s5*(7*x + 2*z) +
									   35*s4*(140*pow(x,2) + 87*x*z - 288*z2) + 280*s3*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
									   21*s2*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
									   56*s*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5) +
									   7*(4*pow(x,6) + 49*pow(x,5)*z - 192*pow(x,4)*z2 + 320*pow(x,3)*z3 - 360*pow(x,2)*z4 + 261*x*z5 -
									      84*z6)) + 15*l4*z5*(8232*s6 + 1008*s5*(21*x - 41*z) +
														    420*s4*(70*pow(x,2) - 215*x*z + 186*z2) + 1120*s3*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
														    42*s2*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
														    168*s*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5) +
														    7*(24*pow(x,6) - 364*pow(x,5)*z + 984*pow(x,4)*z2 - 1155*pow(x,3)*z3 + 820*pow(x,2)*z4 - 414*x*z5 +
														       84*z6)) + 28*l5*r*z6*(c*z2*
																			       (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
																				210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) +
																			       15*r*(420*s5 + 210*s4*(5*x - 4*z) + 420*s*pow(x,2)*pow(x - z,2) + 140*s3*(10*pow(x,2) - 12*x*z + 3*z2) +
																				     210*s2*x*(5*pow(x,2) - 8*x*z + 3*z2) + 7*pow(x,3)*(10*pow(x,2) - 24*x*z + 15*z2))) -
		  12*l7*r*z4*(c*z2*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
						      14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
						      84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) +
					  6*r*(2310*s5 + 105*s4*(55*x - 56*z) + 210*s*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
					       140*s3*(55*pow(x,2) - 84*x*z + 45*z2) + 105*s2*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3) +
					       7*x*(55*pow(x,4) - 168*pow(x,3)*z + 225*pow(x,2)*z2 - 200*x*z3 + 105*z4))) -
		  12*l6*l7*r*(3*r*(180*s3 + 270*s2*x + 180*s*(pow(x,2) - 3*z2) + 45*(pow(x,3) - 6*x*z2 + 8*z3)) -
				  c*(126*s5 + 21*pow(x,5) + 105*s4*(3*x - 2*z) - 42*pow(x,4)*z + 165*pow(x,3)*z2 - 400*pow(x,2)*z3 +
				     510*x*z4 - 396*z5 + 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) +
				     15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
				     6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		  12*l4*l5*r*z2*(15*r*(252*s5 + 210*s4*(3*x - 4*z) + 840*s3*(pow(x,2) - 2*x*z + 2*z2) +
						42*pow(x - z,3)*(pow(x,2) - x*z + 4*z2) + 42*s2*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						84*s*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4)) +
					  c*z2*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
						      210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
						      210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
						      7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))) +
		  8*l5*l6*r*(-(c*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					       30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					       12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					       6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5))) +
				 3*r*(378*s5 + 945*s4*x + 180*s3*(7*pow(x,2) - 15*z2) +
				      135*s2*(7*pow(x,3) - 30*x*z2 + 40*z3) +
				      54*s*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4) +
				      9*(7*pow(x,5) - 75*pow(x,3)*z2 + 200*pow(x,2)*z3 - 255*x*z4 + 168*z5)))))/
	    (2.*l5*pow(l2 - z2,4)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='c'&&label_j=='c'){
	result = 0;
    }
    else if (label_i=='c'&&label_j=='s'){
	result = (x*(36*l*l8*l8*r*(6*s + 3*x - 6*z) + 28*l5*r*z8*
		     (630*s5 + 525*s4*(3*x - 2*z) + 105*pow(x,3)*pow(x - z,2) + 420*s3*(5*pow(x,2) - 5*x*z + z2) +
		      210*s*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 105*s2*x*(15*pow(x,2) - 20*x*z + 6*z2)) -
		     36*l7*l8*r*(40*s3 + 60*s2*(x - z) + 20*s*(2*pow(x,2) - 3*x*z + 3*z2) +
				     10*(pow(x,3) - 2*pow(x,2)*z + 3*x*z2 - 4*z3)) -
		     12*l7*r*z6*(3192*s5 + 420*s4*(19*x - 14*z) + 280*s3*(38*pow(x,2) - 42*x*z + 15*z2) +
					     14*x*pow(x - z,2)*(38*pow(x,2) - 8*x*z + 21*z2) + 420*s2*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
					     84*s*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4)) -
		     12*l6*l7*r*(-126*s5 - 21*pow(x,5) - 105*s4*(3*x - 2*z) + 42*pow(x,4)*z - 165*pow(x,3)*z2 +
				     400*pow(x,2)*z3 - 510*x*z4 + 396*z5 - 60*s3*(7*pow(x,2) - 7*x*z + 11*z2) -
				     15*s2*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) -
				     6*s*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4)) -
		     8*l5*l6*r*z2*(756*s5 + 210*s4*(9*x - 10*z) + 120*s3*(21*pow(x,2) - 35*x*z + 40*z2) +
					     30*s2*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
					     12*s*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4) +
					     6*(21*pow(x,5) - 70*pow(x,4)*z + 200*pow(x,3)*z2 - 330*pow(x,2)*z3 + 270*x*z4 - 126*z5)) +
		     12*l4*l5*r*z4*(2100*s5 + 210*s4*(25*x - 22*z) + 280*s3*(25*pow(x,2) - 33*x*z + 21*z2) +
					     210*s2*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
					     210*s*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4) +
					     7*(50*pow(x,5) - 132*pow(x,4)*z + 210*pow(x,3)*z2 - 240*pow(x,2)*z3 + 135*x*z4 - 30*z5))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else if (label_i=='s'&&label_j=='s'){
	result = (x*(216*c*l*l8*l8*r + 63*l6*l6*(480*s2 + 60*s*(7*x - 6*z) + 120*(pow(x,2) - x*z - z2)) -
		     36*c*l7*l8*r*(120*s2 + 120*s*(x - z) + 20*(2*pow(x,2) - 3*x*z + 3*z2)) +
		     105*s*z4*z5*(4200*s4 + 8400*s3*(x - z) + 420*x*pow(x - z,2)*(2*x - z) + 840*s2*(10*pow(x,2) - 15*x*z + 6*z2) +
				     840*s*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) +
		     210*z4*z5*(840*s5 + 2100*s4*(x - z) + 140*pow(x,2)*pow(x - z,3) + 420*s*x*pow(x - z,2)*(2*x - z) +
				   280*s3*(10*pow(x,2) - 15*x*z + 6*z2) + 420*s2*(5*pow(x,3) - 10*pow(x,2)*z + 6*x*z2 - z3)) -
		     30*l5*l5*(2520*s4 + 420*s3*(11*x - 20*z) + 840*s2*(5*pow(x,2) - 12*x*z + 6*z2) +
				   210*s*(9*pow(x,3) - 24*pow(x,2)*z + 12*x*z2 + 20*z3) + 168*(2*pow(x,4) - 5*pow(x,3)*z + 15*x*z3 - 20*z4)) -
		     14*l2*z7*(105840*s5 + 37800*s4*(6*x - 7*z) + 2100*s3*(120*pow(x,2) - 213*x*z + 116*z2) +
					   2520*s2*(60*pow(x,3) - 145*pow(x,2)*z + 120*x*z2 - 42*z3) +
					   840*pow(x - z,2)*(6*pow(x,3) - 12*pow(x,2)*z + 5*x*z2 - 2*z3) +
					   210*s*(216*pow(x,4) - 675*pow(x,3)*z + 760*pow(x,2)*z2 - 396*x*z3 + 108*z4)) +
		     15*l4*z5*(49392*s5 + 5040*s4*(21*x - 41*z) + 1680*s3*(70*pow(x,2) - 215*x*z + 186*z2) +
					   3360*s2*(21*pow(x,3) - 92*pow(x,2)*z + 120*x*z2 - 66*z3) +
					   84*s*(252*pow(x,4) - 1530*pow(x,3)*z + 2680*pow(x,2)*z2 - 2145*x*z3 + 910*z4) +
					   168*(14*pow(x,5) - 122*pow(x,4)*z + 270*pow(x,3)*z2 - 275*pow(x,2)*z3 + 165*x*z4 - 59*z5)) +
		     60*l6*z3*(8232*s5 + 2520*s4*(7*x + 2*z) + 140*s3*(140*pow(x,2) + 87*x*z - 288*z2) +
					   840*s2*(14*pow(x,3) + 17*pow(x,2)*z - 66*x*z2 + 56*z3) +
					   42*s*(84*pow(x,4) + 195*pow(x,3)*z - 800*pow(x,2)*z2 + 960*x*z3 - 600*z4) +
					   56*(7*pow(x,5) + 33*pow(x,4)*z - 135*pow(x,3)*z2 + 200*pow(x,2)*z3 - 180*x*z4 + 93*z5)) -
		     5*l8*z*(74088*s5 + 7560*s4*(21*x - 13*z) + 840*s3*(210*pow(x,2) - 183*x*z - 2*z2) +
				   5040*s2*(21*pow(x,3) - 22*pow(x,2)*z - 6*x*z2 + 27*z3) +
				   84*s*(378*pow(x,4) - 405*pow(x,3)*z - 460*pow(x,2)*z2 + 1665*x*z3 - 1980*z4) +
				   168*(21*pow(x,5) - 15*pow(x,4)*z - 85*pow(x,3)*z2 + 285*pow(x,2)*z3 - 450*x*z4 + 334*z5)) +
		     28*l5*r*z6*(15*r*(2100*s4 + 840*s3*(5*x - 4*z) + 420*pow(x,2)*pow(x - z,2) +
						   420*s2*(10*pow(x,2) - 12*x*z + 3*z2) + 420*s*x*(5*pow(x,2) - 8*x*z + 3*z2)) +
					     c*z2*(3150*s4 + 2100*s3*(3*x - 2*z) + 1260*s2*(5*pow(x,2) - 5*x*z + z2) +
							 210*pow(x,2)*(3*pow(x,2) - 5*x*z + 2*z2) + 210*s*x*(15*pow(x,2) - 20*x*z + 6*z2))) -
		     12*l7*r*z4*(6*r*(11550*s4 + 420*s3*(55*x - 56*z) + 210*pow(x - z,2)*(11*pow(x,2) - 6*x*z + 7*z2) +
						  420*s2*(55*pow(x,2) - 84*x*z + 45*z2) + 210*s*(55*pow(x,3) - 112*pow(x,2)*z + 90*x*z2 - 40*z3)) +
					     c*z2*(15960*s4 + 1680*s3*(19*x - 14*z) + 840*s2*(38*pow(x,2) - 42*x*z + 15*z2) +
							 840*s*(19*pow(x,3) - 28*pow(x,2)*z + 15*x*z2 - 5*z3) +
							 84*(38*pow(x,4) - 70*pow(x,3)*z + 50*pow(x,2)*z2 - 25*x*z3 + 7*z4))) +
		     12*l4*l5*r*z2*(c*z2*(10500*s4 + 840*s3*(25*x - 22*z) + 840*s2*(25*pow(x,2) - 33*x*z + 21*z2) +
							 420*s*(25*pow(x,3) - 44*pow(x,2)*z + 42*x*z2 - 24*z3) +
							 210*(10*pow(x,4) - 22*pow(x,3)*z + 28*pow(x,2)*z2 - 24*x*z3 + 9*z4)) +
					     15*r*(1260*s4 + 840*s3*(3*x - 4*z) + 2520*s2*(pow(x,2) - 2*x*z + 2*z2) +
						   84*s*(15*pow(x,3) - 40*pow(x,2)*z + 60*x*z2 - 48*z3) +
						   84*(3*pow(x,4) - 10*pow(x,3)*z + 20*pow(x,2)*z2 - 24*x*z3 + 13*z4))) -
		     12*l6*l7*r*(3*r*(540*s2 + 540*s*x + 180*(pow(x,2) - 3*z2)) -
				     c*(630*s4 + 420*s3*(3*x - 2*z) + 180*s2*(7*pow(x,2) - 7*x*z + 11*z2) +
					30*s*(21*pow(x,3) - 28*pow(x,2)*z + 66*x*z2 - 80*z3) +
					6*(21*pow(x,4) - 35*pow(x,3)*z + 110*pow(x,2)*z2 - 200*x*z3 + 170*z4))) +
		     8*l5*l6*r*(3*r*(1890*s4 + 3780*s3*x + 540*s2*(7*pow(x,2) - 15*z2) +
					 270*s*(7*pow(x,3) - 30*x*z2 + 40*z3) + 54*(7*pow(x,4) - 50*pow(x,2)*z2 + 100*x*z3 - 85*z4)) -
				    c*z2*(3780*s4 + 840*s3*(9*x - 10*z) + 360*s2*(21*pow(x,2) - 35*x*z + 40*z2) +
						60*s*(63*pow(x,3) - 140*pow(x,2)*z + 240*x*z2 - 198*z3) +
						12*(63*pow(x,4) - 175*pow(x,3)*z + 400*pow(x,2)*z2 - 495*x*z3 + 270*z4)))))/
	    (12.*l5*pow(l2 - z2,3)*(9*l6 - 63*l4*z2 - 21*l2*z4 - 245*z6));
    }
    else {
        BOOST_ASSERT_MSG(false, "IntegrateRhoDDerivative!");
    }
    result = result*_Rcn;
    return result;
}
double shape::IntegrateRhoDerivativeLowHigh(double z_high, double z_low, char label_i) {
    double result;
    if (label_i == 'l' || label_i == 's') {
        result = (IntegrateRhoDerivative(z_high, label_i) - Rho(z_low) * (-1) -
                                           IntegrateRhoDerivative(z_low, label_i));
    } else {
        result = (IntegrateRhoDerivative(z_high, label_i) - IntegrateRhoDerivative(z_low, label_i));
    }
    return result;
}
double shape::IntegrateRhoDDerivativeLowHigh(double z_high, double z_low, char label_i, char label_j){
    double result;
    if ((label_i == 'l' || label_i == 's')&&(label_j=='l'||label_j=='s')) {
        result = (IntegrateRhoDDerivative(z_high, label_i, label_j) -
				(RhoDerivative(z_low, 'x')+(-1)*RhoDerivative(z_low,label_i)+(-1)*RhoDerivative(z_low,label_j)+IntegrateRhoDDerivative(z_low,label_i,label_j)));
    } else if (label_j == 'l' || label_j == 's') {
        result = (IntegrateRhoDDerivative(z_high, label_i, label_j) - RhoDerivative(z_low, label_i) * (-1) -
                  IntegrateRhoDDerivative(z_low, label_i, label_j));
    } else if (label_i == 'l' || label_i == 's') {
		result = (IntegrateRhoDDerivative(z_high, label_i, label_j) - RhoDerivative(z_low, label_j) * (-1) -
				  IntegrateRhoDDerivative(z_low, label_i, label_j));
	} else{
        result = (IntegrateRhoDDerivative(z_high, label_i, label_j) - IntegrateRhoDDerivative(z_low, label_i, label_j));
    }
    return result;
}
double shape::CenterOfMassDerivative(const char side, const char label) {
    double l = _para_l;
    double r = _para_r;
    double z = _para_z;
    double c = _para_c;
    double s = _para_s;

    double l2=l*2, l3=l*l2, l4=l2*l2, l5=l2*l3, l6=l3*l3, l7=l3*l4, l8=l4*l4;
    double r2=r*2;
    double z2=z*2, z3=z*z2, z4=z2*z2, z5=z2*z3, z6=z3*z3, z7=z3*z4, z8=z4*z4;
    double s2=s*2;

    double result = 0;
    if (side=='L'){
	if(label=='l'){
	    result = ((l2 - z2)*(-1050*l4*l5 + 39*c*l6*l6*r - 210*l*z6*(8*s2 - 4*s*z + z2) +
					     840*l3*z4*(6*s2 - 3*s*z + 2*z2) -
					     1260*l5*z2*(4*s2 - 2*s*z + 3*z2) + 840*l7*(2*s2 - s*z + 4*z2) +
					     14*l6*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					     35*l4*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 66*l5*l5*r*(9*r + 2*c*z*(-2*s + 3*z)) -
					     18*l8*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		     7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		     4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		     12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		     12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3))) -
		((l2 - z2)*(-1050*l4*l5*s + 156*c*l6*l6*r*z + 280*l7*(15*s - 16*z)*z2 +
					24*l3*(175*s - 184*z)*z6 - 14*l*(75*s - 16*z)*z8 - 252*l5*z4*(25*s + 24*z) +
					44*l5*l5*(18 + 27*r2*z - 2*c*r*z3) - 140*l4*z6*(70 + 45*r2*z + c*r*z3) -
					108*l8*z2*(42 + 3*r2*z + 4*c*r*z3) +
					84*l6*z4*(-14 + 99*r2*z + 6*c*r*z3))*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 +
			6*l4*(175*s - 184*z)*z6 - 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 -
			42*l6*z4*(25*s + 24*z) + 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) -
			28*l5*z6*(70 + 45*r2*z + c*r*z3) -
			12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
			12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3),2)) +
		(l*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		    105*l2*z6*(8*s2 - 4*s*z + z2) +
		    210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		    210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		    2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		    7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		    2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		 12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		 12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3));
	}
	else if(label=='r'){
	    result = ((l2 - z2)*(3*c*l6*l7 + 54*l5*l6*r - 54*l4*l5*r*(4*s - z)*z - 18*l7*r*(8*s - 49*z)*z3 -
					     630*l5*r*(4*s - 3*z)*z5 + 2*l7*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					     7*l5*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*(9*r + 2*c*z*(-2*s + 3*z)) -
					     2*l4*l5*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		     7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		     4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		     12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		     12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3))) -
		((l2 - z2)*(12*c*l6*l7*z + 4*l5*l6*(54*r*z - 2*c*z3) -
					28*l5*z6*(90*r*z + c*z3) - 12*l4*l5*z2*(6*r*z + 4*c*z3) +
					12*l7*z4*(198*r*z + 6*c*z3))*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 +
			6*l4*(175*s - 184*z)*z6 - 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 -
			42*l6*z4*(25*s + 24*z) + 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) -
			28*l5*z6*(70 + 45*r2*z + c*r*z3) -
			12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
			12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3),2));
	}
	else if(label=='z'){
	    result = ((l2 - z2)*(840*s*(2*s - z)*z7 - 105*s*z8 - 105*l2*z6*(-4*s + 2*z) +
					     210*l4*z4*(-3*s + 4*z) - 210*l6*z2*(-2*s + 6*z) + 105*l8*(-s + 8*z) -
					     630*l2*z5*(8*s2 - 4*s*z + z2) +
					     840*l4*z3*(6*s2 - 3*s*z + 2*z2) - 420*l6*z*(4*s2 - 2*s*z + 3*z2) +
					     2*l7*r*z3*(441*r + 4*c*(22*s - 7*z)*z - 14*c*z2) -
					     7*l5*r*z5*(-270*r + 2*c*(8*s - 5*z)*z - 5*c*z2) +
					     6*l7*r*z2*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					     35*l5*r*z4*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(6*c*z + 2*c*(-2*s + 3*z)) -
					     2*l4*l5*r*z*(-27*r + 23*c*z2 + 2*c*z*(4*s + 23*z)) - 2*l4*l5*r*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		     7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		     4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		     12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		     12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3))) -
		((l2 - z2)*(12*c*l6*l7*r + 70*l8*(15*s - 16*z)*z - 560*l8*z2 -
					1008*l6*z4 + 36*l4*(175*s - 184*z)*z5 - 1104*l4*z6 -
					56*l2*(75*s - 16*z)*z7 + 112*l2*z8 + 1050*s*z4*z5 -
					168*l6*z3*(25*s + 24*z) + 4*l5*l6*(27*r2 - 6*c*r*z2) -
					28*l5*z6*(45*r2 + 3*c*r*z2) - 12*l4*l5*z2*(3*r2 + 12*c*r*z2) +
					12*l7*z4*(99*r2 + 18*c*r*z2) -
					168*l5*z5*(70 + 45*r2*z + c*r*z3) - 24*l4*l5*z*(42 + 3*r2*z + 4*c*r*z3) +
					48*l7*z3*(-14 + 99*r2*z + 6*c*r*z3))*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 +
			6*l4*(175*s - 184*z)*z6 - 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 -
			42*l6*z4*(25*s + 24*z) + 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) -
			28*l5*z6*(70 + 45*r2*z + c*r*z3) -
			12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
			12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3),2)) -
		(z*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		    105*l2*z6*(8*s2 - 4*s*z + z2) +
		    210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		    210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		    2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		    7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		    2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		 12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		 12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3));
	}
	else if(label=='c'){
	    result = ((l2 - z2)*(3*l6*l7*r + 4*l7*r*(22*s - 7*z)*z5 - 7*l5*r*(8*s - 5*z)*z7 +
					     12*l5*l6*r*z*(-2*s + 3*z) - 2*l4*l5*r*z3*(4*s + 23*z)))/
		(2.*(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		     7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		     4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		     12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		     12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3))) -
		((l2 - z2)*(12*l6*l7*r*z - 8*l5*l6*r*z3 - 48*l4*l5*r*z5 +
					72*l7*r*z7 - 28*l5*r*z4*z5)*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 +
			6*l4*(175*s - 184*z)*z6 - 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 -
			42*l6*z4*(25*s + 24*z) + 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) -
			28*l5*z6*(70 + 45*r2*z + c*r*z3) -
			12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
			12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3),2));
	}
	else if(label=='s'){
	    result = ((l2 - z2)*(105*l8*(4*s - z) - 24*c*l5*l6*r*z - 210*l6*(8*s - 2*z)*z2 +
					     210*l4*(12*s - 3*z)*z4 - 105*l2*(16*s - 4*z)*z6 + 210*s*z8 + 105*(2*s - z)*z8 -
					     2*l4*l5*r*z*(108*r + 4*c*z2) - 7*l5*r*z5*(360*r + 8*c*z2) +
					     2*l7*r*z3*(-72*r + 44*c*z2)))/
		(2.*(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 + 6*l4*(175*s - 184*z)*z6 -
		     7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 - 42*l6*z4*(25*s + 24*z) +
		     4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) - 28*l5*z6*(70 + 45*r2*z + c*r*z3) -
		     12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
		     12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3))) -
		((l2 - z2)*(-105*l5*l5 + 525*l8*z2 - 1050*l6*z4 + 1050*l4*z6 -
					525*l2*z8 + 105*z5*z5)*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
										105*l2*z6*(8*s2 - 4*s*z + z2) +
										210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
										210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
										2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
										7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
										2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(-105*l5*l5*s + 12*c*l6*l7*r*z + 35*l8*(15*s - 16*z)*z2 +
			6*l4*(175*s - 184*z)*z6 - 7*l2*(75*s - 16*z)*z8 + 105*s*z5*z5 -
			42*l6*z4*(25*s + 24*z) + 4*l5*l6*(18 + 27*r2*z - 2*c*r*z3) -
			28*l5*z6*(70 + 45*r2*z + c*r*z3) -
			12*l4*l5*z2*(42 + 3*r2*z + 4*c*r*z3) +
			12*l7*z4*(-14 + 99*r2*z + 6*c*r*z3),2));
	}
	else{
	    BOOST_ASSERT_MSG(false,"CenterOfMassDerivative!");
	}
    }

    else if(side=='R'){
	if(label=='l'){
	    result = -((l2 - z2)*(-1050*l4*l5 + 39*c*l6*l6*r - 210*l*z6*(8*s2 - 4*s*z + z2) +
					      840*l3*z4*(6*s2 - 3*s*z + 2*z2) -
					      1260*l5*z2*(4*s2 - 2*s*z + 3*z2) + 840*l7*(2*s2 - s*z + 4*z2) +
					      14*l6*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					      35*l4*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 66*l5*l5*r*(9*r + 2*c*z*(-2*s + 3*z)) -
					      18*l8*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		     105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		     28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		     12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		     12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3))) +
		((l2 - z2)*(1050*l4*l5*s - 156*c*l6*l6*r*z - 24*l3*(175*s - 184*z)*z6 +
					14*l*(75*s - 16*z)*z8 + 280*l7*z2*(-15*s + 16*z) + 252*l5*z4*(25*s + 24*z) +
					140*l4*z6*(-70 + 45*r2*z + c*r*z3) +
					108*l8*z2*(-42 + 3*r2*z + 4*c*r*z3) -
					84*l6*z4*(14 + 99*r2*z + 6*c*r*z3) + 11*l5*l5*(72 - 108*r2*z + 8*c*r*z3))*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 +
			7*l2*(75*s - 16*z)*z8 - 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) +
			42*l6*z4*(25*s + 24*z) + 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
			12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
			12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3),2))
		- (l*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		      105*l2*z6*(8*s2 - 4*s*z + z2) +
		      210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		      210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		      2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		      7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		      2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		 12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		 12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3));
	}
	else if(label=='r'){
	    result = -((l2 - z2)*(3*c*l6*l7 + 54*l5*l6*r - 54*l4*l5*r*(4*s - z)*z -
					      18*l7*r*(8*s - 49*z)*z3 - 630*l5*r*(4*s - 3*z)*z5 +
					      2*l7*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					      7*l5*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*(9*r + 2*c*z*(-2*s + 3*z)) -
					      2*l4*l5*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		     105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		     28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		     12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		     12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3))) +
		((l2 - z2)*(-12*c*l6*l7*z + 28*l5*z6*(90*r*z + c*z3) +
					12*l4*l5*z2*(6*r*z + 4*c*z3) - 12*l7*z4*(198*r*z + 6*c*z3) +
					l5*l6*(-216*r*z + 8*c*z3))*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
									      105*l2*z6*(8*s2 - 4*s*z + z2) +
									      210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
									      210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
									      2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
									      7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
									      2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 +
			7*l2*(75*s - 16*z)*z8 - 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) +
			42*l6*z4*(25*s + 24*z) + 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
			12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
			12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3),2));
	}
	else if(label=='z'){
	    result = -((l2 - z2)*(840*s*(2*s - z)*z7 - 105*s*z8 - 105*l2*z6*(-4*s + 2*z) +
					      210*l4*z4*(-3*s + 4*z) - 210*l6*z2*(-2*s + 6*z) + 105*l8*(-s + 8*z) -
					      630*l2*z5*(8*s2 - 4*s*z + z2) +
					      840*l4*z3*(6*s2 - 3*s*z + 2*z2) - 420*l6*z*(4*s2 - 2*s*z + 3*z2) +
					      2*l7*r*z3*(441*r + 4*c*(22*s - 7*z)*z - 14*c*z2) -
					      7*l5*r*z5*(-270*r + 2*c*(8*s - 5*z)*z - 5*c*z2) +
					      6*l7*r*z2*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
					      35*l5*r*z4*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(6*c*z + 2*c*(-2*s + 3*z)) -
					      2*l4*l5*r*z*(-27*r + 23*c*z2 + 2*c*z*(4*s + 23*z)) - 2*l4*l5*r*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		     105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		     28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		     12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		     12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3))) +
		((l2 - z2)*(-12*c*l6*l7*r + 560*l8*z2 + 1008*l6*z4 -
					36*l4*(175*s - 184*z)*z5 + 1104*l4*z6 + 56*l2*(75*s - 16*z)*z7 -
					112*l2*z8 - 1050*s*z4*z5 + 70*l8*z*(-15*s + 16*z) + 168*l6*z3*(25*s + 24*z) +
					28*l5*z6*(45*r2 + 3*c*r*z2) + 12*l4*l5*z2*(3*r2 + 12*c*r*z2) -
					12*l7*z4*(99*r2 + 18*c*r*z2) + l5*l6*(-108*r2 + 24*c*r*z2) +
					168*l5*z5*(-70 + 45*r2*z + c*r*z3) + 24*l4*l5*z*(-42 + 3*r2*z + 4*c*r*z3) -
					48*l7*z3*(14 + 99*r2*z + 6*c*r*z3))*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 +
			7*l2*(75*s - 16*z)*z8 - 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) +
			42*l6*z4*(25*s + 24*z) + 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
			12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
			12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3),2))
		+ (z*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		      105*l2*z6*(8*s2 - 4*s*z + z2) +
		      210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		      210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		      2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		      7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		      2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		 12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		 12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3));
	}
	else if(label=='c'){
	    result = -((l2 - z2)*(3*l6*l7*r + 4*l7*r*(22*s - 7*z)*z5 - 7*l5*r*(8*s - 5*z)*z7 +
					      12*l5*l6*r*z*(-2*s + 3*z) - 2*l4*l5*r*z3*(4*s + 23*z)))/
		(2.*(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		     105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		     28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		     12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		     12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3))) +
		((l2 - z2)*(-12*l6*l7*r*z + 8*l5*l6*r*z3 + 48*l4*l5*r*z5 -
					72*l7*r*z7 + 28*l5*r*z4*z5)*
		 (-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
		  105*l2*z6*(8*s2 - 4*s*z + z2) +
		  210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
		  210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
		  2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
		  7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
		  2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 +
			7*l2*(75*s - 16*z)*z8 - 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) +
			42*l6*z4*(25*s + 24*z) + 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
			12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
			12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3),2));
	}
	else if(label=='s'){
	    result = -((l2 - z2)*(105*l8*(4*s - z) - 24*c*l5*l6*r*z - 210*l6*(8*s - 2*z)*z2 +
					      210*l4*(12*s - 3*z)*z4 - 105*l2*(16*s - 4*z)*z6 + 210*s*z8 + 105*(2*s - z)*z8 -
					      2*l4*l5*r*z*(108*r + 4*c*z2) - 7*l5*r*z5*(360*r + 8*c*z2) +
					      2*l7*r*z3*(-72*r + 44*c*z2)))/
		(2.*(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 + 7*l2*(75*s - 16*z)*z8 -
		     105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) + 42*l6*z4*(25*s + 24*z) +
		     28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
		     12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
		     12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3))) +
		((l2 - z2)*(105*l5*l5 - 525*l8*z2 + 1050*l6*z4 - 1050*l4*z6 +
					525*l2*z8 - 105*z5*z5)*(-105*l5*l5 + 3*c*l6*l7*r + 105*s*(2*s - z)*z8 -
										105*l2*z6*(8*s2 - 4*s*z + z2) +
										210*l4*z4*(6*s2 - 3*s*z + 2*z2) -
										210*l6*z2*(4*s2 - 2*s*z + 3*z2) + 105*l8*(2*s2 - s*z + 4*z2) +
										2*l7*r*z3*(-9*r*(8*s - 49*z) + 2*c*(22*s - 7*z)*z2) -
										7*l5*r*z5*(90*r*(4*s - 3*z) + c*(8*s - 5*z)*z2) + 6*l5*l6*r*(9*r + 2*c*z*(-2*s + 3*z)) -
										2*l4*l5*r*z*(27*r*(4*s - z) + c*z2*(4*s + 23*z))))/
		(2.*pow(105*l5*l5*s - 12*c*l6*l7*r*z - 6*l4*(175*s - 184*z)*z6 +
			7*l2*(75*s - 16*z)*z8 - 105*s*z5*z5 + 35*l8*z2*(-15*s + 16*z) +
			42*l6*z4*(25*s + 24*z) + 28*l5*z6*(-70 + 45*r2*z + c*r*z3) +
			12*l4*l5*z2*(-42 + 3*r2*z + 4*c*r*z3) -
			12*l7*z4*(14 + 99*r2*z + 6*c*r*z3) + l5*l6*(72 - 108*r2*z + 8*c*r*z3),2));
	}
	else{
	    BOOST_ASSERT_MSG(false, "CenterOfMassDerivative!");
	}
    }
    else{
        BOOST_ASSERT_MSG(false, "CenterOfMassDerivative!");
    }
    return result;
}
double RhoShape(double zeta, void *params) {
    shape* tmp_shape = (shape*) params;
	double result;
    result = -tmp_shape->Rho(zeta);
//    for(int i=0;i<201;i++){
//        double x=tmp_shape->get_l()*(i-100)/100-tmp_shape->get_s();
//        cout<<x<<' '<<tmp_shape->Rho(x)<<endl;
//    }
//	cout<<zeta<<' '<<result<<endl;
    return result;
}
//  A_i=-\frac{1}{\rho(z;q)^2} \frac{\partial}{\partial q_i} \int_{z_{min}}^z \rho(z';q)^2 dz'
double A(shape shape, double z_high, double z_low, char label_i){
    double result;
    result = -1/shape.Rho(z_high)*shape.IntegrateRhoDerivativeLowHigh(z_high, z_low, label_i);
    return result;
}
double ADerivativeZ(shape shape, double z_high, double z_low, char label_i){
    double result;
    result = 1/pow(shape.Rho(z_high),2)*shape.RhoDerivative(z_high,'x')*shape.IntegrateRhoDerivativeLowHigh(z_high, z_low, label_i)-
	1/shape.Rho(z_high)*shape.RhoDerivative(z_high, label_i);
    return result;
}
double ADerivativeQ(shape shape, double z_high, double z_low, char label_i, char label_l){
    char lrzcs[6]="lrzcs";
    int i,l;
    i = distance(lrzcs, find(lrzcs, lrzcs + 5, label_i));
    l = distance(lrzcs, find(lrzcs, lrzcs + 5, label_l));
    double unswap_result, swap_result, result;
    unswap_result = 1/pow(shape.Rho(z_high),2)*shape.RhoDerivative(z_high, label_l)*shape.IntegrateRhoDerivativeLowHigh(z_high, z_low, label_i);
    if (i>l){
        swap(label_i,label_l);
    }
	swap_result = -1/shape.Rho(z_high)*shape.IntegrateRhoDDerivativeLowHigh(z_high, z_low, label_i, label_l);
    result = unswap_result+swap_result;
    return result;
}
double ADDerivativeZQ(shape shape, double z_high, double z_low, char label_i, char label_l){
    char lrzcs[6]="lrzcs";
    int i,l;
    i = distance(lrzcs, find(lrzcs, lrzcs + 5, label_i));
    l = distance(lrzcs, find(lrzcs, lrzcs + 5, label_l));
    double unswap_result, swap_result, result;
    unswap_result = -2./pow(shape.Rho(z_high),3)*shape.RhoDerivative(z_high, label_l)*shape.RhoDerivative(z_high,'x')*shape.IntegrateRhoDerivativeLowHigh(z_high, z_low, label_i)+
            1/pow(shape.Rho(z_high),2)*shape.RhoDDerivative(z_high, 'x', label_l)*shape.IntegrateRhoDerivativeLowHigh(z_high,z_low,label_i)+
            1/pow(shape.Rho(z_high),2)*shape.RhoDerivative(z_high, label_l)*shape.RhoDerivative(z_high, label_i);
    if (i>l){
        swap(label_i,label_l);
    }
	swap_result = 1/pow(shape.Rho(z_high),2)*shape.RhoDerivative(z_high, 'x')*shape.IntegrateRhoDDerivativeLowHigh(z_high, z_low, label_i, label_l)-
            1/shape.Rho(z_high)*shape.RhoDDerivative(z_high, label_i, label_l);
    result = unswap_result+swap_result;
    return result;
}
