#include <iostream>
#include "../include/shape.h"
#include <gsl/gsl_matrix.h>
#include <cmath>
#include <gsl/gsl_vector.h>

using namespace std;

void input(const string&, int*, double*, double*, int*);
void store(const string&, int*, double*);
double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
void runge_kutta(gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double starting[5], double step_length[5], int steps[5], double storation[],
                 shape shape);
int main()
{
    shape shape;
    int steps[5],gs[5];
    double starting[5],step_length[5];
    input("fort.112",steps,starting,step_length,gs);
    starting[1] = 1 / sqrt(starting[0]) + 0.05;
    starting[3] = -pow(starting[0], -2.5) - 0.5;
    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
    shape.set_steps(steps);
    // 初始化storation为20.
    fill(storation,storation+steps[0]*steps[1]*steps[2]*steps[3]*steps[4],20.);
    store("U236.txt",steps,storation);
    // storation 中间存储的就是势能曲面的值
    // inertia_tensor是惯性张量，dissipative_tenson是耗散张量
    gsl_vector *generalized_coordinates, *generalized_momenta;
    generalized_coordinates = gsl_vector_alloc(5);
    generalized_momenta = gsl_vector_alloc(5);
    gs[0]=1;gs[1]=1;gs[2]=0;gs[3]=0;gs[4]=0;
    for (size_t i=0; i<5; i++){
        gsl_vector_set(generalized_momenta, i, 0);
        gsl_vector_set(generalized_coordinates, i, gs[i]*step_length[i]+starting[i]);
    }

    //这里gamma_tensor的结果就是dissipative_tensor的平方根
    for(int i=0; i<10000; i++) {
        cout<<i<<' '<<storation[0]<<' '<<storation[1]<<endl;
        cout<<starting[0]<<' '<<starting[1]<<' '<<starting[2]<<' '<<endl;
        cout<<step_length[0]<<' '<<step_length[1]<<' '<<step_length[2]<<' '<<endl;
        runge_kutta(generalized_coordinates, generalized_momenta,
                    starting, step_length, steps, storation,
                    shape);
    }
    delete(storation);
    gsl_vector_free(generalized_coordinates);
    gsl_vector_free(generalized_momenta);

    return 0;
}
