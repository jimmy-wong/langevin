#include <iostream>
#include "../include/shape.h"

using namespace std;

void input(const string&, int*, double*, double*, int*);
void store(const string&, int*, double*);
double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
void runge_kutta(gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double starting[5], double step_length[5], double storation[],
                 shape shape);
// 命令行参数是激发能
int main(int argc, char* argv[])
{
    shape shape;
    int steps[5],gs[5],sp[5];
    double starting[5],step_length[5];
    input("U236_Input.dat",steps,starting,step_length,gs);
    gs[0]=1;gs[1]=1;gs[2]=2;gs[3]=0;gs[4]=2;
    starting[1] = 1 / sqrt(starting[0]+gs[0]*step_length[0]) + 0.05;
    starting[3] = -pow(starting[0]+gs[0]*step_length[0], -2.5) - 0.5;
    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
    shape.set_steps(steps);
    shape.set_level_density();
    // 初始化storation为20.
    fill(storation,storation+steps[0]*steps[1]*steps[2]*steps[3]*steps[4],100.);
    store("U236.txt",steps,storation);
    double ground_energy = shape.grid_energy(storation, gs);
    shape.set_ground_energy(ground_energy);
    shape.set_excited_energy(stof(argv[1]));
    // storation 中间存储的就是势能曲面的值
    // inertia_tensor是惯性张量，dissipative_tenson是耗散张量
    gsl_vector *generalized_coordinates, *generalized_momenta;
    generalized_coordinates = gsl_vector_alloc(5);
    generalized_momenta = gsl_vector_alloc(5);

    //这里gamma_tensor的结果就是dissipative_tensor的平方根
    for(int i=0; i<10000; i++) {
        sp[0]=15;sp[1]=13;sp[2]=8;sp[3]=18;sp[4]=11;
        starting[1] = 1 / sqrt(starting[0]+sp[0]*step_length[0]) + 0.05;
        starting[3] = -pow(starting[0]+sp[0]*step_length[0], -2.5) - 0.5;

        for (size_t i=0; i<5; i++){
            gsl_vector_set(generalized_momenta, i, 0);
            gsl_vector_set(generalized_coordinates, i, sp[i]*step_length[i]+starting[i]);
        }

        runge_kutta(generalized_coordinates, generalized_momenta,
                    starting, step_length, storation,
                    shape);
        cout<<shape.AH(generalized_coordinates)<<endl;
    }
    delete(storation);
    gsl_vector_free(generalized_coordinates);
    gsl_vector_free(generalized_momenta);

    return 0;
}
