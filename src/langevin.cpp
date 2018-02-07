#include <iostream>
#include "../include/shape.h"

using namespace std;

void input(const string&, int*, double*, double*);
void store(const string&, const int*, double*);
double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, double density, const char label_i, const char label_j);
int main()
{
    shape shape;
    int steps[5];
    double starting[5],step_length[5];
    input("fort.112",steps,starting,step_length);
    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
    shape.set_steps(steps);
    // 初始化storation为20
    fill(storation,storation+steps[0]*steps[1]*steps[2]*steps[3]*steps[4],20.);
    store("U236.txt",steps,storation);
    // storation 中间存储的就是势能曲面的值
    double inertia_tensor[5][5];
    double dissipative_tensor[5][5];
    char* label="lrzcs";

    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            inertia_tensor[i][j]=inertia(shape,shape.get_density(),label[i],label[j]);
            dissipative_tensor[i][j]=dissipative(shape,label[i],label[j]);
        }
    }

    return 0;
}
