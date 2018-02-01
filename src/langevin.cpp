#include <iostream>
#include "../include/shape.h"

using namespace std;

void input(const string&, int*);
void store(const string&, int*, double*);
int main()
{
    int steps[5];
    input("fort.112",steps);
    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
    // 初始化storation为20
    fill(storation,storation+steps[0]*steps[1]*steps[2]*steps[3]*steps[4],20.);
    store("U236.txt",steps,storation);
    // storation 中间存储的就是势能曲面的值

    return 0;
}
