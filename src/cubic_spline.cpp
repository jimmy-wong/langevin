#include "../include/shape.h"
#include <iostream>
#include <cmath>

// useful for the calculation of Temperature
double hypercubic_interp(shape shape, double* starting, double* step_length, double* storation) {
    starting[1] = 1 / sqrt(starting[0]) + 0.05;
    starting[3] = -pow(starting[0], -2.5) - 0.5;
    int grid[5];
    shape.grid(starting, step_length, grid);
    double d_grid[5][2];
    for(int i=0;i<5;i++){
        d_grid[i][0] = starting[i]+grid[i]*step_length[i];
        d_grid[i][1] = starting[i]+(grid[i]+1)*step_length[i];
    }
    double tmp_res, result = 0;
    int step[5];
    for(int i=0;i<2;i++){
        step[0] = 1-i+grid[0];
        for(int j=0;j<2;j++){
            step[1] = 1-j+grid[1];
            for(int k=0;k<2;k++){
                step[2] = 1-k+grid[2];
                for(int l=0;l<2;l++){
                    step[3] = 1-l+grid[3];
                    for(int m=0;m<2;m++) {
                        step[4] = 1-m+grid[4];
                        tmp_res = (1-2*m)*(shape.get_s()-d_grid[4][m])/(d_grid[4][1]-d_grid[4][0])*
                                (1-2*l)*(shape.get_c()-d_grid[3][l])/(d_grid[3][1]-d_grid[3][0])*
                                (1-2*k)*(shape.get_z()-d_grid[2][k])/(d_grid[2][1]-d_grid[2][0])*
                                (1-2*j)*(shape.get_r()-d_grid[1][j])/(d_grid[1][1]-d_grid[1][0])*
                                (1-2*i)*(shape.get_l()-d_grid[0][i])/(d_grid[0][1]-d_grid[0][0]);
                        result += tmp_res*shape.grid_energy(storation,step);
                    }
                }
            }
        }
    }
    return result;
}

//useful for the calculation of Langevin equation
double hypercubic_interp_df(shape shape, double* starting, double* step_length, double* storation, const char label){
    cout<<"hypercubic_interp"<<endl;
//    cout<<starting[0]<<' '<<starting[1]<<' '<<starting[2]<<' '<<starting[3]<<' '<<starting[4]<<' '<<endl;
//    cout<<step_length[0]<<' '<<step_length[1]<<' '<<step_length[2]<<' '<<step_length[3]<<' '<<step_length[4]<<' '<<endl;
    starting[1] = 1 / sqrt(starting[0]) + 0.05;
    starting[3] = -pow(starting[0], -2.5) - 0.5;
    int grid[5];
    shape.grid(starting, step_length,grid);
//    std::cout<<storation[0]<<' '<<storation[1]<<std::endl;
    cout<<"grid: "<<grid[0]<<' '<<grid[1]<<' '<<grid[2]<<' '<<grid[3]<<' '<<grid[4]<<' '<<endl;
    cout<<"end hypercubic_interp"<<endl;
    double d_grid[5][2];
    for(int i=0;i<5;i++){
        d_grid[i][0] = starting[i]+grid[i]*step_length[i];
        d_grid[i][1] = starting[i]+(grid[i]+1)*step_length[i];
    }
    bool cycle[5] = {false,false,false,false,false};
    if(label=='l'){
        cycle[0] = true;
    }
    if(label=='r'){
        cycle[1] = true;
    }
    if(label=='z'){
        cycle[2] = true;
    }
    if(label=='c'){
        cycle[3] = true;
    }
    if(label=='z'){
        cycle[4] = true;
    }
    double tmp_res_l, tmp_res_r, tmp_res_z, tmp_res_c, tmp_res_s, result = 0;
    int step[5];
    for(int i=0;i<2;i++){
        step[0] = 1-i+grid[0];
        if (cycle[0]){
            tmp_res_l = (1-2*i)/(d_grid[0][1]-d_grid[0][0]);
        }
        else{
            tmp_res_l = (1-2*i)*(shape.get_l()-d_grid[0][i])/(d_grid[0][1]-d_grid[0][0]);
        }
        for(int j=0;j<2;j++){
            step[1] = 1-j+grid[1];
            if (cycle[1]){
                tmp_res_r = (1-2*j)/(d_grid[1][1]-d_grid[1][0]);
            }
            else{
                tmp_res_r = (1-2*j)*(shape.get_r()-d_grid[1][j])/(d_grid[1][1]-d_grid[1][0]);
            }
            for(int k=0;k<2;k++){
                step[2] = 1-k+grid[2];
                if(cycle[2]){
                    tmp_res_z = (1-2*k)/(d_grid[2][1]-d_grid[2][0]);
                }
                else{
                    tmp_res_z = (1-2*k)*(shape.get_z()-d_grid[2][k])/(d_grid[2][1]-d_grid[2][0]);
                }
                for(int l=0;l<2;l++){
                    step[3] = 1-l+grid[3];
                    if (cycle[3]) {
                        tmp_res_c = (1-2*l)/(d_grid[3][1]-d_grid[3][0]);
                    }
                    else{
                        tmp_res_c = (1-2*l)*(shape.get_c()-d_grid[3][l])/(d_grid[3][1]-d_grid[3][0]);
                    }
                    for(int m=0;m<2;m++) {
                        step[4] = 1-m+grid[4];
                        if (cycle[4]) {
                            tmp_res_s = (1 - 2 * m) / (d_grid[4][1] - d_grid[4][0]);
                        }
                        else{
                            tmp_res_s = (1 - 2 * m) * (shape.get_s() - d_grid[4][m]) / (d_grid[4][1] - d_grid[4][0]);
                        }
                        result += tmp_res_l*tmp_res_r*tmp_res_z*tmp_res_c*tmp_res_s*
                                shape.grid_energy(storation,step);
                    }
                }
            }
        }
    }
    return result;
}
