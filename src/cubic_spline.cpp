#include "../include/shape.h"

// useful for the calculation of Temperature
double hypercubic_interp(shape shape, double* starting, double* step_length, double* storation) {

    int grid[5];
    shape.grid(starting, step_length, grid);
    starting[1] = 1/sqrt(starting[0]+grid[0]*step_length[0])+0.05;
    starting[3] = -pow(starting[0]+grid[0]*step_length[0],-2.5)-0.5;
    shape.grid(starting, step_length, grid);
    double d_grid[5][2];
    for(int i=0;i<5;i++){
        d_grid[i][0] = starting[i]+grid[i]*step_length[i];
        d_grid[i][1] = d_grid[i][0]+step_length[i];
    }
    double tmp_res = 0, result = 0;
    double l1 = d_grid[0][1], l0= d_grid[0][0];
    double c1 = -pow(l1, -2.5) - 0.5, c0 = -pow(l0, -2.5) - 0.5;
    double r1 = 1 / sqrt(l1) + 0.05, r0 = 1 / sqrt(l0) + 0.05;
    double slope_c = (c0-c1)/(l0-l1), slope_r = (r0-r1)/(l0-l1);
    while (shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][0]<0){
        d_grid[3][0] = d_grid[3][0]-step_length[3];
        d_grid[3][1] = d_grid[3][1]-step_length[3];
    }
    while ((shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][0])/step_length[1]<0){
        d_grid[1][0] = d_grid[1][0]-step_length[1];
        d_grid[1][1] = d_grid[1][1]-step_length[1];
    }
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
                        tmp_res = (1-2*m)*(shape.get_s()-d_grid[4][m])/(step_length[4])*
                                (1-2*l)*(shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][l])/(step_length[3])*
                                (1-2*k)*(shape.get_z()-d_grid[2][k])/(step_length[2])*
                                (1-2*j)*(shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][j])/(step_length[1])*
                                (1-2*i)*(shape.get_l()-d_grid[0][i])/(step_length[0]);
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
    int grid[5];
    shape.grid(starting, step_length, grid);
    starting[1] = 1/sqrt(starting[0]+grid[0]*step_length[0])+0.05;
    starting[3] = -pow(starting[0]+grid[0]*step_length[0], -2.5)-0.5;
    shape.grid(starting, step_length, grid);

    double d_grid[5][2];
    for(int i=0;i<5;i++){
        d_grid[i][0] = starting[i]+grid[i]*step_length[i];
        d_grid[i][1] = d_grid[i][0]+step_length[i];
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
    if(label=='s'){
        cycle[4] = true;
    }
    double tmp_res_lrc, tmp_res_l, tmp_res_r, tmp_res_z, tmp_res_c, tmp_res_s, result = 0;
    double l1 = d_grid[0][1], l0= d_grid[0][0];
    double c1 = -pow(l1, -2.5) - 0.5, c0 = -pow(l0, -2.5) - 0.5;
    double r1 = 1 / sqrt(l1) + 0.05, r0 = 1 / sqrt(l0) + 0.05;
    double slope_c = (c0-c1)/(l0-l1), slope_r = (r0-r1)/(l0-l1);
    while (shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][0]<0){
        d_grid[3][0] = d_grid[3][0]-step_length[3];
        d_grid[3][1] = d_grid[3][1]-step_length[3];
    }
    while ((shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][0])/step_length[1]<0){
        d_grid[1][0] = d_grid[1][0]-step_length[1];
        d_grid[1][1] = d_grid[1][1]-step_length[1];
    }
    int step[5];
    for(int i=0;i<2;i++){
        step[0] = 1-i+grid[0];
        if(cycle[0]){
            for(int j=0;j<2;j++) {
                step[1] = 1-j+grid[1];
                for(int l=0;l<2;l++) {
                    step[3] = 1-l+grid[3];
                    tmp_res_lrc = (1-2*i)/(step_length[0])*
                            (1-2*j)*(shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][j])/(step_length[1])*
                            (1-2*l)*(shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][l])/(step_length[3])+
                            (1-2*i)*(shape.get_l()-d_grid[0][i])/(step_length[0])*
                            (1-2*j)*(-slope_r)/(step_length[1])*
                            (1-2*l)*(shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][l])/(step_length[3])+
                            (1-2*i)*(shape.get_l()-d_grid[0][i])/(step_length[0])*
                            (1-2*j)*(shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][j])/(step_length[1])*
                            (1-2*l)*(-slope_c)/(step_length[3]);
                    for(int k=0;k<2;k++){
                        step[2] = 1-k+grid[2];
                        tmp_res_z = (1-2*k)*(shape.get_z()-d_grid[2][k])/(step_length[2]);
                        for(int m=0;m<2;m++) {
                            step[4] = 1-m+grid[4];
                            tmp_res_s = (1-2*m)*(shape.get_s()-d_grid[4][m])/(step_length[4]);
                            result += tmp_res_lrc*tmp_res_z*tmp_res_s*
                                      shape.grid_energy(storation,step);
                        }
                    }
                }
            }
        }
        else{
            tmp_res_l = (1-2*i)*(shape.get_l()-d_grid[0][i])/(step_length[0]);
            for(int j=0;j<2;j++) {
                step[1] = 1 - j + grid[1];
                if (cycle[1]) {
                    tmp_res_r = (1-2*j)/(step_length[1]);
                } else {
                    tmp_res_r = (1-2*j)*(shape.get_r()+slope_r*l0-slope_r*shape.get_l()-d_grid[1][j])/(step_length[1]);
                }
                for(int l=0;l<2;l++) {
                    step[3] = 1-l+grid[3];
                    if (cycle[3]) {
                        tmp_res_c = (1-2*l)/(step_length[3]);
                    } else {
                        tmp_res_c = (1-2*l)*(shape.get_c()+slope_c*l0-slope_c*shape.get_l()-d_grid[3][l])/(step_length[3]);
                    }
                    for(int k=0;k<2;k++){
                        step[2] = 1-k+grid[2];
                        if(cycle[2]){
                            tmp_res_z = (1-2*k)/(step_length[2]);
                        }
                        else{
                            tmp_res_z = (1-2*k)*(shape.get_z()-d_grid[2][k])/(step_length[2]);
                        }
                        for(int m=0;m<2;m++) {
                            step[4] = 1-m+grid[4];
                            if (cycle[4]) {
                                tmp_res_s = (1-2*m)/(step_length[4]);
                            }
                            else{
                                tmp_res_s = (1-2*m)*(shape.get_s()-d_grid[4][m])/(step_length[4]);
                            }
                            result += tmp_res_l*tmp_res_r*tmp_res_c*tmp_res_z*tmp_res_s*
                                      shape.grid_energy(storation,step);
                        }
                    }
                }
            }
        }
    }
    return result/shape.get_Rcn();
}
