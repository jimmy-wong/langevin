#include <iostream>
#include "../include/shape.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <cmath>

using namespace std;

void input(const string&, int*, double*, double*);
void store(const string&, const int*, double*);
double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
void runge_kutta(shape shape,
                 gsl_matrix *gamma_tensor, gsl_matrix *inertia_tensor, gsl_matrix *dissipative_tensor,
                 gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double inertia_df_tensor[5][5][5], double starting[5], double step_length[5], double* stortartion);
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
    // inertia_tensor是惯性张量，dissipative_tenson是耗散张量
    gsl_matrix *inertia_tensor, *dissipative_tensor;
    inertia_tensor = gsl_matrix_alloc(5,5);
    dissipative_tensor = gsl_matrix_alloc(5,5);
    gsl_matrix *inertia_inverse, *gamma_tensor;
    inertia_inverse = gsl_matrix_alloc(5,5);
    gamma_tensor = gsl_matrix_alloc(5,5);
    // inertia_df_tensor是惯性张量的微分
    double inertia_df_tensor[5][5][5];
    gsl_vector *generalized_coordinates, *generalized_momenta;
    generalized_coordinates = gsl_vector_alloc(5);
    generalized_momenta = gsl_vector_alloc(5);
    for (size_t i=0; i<5; i++){
        gsl_vector_set(generalized_momenta, i, 0);
        gsl_vector_set(generalized_coordinates, i, gs(i)*step_length[i]+starting[i]);
    }
    string label="lrzcs";

    shape.set_l(gsl_vector_get(generalized_coordinates,0));
    shape.set_r(gsl_vector_get(generalized_coordinates,1));
    shape.set_z(gsl_vector_get(generalized_coordinates,2));
    shape.set_c(gsl_vector_get(generalized_coordinates,3));
    shape.set_s(gsl_vector_get(generalized_coordinates,4));
    shape.efficiency();

    for(size_t i=0;i<5;i++){
        for(size_t j=0;j<=5;j++){
            gsl_matrix_set(inertia_tensor,i,j,inertia(shape,label[i],label[j]));
            gsl_matrix_set(dissipative_tensor,i,j,dissipative(shape,label[i],label[j]));
            for(size_t k=0; k<5; k++){
                inertia_df_tensor[i][j][k] = inertia_df(shape, label[i], label[j], label[k]);
            }
        }
    }

    gsl_permutation * p = gsl_permutation_alloc (5);
    int signum;
    gsl_linalg_LU_decomp (inertia_tensor, p, &signum);
    gsl_linalg_LU_invert (inertia_tensor, p, inertia_inverse);

    gsl_vector *eval = gsl_vector_alloc (5);
    gsl_matrix *evec = gsl_matrix_alloc (5, 5);
    gsl_matrix *eigen = gsl_matrix_alloc (5, 5);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (5);
    gsl_eigen_symmv (dissipative_tensor, eval, evec, w);
    gsl_eigen_symmv_free (w);

    gsl_matrix *evec_invert = gsl_matrix_alloc (5, 5);
    for(size_t i=0; i<5; i++){
        gsl_matrix_set(eigen, i, i, sqrt(gsl_vector_get(eval, i)));
    }
    gsl_linalg_LU_decomp (evec, p, &signum);
    gsl_linalg_LU_invert (evec, p, evec_invert);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, evec_invert, eigen,
                    0.0, gamma_tensor);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, gamma_tensor, evec,
                    0.0, gamma_tensor);
    //这里gamma_tensor的结果就是dissipative_tensor的平方根
    runge_kutta(shape,
                gamma_tensor,inertia_inverse,dissipative_tensor,
                generalized_coordinates,generalized_momenta,
                inertia_df_tensor, starting, step_length, storation);

    return 0;
}
