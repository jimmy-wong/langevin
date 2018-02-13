#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "../include/shape.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/odeint.hpp>
#include <cmath>

using namespace std;

struct para{
    gsl_matrix* gamma;
    gsl_matrix* inertia_inverse;
    gsl_matrix* dissiption;
    gsl_vector* hyper_interpo_df;
    double inertia_df[5][5][5];
    double hyper_interpo;
};

double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
double hypercubic_interp_df(shape shape, double* starting, double* step_length, double* storation, const char label);
double hypercubic_interp(shape shape, double* starting, double* step_length, double* storation);
double random(double rand);
int (* function) (double t, const double y[], double dydt[], void *params);
int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params);
int qp(double t, const double y[], double dydt[], void *params){
    struct para *p = (struct para *) params;

    double tmp1, tmp2;
    for(size_t i=0; i<5; i++){
        tmp1 = 0.;
        tmp2 = 0.;
        for(size_t j=0; j<5; j++){
            tmp1 += gsl_matrix_get(p->inertia_inverse,i,j)*y[j+5];
            for(size_t k=0; k<5; k++){
                tmp2 += -gsl_vector_get(p->hyper_interpo_df,i)
                        -1./2.*(p->inertia_df[j][k][i])*y[j+5]*y[k+5]
                        -gsl_matrix_get(p->dissiption,i,j)*gsl_matrix_get(p->inertia_inverse,j,k)*y[k+5]
                        +gsl_matrix_get(p->gamma,i,j)*random(10.);
            }
        }
        dydt[i] = tmp1;
        dydt[i+5] = tmp2;
    }
    return GSL_SUCCESS;
}
void runge_kutta(shape shape,
                 gsl_matrix *gamma_tensor, gsl_matrix *inertia_inverse, gsl_matrix *dissipative_tensor,
                 gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double inertia_df_tensor[5][5][5], double starting[5], double step_length[5], double* storation){
    double timedelta;
    para* param;
    string label="lrzcs";
    gsl_matrix *inertia_tensor;
    inertia_tensor = gsl_matrix_alloc(5,5);
    gsl_vector* hyperU_df;
    hyperU_df = gsl_vector_alloc(5);

    while (2*gsl_vector_get(generalized_coordinates,0)>11*gsl_vector_get(generalized_coordinates,1)){
        shape.set_l(gsl_vector_get(generalized_coordinates,0));
        shape.set_r(gsl_vector_get(generalized_coordinates,1));
        shape.set_z(gsl_vector_get(generalized_coordinates,2));
        shape.set_c(gsl_vector_get(generalized_coordinates,3));
        shape.set_s(gsl_vector_get(generalized_coordinates,4));
        shape.efficiency();

        for(size_t i=0;i<5;i++){
            gsl_vector_set(hyperU_df,i,hypercubic_interp_df(shape, starting, step_length, storation, label[i]));
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
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_matrix_free(eigen);
        gsl_matrix_free(evec_invert);
        gsl_permutation_free(p);
        param->gamma = gamma_tensor;
        param->inertia_inverse = inertia_inverse;
        param->dissiption = dissipative_tensor;
        param->inertia_df[5][5][5] = inertia_df_tensor[5][5][5];
        param->hyper_interpo = hypercubic_interp(shape, starting, step_length, storation);
        param->hyper_interpo_df = hyperU_df;
    }
}