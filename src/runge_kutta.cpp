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
#include <boost/array.hpp>
#include <cmath>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array<double, 10> state_type;

struct para{
    shape para_shape;
    gsl_matrix* gamma;
    gsl_matrix* inertia_inverse;
    gsl_matrix* dissiption;
    gsl_vector* hyper_interpo_df;
    double inertia_df[5][5][5];
    double starting[5];
    double step_length[5];
    double *storation;
    double gauss_random;
    unsigned long rng_seed = 0;
    };

para param;

double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
double hypercubic_interp_df(shape shape, double starting[], double step_length[], double storation[], const char label);
double hypercubic_interp(shape shape, double* starting, double* step_length, double* storation);
void gauss_rand(double u, unsigned long rng_seed);

void qp(const state_type y, state_type &dydt, double t){

    gsl_matrix *inertia_tensor, *dissipative_tensor;
    inertia_tensor = gsl_matrix_alloc(5,5);
    dissipative_tensor = gsl_matrix_alloc(5,5);

    gsl_matrix *inertia_inverse, *gamma_tensor;
    inertia_inverse = gsl_matrix_alloc(5,5);
    gamma_tensor = gsl_matrix_alloc(5,5);

    char label[6]="lrzcs";

    gsl_vector* hyperU_df;
    hyperU_df = gsl_vector_alloc(5);
    double hyperU;

    // inertia_df_tensor是惯性张量的微分
    double inertia_df_tensor[5][5][5];

    param.para_shape.set_l(y[0]);
    param.para_shape.set_r(y[1]);
    param.para_shape.set_z(y[2]);
    param.para_shape.set_c(y[3]);
    param.para_shape.set_s(y[4]);
    param.para_shape.efficiency();

    double *starting, *step_length;
    starting = param.starting;
    step_length = param.step_length;
    double *storation;
    storation = param.storation;
    cout<<"qp"<<endl;
    cout<<param.starting[0]<<' '<<param.starting[1]<<' '<<param.starting[2]<<' '<<param.starting[3]<<' '<<param.starting[4]<<' '<<endl;

    for(size_t i=0;i<5;i++){
        gsl_vector_set(hyperU_df,i,hypercubic_interp_df(param.para_shape, starting, step_length, storation, label[i]));
        for(size_t j=0;j<5;j++){
//            cout<<label[i]<<label[j]<<' '<<inertia(param.para_shape,label[i],label[j])<<endl;
            gsl_matrix_set(inertia_tensor,i,j,inertia(param.para_shape,label[i],label[j]));
            gsl_matrix_set(dissipative_tensor,i,j,dissipative(param.para_shape,label[i],label[j]));
            for(size_t k=0; k<5; k++){
                inertia_df_tensor[i][j][k] = inertia_df(param.para_shape, label[i], label[j], label[k]);
            }
        }
    }

    gsl_permutation * permu = gsl_permutation_alloc (5);
    int signum;
    cout<<"inertia_tensor"<<endl;
    for(size_t i=0; i<5; i++){
        cout<<i<<' '
            <<gsl_matrix_get(inertia_tensor,i,0)<<' '
            <<gsl_matrix_get(inertia_tensor,i,1)<<' '
            <<gsl_matrix_get(inertia_tensor,i,2)<<' '
            <<gsl_matrix_get(inertia_tensor,i,3)<<' '
            <<gsl_matrix_get(inertia_tensor,i,4)<<' '<<endl;
    }

    gsl_linalg_LU_decomp (inertia_tensor, permu, &signum);
    gsl_linalg_LU_invert (inertia_tensor, permu, inertia_inverse);

    //\frac{\partial M^-1}{\partial q_l} = M^-1*\frac{\partial M}{\partial q_l}M^-1
    gsl_matrix* inertia_df_tensor_tmp;
    inertia_df_tensor_tmp = gsl_matrix_alloc(5, 5);
    for (size_t i=0; i<5; i++){
        for(size_t j=0; j<5; j++){
            for (size_t k=0; k<5; k++){
                gsl_matrix_set(inertia_df_tensor_tmp,j,k,inertia_df_tensor[j][k][i]);
            }
        }
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, inertia_inverse, inertia_df_tensor_tmp,
                        0.0, inertia_df_tensor_tmp);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, inertia_df_tensor_tmp, inertia_inverse,
                        0.0, inertia_df_tensor_tmp);
        for(size_t j=0; j<5; j++){
            for (size_t k=0; k<5; k++){
                inertia_df_tensor[j][k][i] = gsl_matrix_get(inertia_df_tensor_tmp,j,k);
            }
        }
    }
    gsl_matrix_free(inertia_df_tensor_tmp);

    hyperU = hypercubic_interp(param.para_shape, starting, step_length, storation);
    double temperature = sqrt(hyperU);

    gsl_vector *eval = gsl_vector_alloc (5);
    gsl_matrix *evec = gsl_matrix_alloc (5, 5);
    gsl_matrix *eigen = gsl_matrix_alloc (5, 5);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (5);
    gsl_eigen_symmv (dissipative_tensor, eval, evec, w);
    gsl_eigen_symmv_free (w);

    gsl_matrix *evec_invert = gsl_matrix_alloc (5, 5);
    for(size_t i=0; i<5; i++){
        gsl_matrix_set(eigen, i, i, sqrt(gsl_vector_get(eval, i)/temperature));
    }
    gsl_linalg_LU_decomp (evec, permu, &signum);
    gsl_linalg_LU_invert (evec, permu, evec_invert);
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
    gsl_permutation_free(permu);

    param.gamma = gamma_tensor;
    param.inertia_inverse = inertia_inverse;
    param.dissiption = dissipative_tensor;
    param.inertia_df[5][5][5] = inertia_df_tensor[5][5][5];
    param.hyper_interpo_df = hyperU_df;

    gauss_rand(param.gauss_random, param.rng_seed);

    double tmp1, tmp2;
    for(size_t i=0; i<5; i++){
        tmp1 = 0.;
        tmp2 = 0.;
        for(size_t j=0; j<5; j++){
            tmp1 += gsl_matrix_get(param.inertia_inverse,i,j)*y[j+5];
            for(size_t k=0; k<5; k++){
                tmp2 += -gsl_vector_get(param.hyper_interpo_df,i)
                        -1./2.*(param.inertia_df[j][k][i])*y[j+5]*y[k+5]
                        -gsl_matrix_get(param.dissiption,i,j)*gsl_matrix_get(param.inertia_inverse,j,k)*y[k+5]
                        +gsl_matrix_get(param.gamma,i,j)*param.gauss_random;
            }
        }
        dydt[i] = tmp1;
        dydt[i+5] = tmp2;
    }

    gsl_matrix_free(inertia_tensor);
    gsl_matrix_free(dissipative_tensor);
    gsl_matrix_free(inertia_inverse);
    gsl_matrix_free(gamma_tensor);
    gsl_vector_free(hyperU_df);
}
void runge_kutta(gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double starting[5], double step_length[5], int steps[5], double storation[],
                 shape shape){

    copy(starting,starting+5,param.starting);
    copy(step_length,step_length+5,param.step_length);
    param.storation = storation;
    param.para_shape = shape;

    cout<<"rugge kutta"<<endl;
    cout<<starting[0]<<' '<<starting[1]<<' '<<starting[2]<<' '<<endl;
    cout<<param.starting[0]<<' '<<param.starting[1]<<' '<<param.starting[2]<<' '<<endl;
    cout<<step_length[0]<<' '<<step_length[1]<<' '<<step_length[2]<<' '<<endl;
    cout<<param.step_length[0]<<' '<<param.step_length[1]<<' '<<param.step_length[2]<<' '<<endl;
    cout<<param.storation[0]<<' '<<param.storation[1]<<' '<<endl;

    state_type x;
    runge_kutta4< state_type > stepper;
    const double dt=1.e-25;
    // initialize the position and velocity
    for(size_t i=0; i<5; i++){
        x[i] = gsl_vector_get(generalized_coordinates,i);
        x[i+5] = gsl_vector_get(generalized_momenta,i);
    }

    for( double t=0.0 ; t<1.e-20 ; t+= dt ) {
        stepper.do_step(qp, x, t, dt);
        if (2*x[0]>11*x[1])
            break;
    }
    for(size_t i=0; i<5; i++){
        gsl_vector_set(generalized_coordinates,i,x[i]);
    }
}