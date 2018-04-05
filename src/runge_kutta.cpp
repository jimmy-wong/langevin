#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "../include/shape.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array<double, 10> state_type;

struct para{
    shape para_shape;
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
void gauss_rand(double* u, unsigned long* rng_seed);

void qp(const state_type y, state_type &dydt, double t){

    gsl_matrix *inertia_tensor, *dissipative_tensor;
    inertia_tensor = gsl_matrix_alloc(5,5);
    dissipative_tensor = gsl_matrix_alloc(5,5);

    gsl_matrix *inertia_inverse, *gamma_tensor, *gamma_tensor_tmp;
    inertia_inverse = gsl_matrix_alloc(5,5);
    gamma_tensor = gsl_matrix_alloc(5,5);
    gamma_tensor_tmp = gsl_matrix_alloc(5,5);

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

    cout<<"current position: "
        <<param.para_shape.get_l()<<' '
        <<param.para_shape.get_r()<<' '
        <<param.para_shape.get_z()<<' '
        <<param.para_shape.get_c()<<' '
        <<param.para_shape.get_s()<<' '<<endl;

    double *starting, *step_length;
    starting = param.starting;
    step_length = param.step_length;
    double *storation;
    storation = param.storation;
    gsl_matrix* inertia_tensor_tmp = gsl_matrix_alloc(5,5);
    gsl_matrix* dissipative_tensor_tmp = gsl_matrix_alloc(5,5);

    for(size_t i=0;i<5;i++){
        gsl_vector_set(hyperU_df,i,hypercubic_interp_df(param.para_shape, starting, step_length, storation, label[i]));
        for(size_t j=0;j<=i;j++){
            double tmp = inertia(param.para_shape,label[i],label[j]);
            gsl_matrix_set(inertia_tensor,i,j,tmp);
            gsl_matrix_set(inertia_tensor,j,i,tmp);

            gsl_matrix_set(inertia_tensor_tmp,i,j,tmp);
            gsl_matrix_set(inertia_tensor_tmp,j,i,tmp);

            tmp = dissipative(param.para_shape,label[i],label[j]);
            gsl_matrix_set(dissipative_tensor,i,j,tmp);
            gsl_matrix_set(dissipative_tensor,j,i,tmp);

            gsl_matrix_set(dissipative_tensor_tmp,i,j,tmp);
            gsl_matrix_set(dissipative_tensor_tmp,j,i,tmp);

            for(size_t k=0; k<5; k++){
                tmp = inertia_df(param.para_shape, label[i], label[j], label[k]);
                inertia_df_tensor[i][j][k] = tmp;
                inertia_df_tensor[j][i][k] = tmp;
            }
        }
    }

//    cout << "dissipative tensor" << endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(dissipative_tensor, i, 0) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 1) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 2) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 3) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 4) << ' ' << endl;
//    }

    gsl_permutation * permu = gsl_permutation_alloc (5);
    int signum;

    gsl_linalg_LU_decomp (inertia_tensor_tmp, permu, &signum);
    gsl_linalg_LU_invert (inertia_tensor_tmp, permu, inertia_inverse);
    gsl_matrix_free(inertia_tensor_tmp);

    //\frac{\partial M^-1}{\partial q_l} = -M^-1*\frac{\partial M}{\partial q_l}M^-1
    gsl_matrix* inertia_df_tensor_tmp1, *inertia_df_tensor_tmp2;
    inertia_df_tensor_tmp1 = gsl_matrix_alloc(5, 5);
    inertia_df_tensor_tmp2 = gsl_matrix_alloc(5, 5);
    for (size_t i=0; i<5; i++){
        for(size_t j=0; j<5; j++){
            for (size_t k=0; k<5; k++){
                gsl_matrix_set(inertia_df_tensor_tmp1,j,k,-inertia_df_tensor[j][k][i]);
            }
        }
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, inertia_inverse, inertia_df_tensor_tmp1,
                        0.0, inertia_df_tensor_tmp2);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, inertia_df_tensor_tmp2, inertia_inverse,
                        0.0, inertia_df_tensor_tmp1);
        for(size_t j=0; j<5; j++){
            for (size_t k=0; k<5; k++){
                inertia_df_tensor[j][k][i] = gsl_matrix_get(inertia_df_tensor_tmp1,j,k);
            }
        }
    }
    gsl_matrix_free(inertia_df_tensor_tmp1);
    gsl_matrix_free(inertia_df_tensor_tmp2);

    double kinetic_energy;
    gsl_vector* momenta = gsl_vector_alloc(5);
    gsl_vector* momenta_tmp = gsl_vector_alloc(5);
    for(size_t i=0; i<5; i++){
        gsl_vector_set(momenta, i, y[i+5]);
//        cout<<gsl_vector_get(momenta,i)<<endl;
    }
    gsl_blas_dgemv (CblasNoTrans,
                    1.0, inertia_inverse, momenta,
                    0.0, momenta_tmp);
    gsl_blas_ddot(momenta, momenta_tmp, &kinetic_energy);
    kinetic_energy = kinetic_energy/2.;
//    cout<<"momenta="
//        <<gsl_vector_get(momenta,0)<<' '
//        <<gsl_vector_get(momenta,1)<<' '
//        <<gsl_vector_get(momenta,2)<<' '
//        <<gsl_vector_get(momenta,3)<<' '
//        <<gsl_vector_get(momenta,4)<<' '<<endl;
//    cout<<"momenta="<<y[5]<<' '<<y[6]<<' '<<y[7]<<' '<<y[8]<<' '<<y[9]<<' '<<endl;

    gsl_vector_free(momenta);
    gsl_vector_free(momenta_tmp);

    //E_int = E_x - 1/2 M^-1*p*p - U
    hyperU = hypercubic_interp(param.para_shape, starting, step_length, storation);
    double temperature = sqrt((param.para_shape.get_excited_energy()-hyperU-kinetic_energy)/
                                      param.para_shape.get_level_density());

    //gamma tensor: gamma_ij*gamma_jk = T^* D_ik, D_ik is dissipative tensor
    //here we use approximation that T^* = T
    gsl_vector *eval = gsl_vector_alloc (5);
    gsl_matrix *evec = gsl_matrix_alloc (5, 5);
    gsl_matrix *eigen = gsl_matrix_alloc (5, 5);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (5);
    gsl_eigen_symmv (dissipative_tensor_tmp, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_matrix_free(dissipative_tensor_tmp);

//    cout << "dissipative tensor" << endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(dissipative_tensor, i, 0) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 1) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 2) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 3) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 4) << ' ' << endl;
//    }

    gsl_matrix* evec_tmp = gsl_matrix_alloc(5,5);
    for(size_t i=0; i<5; i++){
        for(size_t j=0; j<5; j++){
            gsl_matrix_set(evec_tmp,i,j,gsl_matrix_get(evec,i,j));
        }
    }

//    cout<<"eigen value"<<endl;
    cout<<"temperature="<<temperature<<", "
        <<"U="<<hyperU<<", "
        <<"kinetic energy="<<kinetic_energy<< endl;
//    cout << gsl_vector_get(eval, 0) << ' '
//         << gsl_vector_get(eval, 1) << ' '
//         << gsl_vector_get(eval, 2) << ' '
//         << gsl_vector_get(eval, 3) << ' '
//         << gsl_vector_get(eval, 4) << ' ' << endl;

    gsl_matrix *evec_invert = gsl_matrix_alloc (5, 5);
    for(size_t i=0; i<5; i++){
        if (abs(gsl_vector_get(eval, i))<1e-8){
            gsl_matrix_set(eigen, i, i, 0.);
        }
        else{
            gsl_matrix_set(eigen, i, i, sqrt(abs(gsl_vector_get(eval, i))/temperature));
        }
    }
    gsl_linalg_LU_decomp (evec_tmp, permu, &signum);
    gsl_linalg_LU_invert (evec_tmp, permu, evec_invert);
    gsl_matrix_free(evec_tmp);
//    gsl_matrix* matrix_tmp = gsl_matrix_alloc(5,5);
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
//                    1.0, evec, evec_invert,
//                    0.0, matrix_tmp);

//    cout<<"evec invert"<<endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(evec_invert, i, 0) << ' '
//             << gsl_matrix_get(evec_invert, i, 1) << ' '
//             << gsl_matrix_get(evec_invert, i, 2) << ' '
//             << gsl_matrix_get(evec_invert, i, 3) << ' '
//             << gsl_matrix_get(evec_invert, i, 4) << ' ' << endl;
//    }
//    gsl_matrix_free(matrix_tmp);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, evec_invert, eigen,
                    0.0, gamma_tensor_tmp);
//    cout<<"evec invert"<<endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(evec_invert, i, 0) << ' '
//             << gsl_matrix_get(evec_invert, i, 1) << ' '
//             << gsl_matrix_get(evec_invert, i, 2) << ' '
//             << gsl_matrix_get(evec_invert, i, 3) << ' '
//             << gsl_matrix_get(evec_invert, i, 4) << ' ' << endl;
//    }
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, gamma_tensor_tmp, evec,
                    0.0, gamma_tensor);
    gsl_matrix_free(gamma_tensor_tmp);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(eigen);
    gsl_matrix_free(evec_invert);
    gsl_permutation_free(permu);

//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    gsl_vector* random_vector = gsl_vector_alloc(5);
    for(size_t i=0; i<5; i++){
        gauss_rand(&(param.gauss_random), &(param.rng_seed));
        gsl_vector_set(random_vector,i,param.gauss_random);
    }
//    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;

    double tmp1, tmp2;
//    cout<<"dydt ";
    for(size_t i=0; i<5; i++){
        tmp1 = 0.;
        tmp2 = 0.;
        for(size_t j=0; j<5; j++){
            tmp1 += gsl_matrix_get(inertia_inverse,i,j)*y[j+5];
            for(size_t k=0; k<5; k++){
                tmp2 += -1./2.*(inertia_df_tensor[j][k][i])*y[j+5]*y[k+5]
                        -gsl_matrix_get(dissipative_tensor,i,j)*gsl_matrix_get(inertia_inverse,j,k)*y[k+5];
            }
            tmp2 += gsl_matrix_get(gamma_tensor,i,j)*gsl_vector_get(random_vector,j);
        }
        tmp2 += -gsl_vector_get(hyperU_df,i);
        dydt[i] = tmp1;
        dydt[i+5] = tmp2;
//        cout<<dydt[i]<<' '<<dydt[i+5]<<' ';
    }
//    cout<<endl;
//    cout<<"inertia inverse tensor"<<endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(inertia_inverse, i, 0) << ' '
//             << gsl_matrix_get(inertia_inverse, i, 1) << ' '
//             << gsl_matrix_get(inertia_inverse, i, 2) << ' '
//             << gsl_matrix_get(inertia_inverse, i, 3) << ' '
//             << gsl_matrix_get(inertia_inverse, i, 4) << ' ' << endl;
//    }
//    cout << "dissipative tensor" << endl;
//    for(size_t i=0; i<5; i++) {
//        cout << gsl_matrix_get(dissipative_tensor, i, 0) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 1) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 2) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 3) << ' '
//             << gsl_matrix_get(dissipative_tensor, i, 4) << ' ' << endl;
//    }
//    cout<<"gamma tensor"<<endl;
//    for(size_t i=0; i<5; i++) {
//        cout<<gsl_matrix_get(gamma_tensor,i,0)<<' '
//            <<gsl_matrix_get(gamma_tensor,i,1)<<' '
//            <<gsl_matrix_get(gamma_tensor,i,2)<<' '
//            <<gsl_matrix_get(gamma_tensor,i,3)<<' '
//            <<gsl_matrix_get(gamma_tensor,i,4)<<' '<<endl;
//    }
    gsl_matrix_free(inertia_tensor);
    gsl_matrix_free(dissipative_tensor);
    gsl_matrix_free(inertia_inverse);
    gsl_matrix_free(gamma_tensor);
    gsl_vector_free(hyperU_df);
}
void runge_kutta(gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double starting[5], double step_length[5], double storation[],
                 shape shape){

    copy(starting,starting+5,param.starting);
    copy(step_length,step_length+5,param.step_length);
    param.storation = storation;
    param.para_shape = shape;

    state_type x;
    runge_kutta4< state_type > stepper;
    const double dt=2.5e-4;
    // initialize the position and velocity
    for(size_t i=0; i<5; i++){
        x[i] = gsl_vector_get(generalized_coordinates,i);
        x[i+5] = gsl_vector_get(generalized_momenta,i);
    }

    for( double t=0.0 ; t<1.e1 ; t+= dt ) {
        stepper.do_step(qp, x, t, dt);
        if (2*x[0]>11*x[1])
            break;
    }
    for(size_t i=0; i<5; i++){
        gsl_vector_set(generalized_coordinates,i,x[i]);
    }
}