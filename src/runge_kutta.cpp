#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "../include/shape.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>

using namespace std;
using boost::format;
//using boost::io::group;
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
int counter;

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
//    vector<gsl_matrix>
    double inertia_df_tensor[5][5][5];

    param.para_shape.set_l(y[0]/param.para_shape.get_Rcn());
    param.para_shape.set_r(y[1]/param.para_shape.get_Rcn());
    param.para_shape.set_z(y[2]/param.para_shape.get_Rcn());
    param.para_shape.set_c(y[3]/param.para_shape.get_Rcn());
    param.para_shape.set_s(y[4]/param.para_shape.get_Rcn());
    param.para_shape.coefficiency();

    cout<<"counter: "<<counter<<endl;
    counter += 1;
    cout<<format("current position: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
          %param.para_shape.get_l()
          %param.para_shape.get_r()
          %param.para_shape.get_z()
          %param.para_shape.get_c()
          %param.para_shape.get_s();

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
//    for(size_t i=0; i<5; i++) {
//        cout << format("inertia tensor: %1$+15.3f %2$+15.3f %3$+15.3f %4$+15.3f %5$+15.3f\n")
//                % gsl_matrix_get(inertia_tensor, i, 0)
//                % gsl_matrix_get(inertia_tensor, i, 1)
//                % gsl_matrix_get(inertia_tensor, i, 2)
//                % gsl_matrix_get(inertia_tensor, i, 3)
//                % gsl_matrix_get(inertia_tensor, i, 4);
//    }
//    cout<<"inertia tensor df 1"<<endl;
//    for(size_t i=0; i<5; i++){
//        cout<<label[i]<<endl;
//        for(size_t j=0; j<5; j++){
//            cout << format("inertia tensor df: %1$+15.3f %2$+15.3f %3$+15.3f %4$+15.3f %5$+15.3f\n")
//                    % inertia_df_tensor[j][0][i]
//                    % inertia_df_tensor[j][1][i]
//                    % inertia_df_tensor[j][2][i]
//                    % inertia_df_tensor[j][3][i]
//                    % inertia_df_tensor[j][4][i];
//        }
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
//    cout<<"inertia tensor df 2"<<endl;
//    for(size_t i=0; i<5; i++){
//        cout<<label[i]<<endl;
//        for(size_t j=0; j<5; j++){
//            cout << format("inertia tensor df: %1$+15.10f %2$+15.10f %3$+15.10f %4$+15.10f %5$+15.10f\n")
//                    % inertia_df_tensor[j][0][i]
//                    % inertia_df_tensor[j][1][i]
//                    % inertia_df_tensor[j][2][i]
//                    % inertia_df_tensor[j][3][i]
//                    % inertia_df_tensor[j][4][i];
//        }
//    }

    gsl_matrix_free(inertia_df_tensor_tmp1);
    gsl_matrix_free(inertia_df_tensor_tmp2);

    double kinetic_energy;
    gsl_vector* momenta = gsl_vector_alloc(5);
    gsl_vector* momenta_tmp = gsl_vector_alloc(5);
    for(size_t i=0; i<5; i++){
        gsl_vector_set(momenta, i, y[i+5]);
    }
    gsl_blas_dgemv (CblasNoTrans,
                    1.0, inertia_inverse, momenta,
                    0.0, momenta_tmp);
    gsl_blas_ddot(momenta, momenta_tmp, &kinetic_energy);
    kinetic_energy = kinetic_energy/2.;
//    gsl_blas_dgemv(CblasNoTrans, 1.0, inertia_inverse, momenta,  0.0, momenta_tmp);
//    cout<<format("current dp/dt: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//          %gsl_vector_get(momenta_tmp,0)
//          %gsl_vector_get(momenta_tmp,1)
//          %gsl_vector_get(momenta_tmp,2)
//          %gsl_vector_get(momenta_tmp,3)
//          %gsl_vector_get(momenta_tmp,4);
//
    cout<<format("current momentum: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
          %gsl_vector_get(momenta,0)
          %gsl_vector_get(momenta,1)
          %gsl_vector_get(momenta,2)
          %gsl_vector_get(momenta,3)
          %gsl_vector_get(momenta,4);

    gsl_vector_free(momenta);
    gsl_vector_free(momenta_tmp);

    //E_int = E_x - 1/2 M^-1*p*p - U
    hyperU = hypercubic_interp(param.para_shape, starting, step_length, storation);
    double temperature = param.para_shape.get_excited_energy()-hyperU-kinetic_energy>=0?
                         sqrt((param.para_shape.get_excited_energy()-hyperU-kinetic_energy)/param.para_shape.get_level_density()) :0.002;

    //gamma tensor: gamma_ij*gamma_jk = T^* D_ik, D_ik is dissipative tensor
    //here we use approximation that T^* = T
    gsl_vector *eval = gsl_vector_alloc (5);
    gsl_matrix *evec = gsl_matrix_alloc (5, 5);
    gsl_matrix *eigen = gsl_matrix_alloc (5, 5);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (5);
    gsl_eigen_symmv (dissipative_tensor_tmp, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_matrix_free(dissipative_tensor_tmp);

    gsl_matrix* evec_tmp = gsl_matrix_alloc(5,5);
    for(size_t i=0; i<5; i++){
        for(size_t j=0; j<5; j++){
            gsl_matrix_set(evec_tmp,i,j,gsl_matrix_get(evec,i,j));
        }
    }
    cout<<format("temperature: %1$7.3f, U: %2$7.3f, kinetic energy: %3$7.3f\n")
          %temperature
          %hyperU
          %kinetic_energy;

    gsl_matrix *evec_invert = gsl_matrix_alloc (5, 5);
    gsl_matrix_set_zero(eigen);
    for(size_t i=0; i<5; i++){
        if (abs(gsl_vector_get(eval, i))<1e-8){
            gsl_matrix_set(eigen, i, i, 0.);
        }
        else{
            gsl_matrix_set(eigen, i, i, sqrt(abs(gsl_vector_get(eval, i))));
        }
    }
//    for(size_t i=0; i<5; i++) {
//        cout << format("evec tensor: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//                % gsl_matrix_get(evec, i, 0)
//                % gsl_matrix_get(evec, i, 1)
//                % gsl_matrix_get(evec, i, 2)
//                % gsl_matrix_get(evec, i, 3)
//                % gsl_matrix_get(evec, i, 4);
//    }
//    for(size_t i=0; i<5; i++) {
//        cout << format("eigen tensor: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//                % gsl_matrix_get(eigen, i, 0)
//                % gsl_matrix_get(eigen, i, 1)
//                % gsl_matrix_get(eigen, i, 2)
//                % gsl_matrix_get(eigen, i, 3)
//                % gsl_matrix_get(eigen, i, 4);
//    }

    gsl_linalg_LU_decomp (evec_tmp, permu, &signum);
    gsl_linalg_LU_invert (evec_tmp, permu, evec_invert);
    gsl_matrix_free(evec_tmp);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, evec, eigen,
                    0.0, gamma_tensor_tmp);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, gamma_tensor_tmp, evec_invert,
                    0.0, gamma_tensor);
//    cout << format("eigen: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//            % gsl_vector_get(eval, 0)
//            % gsl_vector_get(eval, 1)
//            % gsl_vector_get(eval, 2)
//            % gsl_vector_get(eval, 3)
//            % gsl_vector_get(eval, 4);
//
//    for(size_t i=0; i<5; i++) {
//        cout << format("gamma tensor: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//                % gsl_matrix_get(gamma_tensor, i, 0)
//                % gsl_matrix_get(gamma_tensor, i, 1)
//                % gsl_matrix_get(gamma_tensor, i, 2)
//                % gsl_matrix_get(gamma_tensor, i, 3)
//                % gsl_matrix_get(gamma_tensor, i, 4);
//    }
//    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,gamma_tensor,gamma_tensor,0.0,gamma_tensor_tmp);
//    for(size_t i=0; i<5; i++) {
//        cout << format("dissipative tensor: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//                % gsl_matrix_get(gamma_tensor_tmp, i, 0)
//                % gsl_matrix_get(gamma_tensor_tmp, i, 1)
//                % gsl_matrix_get(gamma_tensor_tmp, i, 2)
//                % gsl_matrix_get(gamma_tensor_tmp, i, 3)
//                % gsl_matrix_get(gamma_tensor_tmp, i, 4);
//    }
//    for(size_t i=0; i<5; i++) {
//        cout << format("dissipative tensor: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
//                % gsl_matrix_get(dissipative_tensor, i, 0)
//                % gsl_matrix_get(dissipative_tensor, i, 1)
//                % gsl_matrix_get(dissipative_tensor, i, 2)
//                % gsl_matrix_get(dissipative_tensor, i, 3)
//                % gsl_matrix_get(dissipative_tensor, i, 4);
//    }

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
    for(size_t i=0; i<5; i++){
        tmp1 = 0.;
        tmp2 = 0.;
        for(size_t j=0; j<5; j++){
            tmp1 += gsl_matrix_get(inertia_inverse,i,j)*y[j+5];
            for(size_t k=0; k<5; k++){
                tmp2 += -1./2.*(inertia_df_tensor[j][k][i])*y[j+5]*y[k+5]
                        -gsl_matrix_get(dissipative_tensor,i,j)*gsl_matrix_get(inertia_inverse,j,k)*y[k+5];
            }
            tmp2 += gsl_matrix_get(gamma_tensor,i,j)*sqrt(2*temperature/1.e-0)*gsl_vector_get(random_vector,j);
        }
        tmp2 += -gsl_vector_get(hyperU_df,i);
        dydt[i] = tmp1;
        dydt[i+5] = tmp2;
    }
//    double Utmp1, Utmp2;
//    param.para_shape.set_c(param.para_shape.get_c()-0.00001);
//    Utmp1 = hypercubic_interp(param.para_shape, starting, step_length, storation);
//    param.para_shape.set_c(param.para_shape.get_c()+0.00002);
//    Utmp2 = hypercubic_interp(param.para_shape, starting, step_length, storation);
//    cout<<(Utmp1-Utmp2)/(-0.00002*param.para_shape.get_Rcn())<<endl;
//    param.para_shape.set_c(param.para_shape.get_c()-0.00001);
//    cout<<format("dqdt: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n") %dydt[5] %dydt[6] %dydt[7] %dydt[8] %dydt[9];
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
    euler< state_type > stepper;
    const double dt=1.e-0;
    // initialize the position and velocity
    for(size_t i=0; i<5; i++){
        x[i] = gsl_vector_get(generalized_coordinates,i)*shape.get_Rcn();
        x[i+5] = gsl_vector_get(generalized_momenta,i);
    }

    for(double t=0.0; t<1.e4; t+= dt) {
        stepper.do_step(qp, x, t, dt);
        if (2*x[0]>11*x[1]){
            counter = 0;
            break;
        }
    }
    counter = 0;
    for(size_t i=0; i<5; i++){
        gsl_vector_set(generalized_coordinates,i,x[i]);
    }
}