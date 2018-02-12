#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

double random(double rand);
void runge_kutta(gsl_matrix gamma_tensor, gsl_matrix inertia_tensor, gsl_matrix dissipative_tensor,
                 gsl_vector generalized_coordinates, gsl_vector generalized_momenta){
    double timedelta;
    
}