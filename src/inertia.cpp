#include "../include/shape.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <cmath>

double inertia(shape shape, double density, const char label_i, const char label_j){
    // density means nuclear density, a constant
    // label i&j represents M_ij
    double l = shape.get_l();
    double s = shape.get_s();
    double x, w;
    size_t n = 28;
    double result = 0;
    gsl_integration_glfixed_table *table = NULL;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++)
    {
        gsl_integration_glfixed_point(-1., 1., i, &x, &w, table);
        x = l*x + s;
        result += w*(pow(shape.Rho(x),2)/8.*A_derivative(shape, l-s,x,label_i)*A_derivative(shape,l-s,x,label_j)
                     +shape.Rho(x)*A_para(shape, l-s,x,label_i)*A_para(shape, l-s,x,label_j));
    }
    gsl_integration_glfixed_table_free(table);
    result = result*M_PI*density;
    return result;
}
double inertia_df(shape shape, double density, const char label_i, const char label_j){

}