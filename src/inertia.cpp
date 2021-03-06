#include "../include/shape.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <cmath>

double inertia(shape shape, const char label_i, const char label_j){
    // density means nuclear density, a constant
    // label i&j represents M_ij
    double l = shape.get_l();
    double s = shape.get_s();
    double density = shape.get_density();
    double x, w;
    size_t n = 28;
    double result = 0;
    gsl_integration_glfixed_table *table = nullptr;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++)
    {
        gsl_integration_glfixed_point(-l-s, l-s, i, &x, &w, table);
        result += w*(pow(shape.Rho(x),2)/8.*ADerivativeZ(shape, x, -l-s,label_i)*ADerivativeZ(shape, x, -l-s, label_j)
                     +shape.Rho(x)*A(shape, x, -l-s, label_i)*A(shape, x, -l-s, label_j));
    }
    gsl_integration_glfixed_table_free(table);
    result = result*(l*shape.get_Rcn())*M_PI*density;
    return result;
}
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l){
    double l = shape.get_l();
    double s = shape.get_s();
    double density = shape.get_density();
    double x, w;
    size_t n = 28;
    double result = 0;
    gsl_integration_glfixed_table *table = nullptr;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++)
    {
        gsl_integration_glfixed_point(-l-s, l-s, i, &x, &w, table);
        result += w*(shape.Rho(x)/4.*shape.RhoDerivative(x,label_l)*ADerivativeZ(shape, x, -l-s, label_i)*ADerivativeZ(shape, x, -l-s, label_j)+
                     pow(shape.Rho(x),2)/8.*ADDerivativeZQ(shape, x, -l-s, label_i, label_l)*ADerivativeZ(shape, x, -l-s, label_j)+
                     pow(shape.Rho(x),2)/8.*ADerivativeZ(shape, x, -l-s, label_i)*ADDerivativeZQ(shape, x, -l-s, label_j, label_l)+
                     shape.RhoDerivative(x, label_l)*A(shape, x, -l-s, label_i)*A(shape, x, -l-s, label_j)+
                     shape.Rho(x)*ADerivativeQ(shape, x, -l-s, label_i, label_l)*A(shape, x, -l-s, label_j)+
                     shape.Rho(x)*A(shape, x, -l-s, label_i)*ADerivativeQ(shape, x, -l-s, label_j, label_l));
    }
    gsl_integration_glfixed_table_free(table);
    result = result*(l*shape.get_Rcn())*M_PI*density;
    return result;
}