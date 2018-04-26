#include <string>
#include <gsl/gsl_integration.h>
#include "../include/shape.h"
#include <gsl/gsl_min.h>

double dissipative_wall(shape shape, char label_i, char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double x, w;
    size_t n = 28;
    double tmp_result1 = 0, tmp_result2 = 0, result = 0;
    gsl_integration_glfixed_table *table = nullptr;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++) {
        gsl_integration_glfixed_point(-l - s, z - s, i, &x, &w, table);
        tmp_result1 += (shape.RhoDerivative(x, label_i) * shape.RhoDerivative(x, label_j) /
                   sqrt(4.*shape.Rho(x)+pow(shape.RhoDerivative(x,'x'),2.))) * w;
    }
    tmp_result1 = tmp_result1*((z+l)/2*shape.get_Rcn());
    for(size_t i=0; i<n; i++){
        gsl_integration_glfixed_point(z - s, l - s, i, &x, &w, table);
        tmp_result2 += (shape.RhoDerivative(x,label_i)*shape.RhoDerivative(x,label_j)/
                   sqrt(4*shape.Rho(x)+pow(shape.RhoDerivative(x,'x'),2.)))*w;
    }
    tmp_result2 = tmp_result2*((l-z)/2*shape.get_Rcn());
    gsl_integration_glfixed_table_free(table);
    result = (tmp_result2+tmp_result1)*M_PI*shape.get_density()*shape.get_average_v()*3./4.;
    return result;
}
double dissipative_wall2(shape shape, char label_i, char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double x, w;
    size_t n = 28;
    double tmp_result1 = 0, tmp_result2 = 0, result = 0;
    gsl_integration_glfixed_table *table = nullptr;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++) {
        gsl_integration_glfixed_point(-l-s, z-s, i, &x, &w, table);
        tmp_result1 += ((shape.RhoDerivative(x, label_i) +
                         shape.RhoDerivative(x, 'x') * shape.CenterOfMassDerivative('L', label_i)) *
                        (shape.RhoDerivative(x, label_j) +
                         shape.RhoDerivative(x, 'x') * shape.CenterOfMassDerivative('L', label_j)) /
                        sqrt(4. * shape.Rho(x) + pow(shape.RhoDerivative(x, 'x'), 2.))) * w;
    }
    tmp_result1 = tmp_result1*(l+z)/2*shape.get_Rcn();
    for(size_t i=0; i<n; i++) {
        gsl_integration_glfixed_point(z-s, l-s, i, &x, &w, table);
        tmp_result2 += ((shape.RhoDerivative(x,label_i)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('R',label_i))*
                   (shape.RhoDerivative(x,label_j)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('R',label_j))/
                   sqrt(4.*shape.Rho(x)+pow(shape.RhoDerivative(x,'x'),2.)))*w;
    }
    tmp_result2 = tmp_result2*(l-z)/2*shape.get_Rcn();
    gsl_integration_glfixed_table_free(table);
    result = (tmp_result1+tmp_result2)*M_PI*shape.get_average_v()/2.;
    return result;
}
double dissipative_window(shape shape, const char label_i, const char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double sigma = M_PI*pow(shape.get_r()*shape.get_Rcn(),2);
    double result;
    result = sigma*((shape.CenterOfMassDerivative('R',label_i)-shape.CenterOfMassDerivative('L',label_i))*
            (shape.CenterOfMassDerivative('R',label_j)-shape.CenterOfMassDerivative('L',label_j)))+
            32./9./sigma*(M_PI*(-1)*A(shape,z-s,-l-s,label_i)*shape.Rho(z-s))*(M_PI*(-1)*A(shape,z-s,-l-s,label_j)*shape.Rho(z-s));
    result = result*shape.get_density()*shape.get_average_v()/2.;
    return result;
}
double gsl_mini(shape shape, double fn1(double, void*), double xguess, double x_low, double x_high)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;
    F.function = fn1;
    F.params = &shape;

    T = gsl_min_fminimizer_quad_golden;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, xguess, x_low, x_high);
    do
    {
        iter++;
        xguess = gsl_min_fminimizer_x_minimum (s);
        x_low = gsl_min_fminimizer_x_lower (s);
        x_high = gsl_min_fminimizer_x_upper (s);
        status = gsl_min_test_interval (x_low, x_high, 1e-6, 0.0);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_min_fminimizer_free (s);
    return xguess;
}
double dissipative(shape shape, const char label_i, const char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double r = shape.get_r();
    double c = shape.get_c();
    double gamma_wall, gamma_ww;
    double result;
    double Rmin;
    double ks=0.27;
    gamma_wall = dissipative_wall(shape, label_i, label_j);
    gamma_ww = dissipative_wall2(shape, label_i, label_j)+dissipative_window(shape, label_i, label_j);
    if (c>0.) {
        Rmin = max(shape.Rho(gsl_mini(shape, RhoShape, z-s+1e-6, z-s, l-s)),
                   shape.Rho(gsl_mini(shape, RhoShape, z-s-1e-6, -l-s, z-s)))/pow(shape.get_Rcn(),2);
        result = ks*(pow(sin(M_PI*pow(r,2)/Rmin/2.),2)*gamma_wall +
                 pow(cos(M_PI*pow(r,2)/Rmin/2.),2)*gamma_ww);
    }
    else{
        result = ks*gamma_wall;
    }
    return result;
}