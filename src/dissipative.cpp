#include <string>
#include <gsl/gsl_integration.h>
#include "../include/shape.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

double dissipative_wall(shape shape, char label_i, char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double x, w;
    size_t n = 28;
    double result = 0;
    gsl_integration_glfixed_table *table = nullptr;
    table = gsl_integration_glfixed_table_alloc(n);
    for(size_t i=0; i<n; i++)
    {
        gsl_integration_glfixed_point(-1., 1., i, &x, &w, table);
        x = (l + z)/2. + s + (-l + z)/2.;
        result += ((shape.RhoDerivative(x,label_i)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('L',label_i))*
                   (shape.RhoDerivative(x,label_j)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('L',label_j))*
                   sqrt(4.*shape.Rho(x)+pow(shape.RhoDerivative(x,'x'),2)))*w;
        x = (l - z)/2. - s + (l + z)/2.;
        result += ((shape.RhoDerivative(x,label_i)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('L',label_i))*
                   (shape.RhoDerivative(x,label_j)+shape.RhoDerivative(x,'x')*shape.CenterOfMassDerivative('L',label_j))*
                   sqrt(4*shape.Rho(x)+pow(shape.RhoDerivative(x,'x'),2)))*w;
    }
    gsl_integration_glfixed_table_free(table);
    result = result*M_PI*shape.get_average_v()/2.;
    return result;
}
double dissipative_window(shape shape, const char label_i, const char label_j){
    double l = shape.get_l();
    double z = shape.get_z();
    double s = shape.get_s();
    double sigma = M_PI*pow(shape.get_r(),2);
    double result;
    result = sigma*((shape.CenterOfMassDerivative('R',label_i)-shape.CenterOfMassDerivative('R',label_i))*
            (shape.CenterOfMassDerivative('R',label_j)-shape.CenterOfMassDerivative('L',label_j)))+
            32./(9.*sigma)*(M_PI*A_para(shape, z-s,-l-s,label_i)*shape.Rho(-l-s))*(M_PI*A_para(shape, z-s,-l-s,label_j)*shape.Rho(-l-s));
    result = result*shape.get_density()*shape.get_average_v()/2.;
    return result;
}

double gsl_mini(double fn1(double, void*), double xguess, double x_low, double x_high)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;
    F.function = fn1;
    T = gsl_min_fminimizer_quad_golden;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, xguess, x_low, x_high);
    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);
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
    double result;
    double Rmin;
    double ks=0.27;
    Rmin = min(gsl_mini(RhoShape,(-l+z-2*s)/2., -l-s, z-s),gsl_mini(RhoShape, (l-s+z-s)/2., z-s, l-s));
    result = pow(sin(M_PI*pow(r/Rmin,2)/2),2)*dissipative_wall(shape,label_i,label_j)+
            ks*pow(cos(M_PI*pow(r/Rmin,2)/2),2)*dissipative_window(shape,label_i,label_j);
    return result;
}



























