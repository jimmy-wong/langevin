#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>

void gauss_rand(double* u, unsigned long* rng_seed)
{
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_taus;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, *rng_seed);
    *rng_seed = gsl_rng_get(r);
    *u = gsl_ran_gaussian(r, 1.);
    gsl_rng_free (r);
    std::ofstream myfile;
    myfile.open ("example.txt",std::ios_base::app);
    myfile<<*u<<"\n";
    myfile.close();
}
