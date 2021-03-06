#include <iostream>
#include "../include/shape.h"
#include <boost/format.hpp>
#include <fstream>
#include <omp.h>
#include <mpi.h>

using namespace std;
using namespace boost;

void input(const string&, int*, double*, double*, int*);
void store(const string&, int*, double*);
double dissipative(shape shape, const char label_i, const char label_j);
double inertia(shape shape, const char label_i, const char label_j);
double inertia_df(shape shape, const char label_i, const char label_j, const char label_l);
void runge_kutta(gsl_vector *generalized_coordinates, gsl_vector *generalized_momenta,
                 double starting[5], double step_length[5], double storation[],
                 shape shape, int rank);
// 命令行参数是激发能
int main(int argc, char* argv[])
{
    int rank, size, ierr;
    MPI_Status status;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    shape shape;
    int steps[5],gs[5],sp[5];
    double starting[5],step_length[5];
    input("U236_Input.dat",steps,starting,step_length,gs);
    gs[0]=1;gs[1]=1;gs[2]=2;gs[3]=0;gs[4]=2;
    starting[1] = 1 / sqrt(starting[0]+gs[0]*step_length[0]) + 0.05;
    starting[3] = -pow(starting[0]+gs[0]*step_length[0], -2.5) - 0.5;
    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
    shape.set_gs(gs);
    shape.set_steps(steps);
    shape.set_level_density();
    // 初始化storation为20.
    fill(storation,storation+steps[0]*steps[1]*steps[2]*steps[3]*steps[4],10.);
    store("U236.txt",steps,storation);
    double ground_energy = shape.grid_energy(storation, gs);
    shape.set_ground_energy(ground_energy);
    shape.set_excited_energy(stof(argv[1]));
    // storation 中间存储的就是势能曲面的值
    // inertia_tensor是惯性张量，dissipative_tenson是耗散张量
    gsl_vector *generalized_coordinates, *generalized_momenta;
    generalized_coordinates = gsl_vector_alloc(5);
    generalized_momenta = gsl_vector_alloc(5);

    FILE *myfile;
    string Langevin_Result;
    Langevin_Result = "Langevin_Result_"+to_string(rank)+".txt";
    myfile = fopen(Langevin_Result.c_str(),"w");
//    myfile.close();

    //这里gamma_tensor的结果就是dissipative_tensor的平方根
    int i;
    omp_set_num_threads(2);
    #pragma omp parallel for default(none) schedule(dynamic) firstprivate(myfile, shape, generalized_coordinates, generalized_momenta) private(i) shared(rank, size, cout, sp, starting, step_length, storation)
    for(i=0+rank; i<10000; i=i+size) {
        sp[0]=10;sp[1]=9;sp[2]=16;sp[3]=11;sp[4]=7;
        starting[1] = 1 / sqrt(starting[0]+sp[0]*step_length[0]) + 0.05;
        starting[3] = -pow(starting[0]+sp[0]*step_length[0], -2.5) - 0.5;

        for (size_t j=0; j<5; j++){
            gsl_vector_set(generalized_momenta, j, 0);
            gsl_vector_set(generalized_coordinates, j, sp[j]*step_length[j]+starting[j]);
        }
//        #pragma omp critical
//        if (rank==0) {
//            fprintf(myfile,"%d\n",omp_get_thread_num()+100);
//        }
//        else{
//            fprintf(myfile,"%d\n",omp_get_thread_num());
//        }
        runge_kutta(generalized_coordinates, generalized_momenta,
                    starting, step_length, storation,
                    shape, rank);
        // cout<<format("generalized coordinates: %1$+7.3f %2$+7.3f %3$+7.3f %4$+7.3f %5$+7.3f\n")
        //       %gsl_vector_get(generalized_coordinates,0)
        //       %gsl_vector_get(generalized_coordinates,1)
        //       %gsl_vector_get(generalized_coordinates,2)
        //       %gsl_vector_get(generalized_coordinates,3)
        //       %gsl_vector_get(generalized_coordinates,4);
//        #pragma omp critical
        fprintf(myfile,"%d %7.3f %7.3f %7.3f %7.3f %7.3f %10.3f\n",
                i,
                gsl_vector_get(generalized_coordinates,0),
                gsl_vector_get(generalized_coordinates,1),
                gsl_vector_get(generalized_coordinates,2),
                gsl_vector_get(generalized_coordinates,3),
                gsl_vector_get(generalized_coordinates,4),
                shape.AH(generalized_coordinates));
    }
    fclose(myfile);
    delete(storation);
    gsl_vector_free(generalized_coordinates);
    gsl_vector_free(generalized_momenta);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Finalize();
    return 0;
}
