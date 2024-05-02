/**
# intialised third order stokes wave with periodic boundary conditions

based on bar.c and breaking.c
the parameters for the wave are from 
  Jensen, A., Sveen, J. K., Grue, J., Richon, J. B., & Gray, C. (2001). 
  Accelerations in water waves by extended particle image velocimetry. 
  Experiments in Fluids, 30(5), 500–510. https://doi.org/10.1007/s003480000229
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

#include "layered/check_eta.h"
#include "layered/perfs.h"

double Tend = 30;    //the end time of the simulation
double Lx = 24.6; //The length of the simulation domain
double g = 9.81;
double omega = 8.95;
double piston_amplitude = 0.0128;
int LEVEL = 9;      //the grid resolution in x direction Nx = 2**LEVEL
double rmin = 0.5;  //rmin the relative height of the top layer compared to 
                    //a regular distribution. the rest fo the layers follow a geometric distribution.
double h_ = 0.6;                   
char save_location[] = "./";

#define nl_ 5  //the default number of layers if none are given as command line arguments
#define g_ g

//U=d/dt (piston_amplitude*tanh(t)*sin(omega*t))
#define U (piston_amplitude*(sin(omega*t)/cosh(t)/cosh(t)+tanh(t)*omega*cos(omega*t)))

/**
We use Stokes third order wave from src/test/stokes.h to initialise the surface elevation 
and the velocity field*/

event init (i = 0)
{
  u.n[left]  = radiation (U);
  u.n[right] = + radiation (0);

  //geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
  foreach() {
    zb[] = -h_;
    foreach_layer(){
      h[] = h_/nl;
    }
  }
#if _MPI
  fprintf(stderr, "mpi\n");
#elif _OPENMP
  int num_omp = omp_get_max_threads();
  fprintf(stderr, "number of openmp threads:%d\n", num_omp);
#endif

}

int main(int argc, char *argv[])
  {
  if (argc>1)
    LEVEL = atoi(argv[1]);
    if (argc>2)
      nl = atoi(argv[2]);
    else
      nl = nl_;
  
  //origin (0, h_);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  breaking = 0.1;
  CFL_H = .5;
  TOLERANCE = 10e-5;
  //theta_H = 0.51;
  run();
}

event save_velocity(t += Tend; t<=Tend)
{
  char filename[200];
  sprintf(filename, "%svelocities_nx%d_nl%d_timestep_%d.csv", save_location, 1<<LEVEL, nl, i);
  fprintf(stderr, "saving results to:%s\n", filename);
  FILE *fp = fopen(filename, "w"); //if at the first timestep overwrite the previous file, can later add run parameters here
  fprintf(fp, "\"Time\",\"layer\",\"Points:0\",\"Points:1\",\"Points:2\",\"u.x\",\"u.z\",\"eta\"\n");
  foreach(serial){
    double z = zb[];
    foreach_layer(){
      z += h[]/2;
      fprintf(fp, "%f, %d, %f, %f, %f, %f, %f, %f\n", t, point.l, x, 0., z, u.x[], w[], eta[]);
      z += h[]/2;
    }
  }
  fclose(fp); 
}

event save_energy(i++)
{
  char filename[200];
  sprintf(filename, "%senergy_nx%d_nl%d.csv", save_location, 1<<LEVEL, nl);
  static FILE * fp = fopen(filename, "w");
  double gpe = 0.;
  double ke = 0.;

  foreach(reduction(+:ke) reduction(+:gpe)){
    double z = zb[]*0;
    foreach_layer(){
      double norm2 = sq(w[]);
      foreach_dimension()
	      norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (z + h[]/2)*h[]*dv();
      z += h[];
    }
  }
  if (i == 0)
    fprintf (fp, "ke, gpe, t\n");
  fprintf(fp, "%f, %f, %f\n", ke, gpe, t);
}

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=Tend$.*/

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [0.0:%f][-0.8:0.6]'-' u 1:3:2 w filledcu lc 3 t ''\n", t, Lx);
  foreach(serial)
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}


event profiles (t += 0.05)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	      norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n");
  plot_profile (t, fp);
}

event gnuplot (t = Tend) {
  FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set term pngcairo enhanced size 640,640 font \",8\"\n"
           "set output 'snapshot.png'\n");
  plot_profile (t, fp);
}

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%g, ", t, i, dt);
  printf("%.2f%%\n", progress*100);
}


//gauges to compare the surface elevation
Gauge gauges[] = {
  {"X_0",  8.009},
  {"X_1",  10.048},
  {"X_2",  10.745},
  {"X_3",  11.498},
  {NULL}
  };


event output (t += 0.01; t <= Tend)
  output_gauges (gauges, {eta});
