/**
# intialised third order stokes wave with periodic boundary conditions

based on bar.c
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

double ak = 0.16;
double Tend = 5;
double Lx = 0.7903377744879982;
double k = 7.95;
double g = 9.81;
int LEVEL = 7;
double rmin = 0.5;

#define nl_ 10
#define k_ k
#define h_ 0.6
#define g_ g
#define T0 0.7115297011824904
#include "test/stokes.h"


/**
We use Stokes third order wave from src/test/stokes.h to initialise the surface elevation 
and the velocity field*/

event init (i = 0)
{
  geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
  //set the surface of the wave
  foreach() {
    zb[] = -h_;
    foreach_layer(){
      h[] = (max(- zb[], 0.)+wave(x, y))*beta[point.l];
    }
  }
  /*set the velocity field of the wave*/
  foreach(){
    double z = zb[];
    foreach_layer(){
      z += h[]/2;
      u.x[] = u_x(x,z);
      w[] = u_y(x,z);
      z += h[]/2;
    }
  }
}

int main(int argc, char *argv[])
  {
  if (argc>1)
    LEVEL = atoi(argv[1]);
    if (argc>2)
      nl = atoi(argv[2]);
    else
      nl = nl_;
  
  origin (-Lx/2., h_);
  periodic(right);
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
  sprintf(filename, "/home/martin/Documents/master/matlab/PIV_basilisk/basilisk_results/velocities_nx%d_nl%d_timestep_%d.csv",1<<LEVEL, nl, i);
  fprintf(stderr, "saving results to:%s\n", filename);
  FILE *fp = fopen(filename, "w"); //if at the first timestep overwrite the previous file, can later add run parameters here
  fprintf(fp, "\"Time\",\"layer\",\"Points:0\",\"Points:1\",\"Points:2\",\"u.x\",\"u.z\",\"eta\"\n");
  foreach(){
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
  sprintf(filename, "/home/martin/Documents/master/matlab/PIV_basilisk/basilisk_results/energy_nx%d_nl%d.csv",1<<LEVEL, nl);
  static FILE * fp = fopen(filename, "w");
  double ke = 0.;
  double gpe = 0.;
  foreach(){
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
	   "p [-%f:%f][-0.6:0.04]'-' u 1:3:2 w filledcu lc 3 t ''\n", t, Lx/2,Lx/2);
  foreach()
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

// gauges to compare the surface elevation
// Gauge gauges[] = {
//   {"X_0",  0},
//   {NULL}
// };


// event output (i += 2; t <= Tend)
//   output_gauges (gauges, {eta});
