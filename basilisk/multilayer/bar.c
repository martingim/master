/**
# Sinusoidal wave propagation over a bar

[Beji and Battjes, 1993](/src/references.bib#beji1993) and [Luth et
al, 1994](/src/references.bib#luth1994) studied experimentally the
transformation of sinusoidal waves propagating over a submerged bar
(or reef). This is a good test case for dispersive models as higher
harmonics are nonlinearly generated and released with phase shifts
corresponding to the dispersion relation. 

This test case is discussed in [Popinet
(2020)](/Bibliography#popinet2020) for the layered version. */

#include "grid/multigrid1D.h"

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

#include "layered/check_eta.h"
#include "layered/perfs.h"


#define ML 1

double ak = 0.16;
double Tend = 100;
double Lx = 1.5806755489759965;
double Ly = 1;
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
  //geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
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


void move_origin(){
  if (Lx>Ly){  
  size(Lx);
  }
  else{
  size(Ly);
  }
  origin (-Lx/2., h_);
}

int main(int argc, char *argv[])
  {
  if (argc>1)
    LEVEL = atoi(argv[1]);
    if (argc>2)
      nl = atoi(argv[2]);
    else
      nl = nl_;
  move_origin();
  periodic(right);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  breaking = 0.1;
  CFL_H = 0.5;
  run();
}

event save_velocity(t += Tend/10.; t<=Tend)
{
  char filename[200];
  sprintf(filename, "/home/martin/Documents/master/matlab/ex2/basilisk_results/velocities_nx%d_nl%d_timestep_%d.csv",1<<LEVEL, nl, i);
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

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=40$. 

![Snapshot of waves. The top of the bar is seen in white.](bar/snapshot.png)
*/

void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [-0.8:0.8][-0.6:0.04]'-' u 1:3:2 w filledcu lc 3 t ''\n", t);
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