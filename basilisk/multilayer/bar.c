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


#include "vtk.h"
#define ML 1

double ak = 0.16;
double Tend = 10;
double Lx = 1.5806755489759965;
double Ly = 1;
double k = 7.95;
double g = 9.81;

int LEVEL = 6;
#define nl_ 10
#define rmin 0.5
#define k_ k
#define h_ 0.6
#define g_ g
#define T0 0.7115297011824904
#include "test/stokes.h"

void move_origin(){
  if (Lx>Ly){  
  size(Lx);
  }
  else{
  size(Ly);
  }
  origin (-Lx/2., h_);
}




int main() {
  move_origin();
  periodic(right);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  nl = nl_;  
  breaking = 0.1;
  CFL_H = 0.5;
  run();
}

/**
We use ["radiation"
conditions](/src/elevation.h#radiation-boundary-conditions) at the
inlet and outlet. At the inlet (on the left), we try to impose the
desired sinusoidal wave form. We have to tune the amplitude to obtain
the required amplitude as measured in the experiment at gauge 4. The
period of 2.02 seconds matches that of the experiment. */

event init (i = 0)
{
  geometric_beta(rmin, true);
  /**
  set the surface of the wave */

  foreach() {
    zb[] = -h_;
    foreach_layer()
      h[] = (max(- zb[], 0.)+wave(x, y))*beta[point.l];
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
  //fprintf (stderr, "%g %f %g %g\n", t, interpolate (eta, 17.3, 0.), ke, gpe);
  
  #if 0
  //print the height of the top of the layers for x = 1
  double layer_z = 0;
  foreach_layer() 
  {
    layer_z += interpolate(h, 0);
    fprintf (stderr, "%f, %d\n", layer_z, _layer);
  }
  #endif
}

/**
This optionally displays consistency between `res_eta` and `deta`
(corresponding to the `check_eta.h` option above). */


event gnuplot (t = Tend) {
  FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set term pngcairo enhanced size 640,640 font \",8\"\n"
           "set output 'snapshot.png'\n");
  plot_profile (t, fp);
}

event save_velocity(t += 1; t<=Tend)
{
  char filename[200];
  sprintf(filename, "/home/martin/Documents/master/matlab/ex2/basilisk_results/velocities_nx%d_nl%d.csv",1<<LEVEL, nl_);
  FILE *fp;
  if (t==0)
    fp = fopen(filename, "w"); //if at the first timestep overwrite the previous file, can later add run parameters here
  else
    fp = fopen(filename, "a");

  foreach(){
    double z = zb[];
    foreach_layer(){
      z += h[]/2;
      fprintf(fp, "%f, %d, %f, %f, %f, %f, %f,\n", t, point.l, x, z, u.x[], w[], eta[]);
      z += h[]/2;
    }
  }
  fclose(fp); 
}


// gauges to compare the surface elevation
// Gauge gauges[] = {
//   {"X_0",  0},
//   {NULL}
// };


// event output (i += 2; t <= Tend)
//   output_gauges (gauges, {eta});
