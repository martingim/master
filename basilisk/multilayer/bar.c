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
int LEVEL = 6;
double Tend = 10;
double Lx = 1.5806755489759965;
double Ly = 1;
double k = 7.95;
double g = 9.81;


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
#if ML
  nl = 20;  
  breaking = 0.1;
  CFL_H = 0.5;
#endif
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
  
  /**
  set the surface of the wave */

  foreach() {
    zb[] = -h_;
    foreach_layer()
      h[] = (max(- zb[], 0.)+wave(x, y))/nl;
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


void velocity_profile()
  {
  char filename[40];
  sprintf(filename, "velocities.csv");
  FILE *fp = fopen(filename, "a");
  foreach(){
    double z = zb[];
    foreach_layer(){
      z += h[]/2;
      fprintf(fp, "%f, %d, %f, %f, %f, %f, %f,\n", t, point.l, x, z, u.x[], w[], eta[]);
      z += h[]/2;
    }
  }
  /*
  foreach(){
    double z_layer = zb[];
    double prev_h = zb[]*0;
  }
    
  foreach_layer() 
  {
    z_layer += prev_h/2;
    z_layer += h[]/2;
    prev_h = h[];
    fprintf (fp, "%f, %f, %f, %d", t, interpolate(z_layer,x), x, _layer);
    fprintf(fp, "\n");
  }
  */
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
#if ML      
      double norm2 = sq(w[]);
#else
      double norm2 = 0.;
#endif
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

/**
The location of the gauges is difficult to find in the litterature, we
used a combination of [Yamazaki et al,
2009](/src/references.bib#yamazaki2009) and [Dingemans,
1994](/src/references.bib#dingemans1994). */

Gauge gauges[] = {
  {"X_0",  0},
  {NULL}
};


event output (i += 2; t <= 40)
  output_gauges (gauges, {eta});


event save_velocity(t += 1; t<=Tend)
  velocity_profile(); 
//event velocity_profile(i += 2; t <= Tend)
  //double layer_z = 0;
  //foreach_layer() 
  //{
  //  layer_z += interpolate(h, 0.0);
  //  fprintf (stderr, "%f, %d\n", layer_z, _layer);
  //}
  //double X_0 = 0;
  //char filename[40];
  //sprintf(filename, "velocity_profile_at_x_%f.txt", X_0);
  //FILE *fp = fopen(filename, "w");
  //double z = 0;
  //foreach_layer()
    //z += interpolate(h, X_0);
    //fprintf(fp, "%f, %f, %f, %f, %f\n", t, X_0, z, u.x, u.y);
  
//save ordered mesh
event save_vtk(t+=0.2;t<10){
  char filename[40];
  sprintf(filename, "vtk/TIME-%06g.vtk", (t*100));
  FILE *fp = fopen(filename, "w");
  output_vtk((scalar *) {eta}, 100, fp, true);
}
// //save unordered mesh
// event vtu(t+=0.2;t<10){
//   char filename[40];
//   sprintf(filename, "vtu/TIME-%06g", (t*100));
//   output_vtu((scalar *) {eta}, (vector *) {u}, filename);
// }


/**
The modelled and experimental (circles) timeseries compare quite
well. The agreement is significantly better than that in [Yamazaki et
al, 2009](/src/references.bib#yamazaki2009) (figure 4) in particular
for gauge 9, but probably not as good as that in [Lannes and Marche,
2014](/src/references.bib#lannes2014) (figure 12), who used a
higher-order scheme, and a three-parameter optimised dispersion
relation. Note that using the optimised dispersion relation (with
$\alpha_d=1.153$) is necessary to obtain such an agreement.

~~~gnuplot Comparison of experimental and numerical timeseries
set term svg enhanced size 640,480 font ",10"
set xrange [33:39]
set yrange [-2:4]
set ytics -2,2,4
set key top left
set multiplot layout 4,2 scale 1.05,1.1
set rmargin 2.
set tmargin 0.5
unset xtics
# t0 is a tunable parameter
t0 = -0.24
plot 'WG4' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 4', \
     '../gauge-4' pt 6 lc -1 t ''
unset ytics
plot 'WG5' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 5', \
     '../gauge-5' pt 6 lc -1 t ''
set ytics -2,2,6
plot 'WG6' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 6', \
     '../gauge-6' pt 6 lc -1 t ''
unset ytics
plot 'WG7' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 7', \
     '../gauge-7' pt 6 lc -1 t ''
set ytics -2,2,6
plot 'WG8' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 8', \
     '../gauge-8' pt 6 lc -1 t ''
unset ytics
plot 'WG9' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 9', \
     '../gauge-9' pt 6 lc -1 t ''
set xtics
set ytics -2,2,6
plot 'WG10' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 10', \
     '../gauge-10' pt 6 lc -1 t ''
unset ytics
plot 'WG11' u ($1+t0):($2*100.) w l lc -1 lw 2 t 'gauge 11', \
     '../gauge-11' pt 6 lc -1 t ''
unset multiplot
~~~
*/
