/**
# Breaking Stokes wave

We solve the two-phase Navier--Stokes equations using a
momentum-conserving transport of each phase. Gravity is taken into
account using the "reduced gravity approach". */

#include <sys/stat.h>
#include "utils.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "vtk.h"
#include "output_vtu_foreach.h"


/**
The primary parameters are the wave steepness $ak$ and the Reynolds
number. */
double ak = 0.35;
double RE = 40000.;
int LEVEL = 5;
int maxlevel = 9;
double Tend = 20;
double lx = 20;
double Ly = 1.;
double water_depth = 0.6;
/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
The density and viscosity ratios are those of air and water. */
#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. *T0* is the wave period. */

#define k_  (2.*pi)
#define h_   0.6
#define g_   1.
#define T0  (k_/sqrt(g_*k_))


/**
The program takes optional arguments which are the level of
refinement, steepness and Reynolds numbers. */
void move_origin(){
  origin (-lx/2., -water_depth);
}
void mask_domain(){
  mask(y > Ly-water_depth ? top : none);
}

int main (int argc, char * argv[])
{
  //Make folders for saving data
  mkdir("./vtk",0755);
  mkdir("./vtu",0755);
  mkdir("./surface_profiles", 0755);

  //Read ak and LEVEL from command line args.
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);

  L0 = lx;
  move_origin();

  //periodic (right);
  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  G.y = -g_;

  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
  N = 1 << LEVEL;
  DT = 1e-2;
  run();
}



/**
## Initial conditions
We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */

#include "test/stokes.h"

event init (i = 0)
{

  if (!restore ("restart")) {

    /**
    We need to make sure that fields are properly initialised before
    refinement below, otherwise a
    [-catch](/src/README#tracking-floating-point-exceptions) exception
    will be triggered when debugging. */

    event ("properties");
    mask_domain();
    do {
      fraction (f, wave(x,y));
      foreach()
	    foreach_dimension()
	    u.x[] = u_x(x,y) * f[];
    }

    /**
    We repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

    while (adapt_wavelet({f,u},
			  (double[]){0.01,uemax,uemax,uemax}, maxlevel, LEVEL).nf);
  }
}


//save the surface profile
event surface_profile (t += Tend/4.; t <= Tend) {
  char name[80];
  FILE * fp ;
  sprintf(name, "surface_profiles/res_%4.4f.txt",t);
  fp = fopen(name, "w"); 
  output_facets (f, fp);
  fprintf (fp, "\n");
  fclose (fp);
}


//save ordered mesh
event vtk(t+=0.2;t<Tend){
  char filename[40];
  sprintf(filename, "vtk/TIME-%06g.vtk", (t*100));
  FILE * fp = fopen(filename, "w");
  output_vtk(all, 100 , fp, true);
}

//save unordered mesh
event vtu(t+=0.2;t<Tend){
  char filename[40];
  sprintf(filename, "vtu/TIME-%06g", (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}


/**
## Mesh adaptation. We adapt the mesh according to the error on volume fraction
and velocity. */


event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){0.01,uemax,uemax,uemax}, maxlevel, LEVEL);
}




/*
event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += norm2*f[]*dv();
    gpe += y*f[]*dv();
  }
  printf ("%g %g %g\n", t/(k_/sqrt(g_*k_)), rho1*ke/2., rho1*g_*gpe + 0.125);
}
*/



event progressbar(i+=10)
{
  float progress = 0;
  int barWidth = 80;
  char filled[] = "===========================================================================================";
  char empy[] =   "                                                                                           ";

  progress = t /Tend;
  int filled_segments = barWidth*progress;
  int empty_segments = barWidth - filled_segments;

  printf("\r t=%.3f ", t);
  printf("[%.*s>", filled_segments, filled);
  printf("%.*s]\r", empty_segments, empy);
}


