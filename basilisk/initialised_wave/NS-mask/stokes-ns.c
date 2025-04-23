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
double ak = 0.16;
double RE = 45491.;
int LEVEL = 5;
int maxlevel = 8;
double Tend = 5;
double Ly = 1.; //the height of the waveflume
double water_depth = 0.6;
double k = 7.95;
double lx = 0.7903377744879982*2; //must be larger than Ly to get the periodic boundary right

/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/*the error in the VOF field*/
double femax = 0.0001;
/**
The density and viscosity ratios are those of air and water. */
#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. *T0* is the wave period. */

double k_  = 7.9583;
double h_  = 0.6;
double g_  =  9.81;
#define T0  0.7115297011824904


/**
The program takes optional arguments which are the level of
refinement, steepness and Reynolds numbers. */
void move_origin(){
  if (lx>Ly){  
  size(lx);
  }
  else{
  size(Ly);
  }
  origin (-lx/2., -water_depth);
}

void mask_domain(){
  if (lx>Ly){  
    mask(y > Ly-water_depth ? top : none);
  }
  else{
    mask(x >= lx/2. ? right : none);
  }
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
  periodic (right);
  
  
  rho1 = 1025;
  rho2 = 1.225;
  mu1 = 8.9e-4; 
  mu2 = 17.4e-6;
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
			  (double[]){femax,uemax,uemax,uemax}, maxlevel, LEVEL).nf);
  }
  
}

/**
## Mesh adaptation. We adapt the mesh according to the error on volume fraction
and velocity. */


event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){femax,uemax,uemax,uemax}, maxlevel, LEVEL);
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
event vtk(t+=1;t<Tend){
  char filename[40];
  sprintf(filename, "vtk/TIME-%06g.vtk", (t*100));
  FILE * fp = fopen(filename, "w");
  output_vtk(all, 100 , fp, true);
}

//save unordered mesh
event vtu(t+=1;t<Tend){
  char filename[40];
  sprintf(filename, "vtu/TIME-%06g", (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}






event logfile (i++)
{
  char filename[200];
  sprintf(filename, "/home/martin/Documents/master/matlab/PIV_basilisk/basilisk_results/energy_ns_coarseN%d_fineN%d.csv",1<<LEVEL, 1<<maxlevel);
  static FILE * fp = fopen(filename, "w");
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += norm2*f[]*dv();
    gpe += y*f[]*dv();
  }
  if (i == 0)
    fprintf (fp, "ke, gpe, t\n");
  fprintf(fp, "%f, %f, %f\n", rho1*ke/2., rho1*g_*gpe, t);
}




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


