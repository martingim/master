/**
# Breaking Stokes wave

We solve the two-phase Navier--Stokes equations using a
momentum-conserving transport of each phase. Gravity is taken into
account using the "reduced gravity approach". */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "output_vtu_foreach.h"
#define TREE 0
/**
The wave number, fluid depth and acceleration of gravity are set to
these values. */
  
double k_ = 7.9583, h_ = 0.6, g_ = 9.81;

/**
The primary parameters are the wave steepness $ak$ and the Reynolds
number. */

double ak = 0.16;
double RE = 40000.;
int LEVEL = 9;
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files
char energy_file[60];
/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
The density and viscosity ratios are those of air and water. */

#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
*T0* is the wave period. */

#define T0  (k_*L0/sqrt(g_*k_))

/**
The program takes optional arguments which are the level of
refinement, steepness and Reynolds numbers. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);

  /**
  The domain is a cubic box centered on the origin and of length
  $L0=1$, periodic in the x-direction. */
   
  origin (-L0/2, -L0/2);
  periodic (right);

  sprintf(results_folder, "results/LEVEL%d", LEVEL);
  sprintf(vtu_folder, "%s/vtu", results_folder);
  sprintf(energy_file, "%s/energy.txt", results_folder);


  char remove_old_results[100];
  sprintf(remove_old_results, "rm -r %s", results_folder);
  if (system(remove_old_results)==0){
    printf("removed old results in:%s\n", results_folder);
  }

  char make_results_folder[100];
  sprintf(make_results_folder, "mkdir -p %s", vtu_folder);
  if (system(make_results_folder)==0){
    printf("made results folder:%s\n", results_folder);
  }
  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 997.;
  rho2 = 1.204;
  mu1 = 8.9e-4; 
  mu2 = 17.4e-6;
  G.y = -g_;
  L0 = 0.7895;
  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
#if TREE  
  N = 32;
#else
  N = 1 << LEVEL;
#endif
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
    
    do {
      fraction (f, wave(x,y));
      foreach()
	foreach_dimension()
	  u.x[] = u_x(x,y) * f[];
    }

    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

#if TREE
    while (adapt_wavelet ({f,u},
			  (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5).nf);
#else
    while (0); // to match 'do' above
#endif
  }
}

///event profiles (t += T0/4.; t <= 2.5*T0) {
  //output_facets (f, stderr);
  //fprintf (stderr, "\n");
//}

event logfile (i++)
{
  FILE *fp = fopen(energy_file, "a");
  if (i==0){
    fprintf(fp, "t, ke, gpe\n");
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += norm2*f[]*dv();
    gpe += y*f[]*dv();
  }
  fprintf (fp,"%f %.5f %.5f\n", t, rho1*ke/2., rho1*g_*gpe);
  fclose(fp);
}
//save unordered mesh
event vtu(t+=.1, t<10){
  //printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}
/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}
#endif

