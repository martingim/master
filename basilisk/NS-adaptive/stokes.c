/**
# intialised Stokes wave

We solve the two-phase Navier--Stokes equations using a
momentum-conserving transport of each phase. Gravity is taken into
account using the "reduced gravity approach". */
#include "adapt_wavelet_leave_interface.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "output_vtu_foreach.h"
#include "heights.h"
int vtucount = 0;

/**
The primary parameters are the wave steepness $ak$ and the wavenumber
*/
double ak = 0.17;
int LEVEL = 5;
int max_LEVEL = 8;
int padding =  2;
double Lx = 0.7903377744879982;

vector h[]; //scalar field of the distance from the surface, using heights.h
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files

int n_probes = 1;
double probe_positions[] = {0.00001};
/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.05;

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. *T0* is simulation duration */

#define k_  7.95
#define h_   0.6
#define g_   9.81
#define T0  10

/**
The program takes optional arguments which are the level of
refinement, steepness and Reynolds numbers. */

int main (int argc, char * argv[])
{
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) 
    {                 
      max_LEVEL = atoi(argv[j + 1]);
    }
  }
  sprintf(results_folder, "results/LEVEL%d", LEVEL);
  sprintf(vtu_folder, "%s/vtu", results_folder);


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
  

  L0 = Lx;   
  origin (-L0/2, -h_);
  periodic (right);

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 997;
  rho2 = 1.204;
  mu1 = 8.9e-4;
  mu2 = 17.4e-6;
  G.y = -g_;

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
    
  fraction (f, wave(x, y)); //set the water depth
  while (adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax,uemax,uemax, uemax},max_LEVEL, LEVEL,padding).nf){  //for adapting more around the piston interface
    fraction (f, wave(x, y)); //set the water level on the refined mesh
  }
  foreach()
    foreach_dimension()
      u.x[] = u_x(x,y) * f[];
  

  /**
  On trees, we repeat this initialisation until mesh adaptation does
  not refine the mesh anymore. */


}

// event profiles (t += T0/4.; t <= 2.5*T0) {
//   output_facets (f, stderr);
//   fprintf (stderr, "\n");
// }

// event logfile (i++)
// {
//   double ke = 0., gpe = 0.;
//   foreach (reduction(+:ke) reduction(+:gpe)) {
//     double norm2 = 0.;
//     foreach_dimension()
//       norm2 += sq(u.x[]);
//     ke += norm2*f[]*dv();
//     gpe += y*f[]*dv();
//   }
//   printf ("%g %g %g\n", t/(k_/sqrt(g_*k_)), rho1*ke/2., rho1*g_*gpe + 0.125);
// }

//save unordered mesh
event vtu(t+=.1, last){
  printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}

event adapt (i++) {
  adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax, uemax, uemax, uemax}, max_LEVEL, LEVEL,padding);
}


event surface_probes(t+=0.01, t<T0){
  char filename[100];
  sprintf(filename, "%s/surface_probes.csv", results_folder);
  FILE *fp = fopen(filename, "a");
    if (i==0){
      fprintf(fp, "time");
      for (int j = 0;j<n_probes;j++){
        fprintf(fp, ", %f", probe_positions[j]);    
      }
      fprintf(fp, "\n");
  }
  fprintf(fp, "%f", t);
  heights(f, h);
  double min_height; //vertical distance from the center of the cell to the air water interface
  double surface_elevation=0;
  for (int probe_n=0;probe_n<n_probes;probe_n++){
    min_height=1;
    foreach(serial){
      if((fabs(x-probe_positions[probe_n])<Delta/2.0+1e-9)&&(fabs(y)<0.1)){
        if(h.y[] != nodata){
          if (fabs(min_height)>fabs(height(h.y[])*Delta)){
            min_height=height(h.y[])*Delta;
            surface_elevation = y + height(h.y[])*Delta;
          }
        }
      }
    }
    fprintf(fp, ", %f", surface_elevation);
  }
  fprintf(fp, "\n");
  fclose(fp);
}
