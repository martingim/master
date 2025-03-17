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
#include "heights.h"
int vtucount = 0;
int set_n_threads = 2;

/**
The primary parameters are the wave steepness $ak$ and the Reynolds
number. */
#define TREE 0
double ak = 0.16314515;
int LEVEL = 7;
double Lx = 0.7895;//set to exactly one wavelength

vector h[]; //scalar field of the distance from the surface, using heights.h
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files
char energy_file[60];

int n_probes = 1;
double probe_positions[] = {0.00001};

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. *T0* is simulation duration */

#define k_  7.9583
#define h_   0.6
#define g_   9.81
#define T0  10

/**
The program takes an optional argument which is the level of
refinement */

int main (int argc, char * argv[])
{
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) 
    {                 
      LEVEL = atoi(argv[j + 1]);
    }
  }
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

  L0 = Lx;   
  origin (-L0/2, -h_);
  periodic (right);

  rho1 = 997;
  rho2 = 1.204;
  mu1 = 8.9e-4;
  mu2 = 17.4e-6;
  G.y = -g_;
  N = 1 << LEVEL;
  DT = 1e-2;
#if _OPENMP
  int num_omp = omp_get_max_threads();
  fprintf(stderr, "max number of openmp threads:%d\n", num_omp);
  if (set_n_threads){  //set number of omp threads
    omp_set_num_threads(set_n_threads);
  }
  fprintf(stderr, "set openmp threads:%d\n", omp_get_max_threads());
#endif
  run();
}

/**
Initialise the wave
using the third-order Stokes wave solution. 
*/
#include "test/stokes.h"

event init (i = 0)
{
  
  fraction (f, wave(x,y));
  foreach()
	  foreach_dimension()
	    u.x[] = u_x(x,y) * f[];
}

// event profiles (t += T0/4.; t <= 2.5*T0) {
//   output_facets (f, stderr);
//   fprintf (stderr, "\n");
// }

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

  fprintf(fp, "%g, %e, %e\n", t, rho1*ke/2., rho1*g_*gpe);
  fclose(fp);
}

//save unordered mesh
event vtu(t+=.1, last){
  //printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}

// void save_vtu ( int nf, int j)
// {
//   char name[80];
//   FILE * fp ;
//   sprintf(name, "RES_VTK/res_%4.4d.vtu",j);
//   //nf > 0 ? sprintf(name, "RES_VTK/res_n%3.3d_%4.4d.vtu",pid(),j) : sprintf(name, "RES_VTK/res_%4.4d.vtu",j);
//   fp = fopen(name, "w"); 
//   output_vtu_bin_foreach ((scalar *) {f,p}, (vector *) {u}, N, fp, false); 
//   fclose (fp);
// }


// event logfilevtu (t=0.0;t<=10;t+=0.04) 
// {
//   save_vtu(0,vtucount);
//   vtucount += 1;
// }


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
