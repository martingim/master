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
#include "output_pvd.h"


int set_n_threads = 2; //set to 0 to use all available threads (I got best results with 1 thread per core)
double Tend = 10;     //the end time of the simulation
#define nl_ 10  //the default number of layers if none are given as command line arguments
double rmin = 0.3;  //rmin the relative height of the top layer compared to a regular distribution. the rest fo the layers follow a geometric distribution.
int LEVEL = 7;      //the grid resolution in x direction Nx = 2**LEVEL
double rho = 997;
double k_ = 7.9596;    //the wavenumber
double a = 0.0205;    //the amplitude of the wave
double h_ = 0.6;      //water depth

int n_waves = 1;      //number of waves to fit the domain to 
#define Lx 2*pi/k_*n_waves  // the length of the simulation domain
#define ak a*k_  //wave steepness
                    
char results_folder[40]; //the location to save the results
char vts_folder[50]; //the locaton to save the vtu files
char guage_name[50];

#define g_ 9.81
#include "test/stokes.h" //third order stokes wave


/**
We use Stokes third order wave from src/test/stokes.h to initialise the surface elevation 
and the velocity field*/

event init (i = 0)
{
  geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
  foreach() {
    zb[] = -h_;
    foreach_layer(){
      h[] = (max(- zb[], 0.)+wave(x, z))*beta[point.l]; //set the thickess of the layer based on the surface of the wave
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
  nl=nl_;
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) 
    {                 
      LEVEL = atoi(argv[j + 1]);
    }
    if (strcmp(argv[j], "-nl") == 0)
    {                
      nl = atoi(argv[j + 1]); 
    }
    if (strcmp(argv[j], "-n") == 0)
    {                
      n_waves = atoi(argv[j + 1]);
    }
  }

  sprintf(results_folder, "results/LEVEL%d_layers%d", LEVEL, nl);
  sprintf(vts_folder, "%s/vts", results_folder);
  
  char remove_old_results[100];
  sprintf(remove_old_results, "rm -r %s", results_folder);
  if (system(remove_old_results)==0){
    printf("removed old results in:%s\n", results_folder);
  }

  char make_results_folder[100];
  sprintf(make_results_folder, "mkdir -p %s", vts_folder);
  if (system(make_results_folder)==0){
    printf("made results folder:%s\n", results_folder);
  }
  sprintf(guage_name, "%s/X_0", results_folder);
  
  //copy the script to the results folder for later incpection if needed
  char copy_script[100];
  sprintf(copy_script, "cp multilayer.c %s/multilayer.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }



  origin (-Lx/2., h_);
  periodic(right);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g_;
  breaking = 0.1;
  CFL_H = .1;
  TOLERANCE = 10e-5;
  theta_H = 0.51;
  nu = 1e-6;
  printf("nu:%e\n", nu);
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

event save_velocity(t += Tend; t<=Tend)
{
  char filename[200];
  sprintf(filename, "%s/velocities_nx%d_nl%d_timestep_%d.csv", results_folder, 1<<LEVEL, nl, i);
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
  sprintf(filename, "%s/energy_nx%d_nl%d.csv", results_folder, 1<<LEVEL, nl);
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
    fprintf (fp, "t, ke, gpe\n");
  fprintf(fp, "%f, %e, %e\n", t, rho*ke/2., rho*g_*gpe);
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


// event profiles (t += 0.05)
// {
//   double ke = 0., gpe = 0.;
//   foreach (reduction(+:ke) reduction(+:gpe)) {
//     double zc = zb[];
//     foreach_layer() {
//       double norm2 = sq(w[]);
//       foreach_dimension()
// 	      norm2 += sq(u.x[]);
//       ke += norm2*h[]*dv();
//       gpe += (zc + h[]/2.)*h[]*dv();
//       zc += h[];
//     }
//   }
//   static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
//   if (i == 0)
//     fprintf (fp, "set term x11\n");
//   plot_profile (t, fp);
// }

// event gnuplot (t = Tend) {
//   FILE * fp = popen ("gnuplot", "w");
//   fprintf (fp,
//            "set term pngcairo enhanced size 640,640 font \",8\"\n"
//            "set output 'snapshot.png'\n");
//   plot_profile (t, fp);
// }

event output_field (t <= Tend; t+=.05)
{
    fprintf(stdout, "field vts output at step: %d, time: %.2f \n", i, t);
    static int j = 0;
    char name[100];
    sprintf(name, "%s/field_%.6i.vts", vts_folder, j++);
    fprintf(stdout, "written to: %s\n", name);
    FILE* fp = fopen(name, "w");
    output_vts_ascii_all_layers(fp, {eta,h,u,w}, N);
    fclose(fp);
    #if _OPENMP
    omp_set_num_threads(6);
    #endif
}
//gauges to compare the surface elevation
Gauge gauges[] = {
  {guage_name,  0},
  {NULL}
};


event output (t +=0.01 ; t <= Tend)
  output_gauges (gauges, {eta});
