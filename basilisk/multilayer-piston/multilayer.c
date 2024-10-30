/**

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
//#include "view.h"
#include "output_pvd.h"

//#include "display.h"
double Tend = 10;    //the end time of the simulation
double Lx = 24.6; //The length of the simulation domain
double g = 9.81;
int LEVEL = 9;      //the grid resolution in x direction Nx = 2**LEVEL
double rmin = 0.5;  //rmin the relative height of the top layer compared to 
                    //a regular distribution. the rest fo the layers follow a geometric distribution.
double h_ = 0.6;    //water depth
int set_n_threads = 6; //set number of threads manually
#define nl_ 10  //the default number of layers if none are given as command line arguments
#define g_ g

char results_folder[40]; //the location to save the results
char vts_folder[50]; //the locaton to save the vtu files
//surface probes
double probe_positions[144] = {8.0,10.048,10.745,11.498};
int n_probes = 144;

event setup_probe_positions(i=0){
  double probe_x = 0.1;
  for(int j=4;j<n_probes;j++){
    probe_positions[j] =probe_x;
    probe_x = probe_x+0.1;
  }
}

//piston parameters
int run_number = 1;
double piston_f = 1.425;//the frequency of the piston movement
double piston_amplitude = 0.042*0.308;
double U_X = 0.; //the speed of the piston

event piston_update(i++){
  if (t<0.1){
    U_X = 0;
  }
  else{
  U_X = piston_amplitude*(4/cosh(2*t-0.2)/cosh(2*t-0.2)*tanh(2*t-0.2)*sin(2*pi*piston_f*t-0.34) + 2*piston_f*pi*cos(2*pi*piston_f*t-0.34)*tanh(2*t-0.2)*tanh(2*t-0.2));
  }
  u.n[left] = dirichlet(U_X);
}

event init (i = 0)
{
  #if _OPENMP
  fprintf(stderr, "current threads:%d\n", omp_get_num_threads());
  fprintf(stderr, "max number of openmp threads:%d\n", omp_get_max_threads());
  #endif
  // u.n[right] = dirichlet(0);
  //u.n[left] = neumann(0);
  geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
  foreach() {
    zb[] = -h_;
    foreach_layer(){
      h[] = (max(- zb[], 0.))*beta[point.l];
    }
  }
  
  // foreach() {
  //  zb[] = -h_;
  //  foreach_layer(){
  //    h[] = h_/nl;
  //  }
  // }

}

int main(int argc, char *argv[])
  {
    //set LEVEL, layers and run number from command line args
  nl = nl_;
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) 
        {                 
            LEVEL = atoi(argv[j + 1]);
        }
    if (strcmp(argv[j], "-nl") == 0)
    {                
      nl = atoi(argv[j + 1]); 
    }
  }
  //make folders for saving the results
  sprintf(results_folder, "results/run%d/LEVEL%d_layers%d", run_number, LEVEL, nl);
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
  
  
  //copy the script to the results folder for later incpection if needed
  char copy_script[100];
  sprintf(copy_script, "cp multilayer.c %s/multilayer.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }

  origin (0, h_);
  
  //origin (0, h_);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  breaking = 0.1;
  CFL_H = .5;
  TOLERANCE = 10e-5;
  //theta_H = 0.51;
  #if _OPENMP
  if (set_n_threads>0)
  {
  int num_omp = omp_get_max_threads();
  fprintf(stderr, "max number of openmp threads:%d\n", num_omp);
  omp_set_num_threads(set_n_threads);   
  fprintf(stderr, "set openmp threads:%d\n", 4);
  }
  #elif _MPI
  fprintf(stderr, "npe:%d\n", npe());
  #endif

  run();
}

#if _OPENMP
event output_field (t <= Tend; t += .1)
{
    fprintf(stdout, "field vts output at step: %d, time: %.2f \n", i, t);
    static int j = 0;
    char name[100];
    sprintf(name, "%s/field_%.6i.vts",vts_folder, j++);
    fprintf(stdout, "written to: %s\n", name);
    FILE* fp = fopen(name, "w");
    output_vts_ascii_all_layers(fp, {eta,h,u}, N);
    fclose(fp);
    #if _OPENMP
    omp_set_num_threads(set_n_threads);
    #endif
}
#endif
// event movie (t += 1./25., t <= Tend)
// {
//   view (fov = 17.3106,quat = {0.549, -0.058, -0.101, 0.828},
//       tx = -0.356, ty = -0.266, tz = -1.720,
//       width = 1200, height = 768);
//   char s[80];
//   sprintf (s, "t = %.2f T0", t);
//   draw_string (s, size = 80);
//   squares ("eta", linear = true, z = "eta", min = -0.05, max = 0.06);
//   save ("movie.mp4");
// }

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=Tend$.*/
#if 1
// void plot_profile (double t, FILE * fp)
// {
//   fprintf (fp,
// 	   "set title 't = %.2f'\n"
// 	   "p [0.0:%f][-0.8:0.6]'-' u 1:3:2 w filledcu lc 3 t ''\n", t, Lx);
//   foreach(serial)
//     fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
//   fprintf (fp, "e\n\n");
//   fflush (fp);
// }


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
//   fprintf (fp,pstn1[] =  piston_function(x, y, 0., 0.5, Delta);  
// }

#endif

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%g, ", t, i, dt);
  printf("%.2f%%\n", progress*100);
}

event surface_probes(t+=0.01){
  char filename[100];
  sprintf(filename, "%s/surface_probes.csv", results_folder);
  FILE *fp = fopen(filename, "a");
    if (i==0){
      fprintf(fp, "time, U_X");
      for (int j = 0;j<n_probes;j++){
        fprintf(fp, ", %f", probe_positions[j]);    
      }
      fprintf(fp, "\n");
  }
  fprintf(fp, "%f, %f", t, U_X);
  for (int probe_n=0;probe_n<n_probes;probe_n++){
    
    fprintf(fp, ", %f", interpolate(eta, probe_positions[probe_n],-h_));
  }
  fprintf(fp, "\n");
  fclose(fp);
}