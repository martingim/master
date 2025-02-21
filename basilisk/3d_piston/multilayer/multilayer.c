/**

based on bar.c and breaking.c
the parameters for the wave are from 
  Jensen, A., Sveen, J. K., Grue, J., Richon, J. B., & Gray, C. (2001). 
  Accelerations in water waves by extended particle image velocimetry. 
  Experiments in Fluids, 30(5), 500–510. https://doi.org/10.1007/s003480000229
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/check_eta.h"
#include "layered/perfs.h"

#include "view.h"
#if _OPENMP
#include "output_pvd.h"
#endif
//#include "display.h"
double Tend = 20;    //the end time of the simulation
double Lx = 7; //The length of the simulation domain
double g = 9.81;
int LEVEL = 8;      //the grid resolution in x direction Nx = 2**LEVEL
double rmin = 0.5;  //rmin the relative height of the top layer compared to 
                    //a regular distribution. the rest fo the layers follow a geometric distribution.
double h_ = 0.6;                   
char save_location[] = "./";

int set_n_threads = 2; //set number of threads manually
#define nl_ 10  //the default number of layers if none are given as command line arguments
#define g_ g

double t_dim = 1;

char results_folder[40]; //the location to save the results
char vts_folder[50]; //the locaton to save the vtu files

//piston file 
int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/fil3.dat";
//piston parameters 
double piston_width = 0.5;
int file_samplerate = 100; //the samplerate of the piston speed file
#define n_pistons 14 //number of pistons in the tank
#define piston_timesteps 10000//the max number of timesteps in the piston file
char piston_folder[40];
double piston_ux[n_pistons][piston_timesteps];

void read_piston_data(){
  char piston_file[60];
  for (int n=0;n<n_pistons;n++){
    int count = 0;
    FILE *file;
    sprintf(piston_file, "%spiston_speed_%02d.dat", piston_folder, n);
    printf("opening piston file:%s\n", piston_file);

    file = fopen(piston_file, "r");
    if(!file)
      {
        perror("Error opening piston position file");
      }

      int _running=1;
      while(_running && count< piston_timesteps ){
        _running = fscanf(file, "%lf", &(piston_ux[n][count]));
        count++;
      }
      while (count<piston_timesteps){
        piston_ux[n][count]=0;
        count++;
      }
    
    fclose(file);
  }
  for (int n=0;n<piston_timesteps;n++){
    piston_ux[4][n] = -piston_ux[4][n];
  }
}


double Wave_VeloX(double x, double y, double z, double t){
  int piston_number = 0;
  piston_number = (int) floor(y/piston_width);
  return piston_ux[piston_number][(int) floor(t*file_samplerate)];
  
}

event init (i = 0)
{
  #if _OPENMP
  fprintf(stderr, "current threads:%d\n", omp_get_num_threads());
  fprintf(stderr, "max number of openmp threads:%d\n", omp_get_max_threads());
  #endif
  // u.n[right] = dirichlet(0);

  geometric_beta(rmin, true); //set the layer thickness smaller nearer the surface
  foreach() {
    zb[] = -h_;
    foreach_layer(){
      h[] = (max(- zb[], 0.))*beta[point.l];
    }
  }
  //# define neumann_pressure_variable(i) ((neumann_pressure_function(t)))
  # define neumann_pressure_variable(i) (((Wave_VeloX(x, y, z, t+dt)-Wave_VeloX(x, y, z, t))/dt))
  u.n[left] = dirichlet(Wave_VeloX(x,y,z,t));
  phi[left] = neumann(neumann_pressure_variable(0));
  
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
    if (strcmp(argv[j], "-r") == 0) 
        {                 
            run_number = atoi(argv[j + 1]);
        }  
    if (strcmp(argv[j], "-nl") == 0)
    {                 
      nl = atoi(argv[j + 1]); 
    }
  }
  #if _MPI
  dimensions(3);
  #endif
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
  
  //piston data
  sprintf(piston_folder, "piston_files/%d/", run_number);
  printf("%s\n", piston_folder);
  read_piston_data();
  //copy the script to the results folder for later inspection if needed
  char copy_script[100];
  sprintf(copy_script, "cp multilayer.c %s/multilayer.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }

  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  breaking = 0.1;
  CFL_H = .5;
  TOLERANCE = 10e-5;
  DT = 0.01;
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
event output_field (t <= Tend; t += 1)
{
    fprintf(stdout, "field vts output at step: %d, time: %.2f \n", i, t);
    static int j;
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
/**
event movie (t += 1./25., t <= Tend)
{
  view (fov = 17.3106,quat = {0.549, -0.058, -0.101, 0.828},
      tx = -0.356, ty = -0.266, tz = -1.720,
      width = 1200, height = 768);
  char s[80];
  sprintf (s, "t = %.2f T0", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.05, max = 0.06);
  save ("movie.mp4");
}
**/
/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=Tend$.*/
#if 0
void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [0.0:%f][-0.8:0.6]'-' u 1:3:2 w filledcu lc 3 t ''\n", t, Lx);
  foreach(serial)
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
  fprintf (fp,pstn1[] =  piston_function(x, y, 0., 0.5, Delta);
  
}

#endif

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%g, ", t, i, dt);
  printf("%.2f%%\n", progress*100);
}


//gauges to compare the surface elevation
Gauge gauges[] = {
  {"X_0",  1.5, 3.5},
  {"X_1",  10.048, 3.5},
  {"X_2",  10.745, 3.5},
  {"X_3",  11.498, 3.5},
  {NULL}
  };


event output (t += 0.01; t <= Tend)
 output_gauges (gauges, {eta});