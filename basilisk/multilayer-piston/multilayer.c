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
#include "output_pvd.h"
double Tend = 30;    //the end time of the simulation
double Lx = 7; //The length of the simulation domain
double g = 9.81;
int LEVEL = 7;      //the grid resolution in x direction Nx = 2**LEVEL
double rmin = 0.5;  //rmin the relative height of the top layer compared to 
                    //a regular distribution. the rest fo the layers follow a geometric distribution.
double h_ = 0.6;                   
char save_location[] = "./";

int set_n_threads = 6; //set number of threads manually
#define nl_ 5  //the default number of layers if none are given as command line arguments
#define g_ g

double t_dim = 1;

char results_folder[40]; //the location to save the results
char vts_folder[50]; //the locaton to save the vtu files

float piston_function(double y, double y0, double y1, double Delta){
    if ((y<y0-Delta/2.)||(y>y1+Delta/2.)){
      return 0.;
    }
    else if ((y-y0)<Delta/2.) {
      return 0.5-(y0-y)/Delta;
    }
    else if ((y1-y)<Delta/2.) {
      return 0.5-(y-y1)/Delta;
    }
    else{
      return 1.;
    }
}


//piston file 
int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/fil3.dat";
char piston_file[40];
int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
int piston_counter;
double piston_positions[piston_timesteps];
double piston_position = 0; //the starting position of the piston
double piston_position_p; //piston position at the previous timestep
double U_X = 0.; //the speed of the piston
//piston parameters 

scalar pstn1[] = 0;
scalar pstn2[] = 0;
scalar pstn3[] = 0;
scalar pstn4[] = 0;
scalar pstn5[] = 0;
scalar pstn6[] = 0;
scalar pstn7[] = 0;
scalar pstn8[] = 0;
scalar pstn9[] = 0;
scalar pstn10[] = 0;
scalar pstn11[] = 0;
scalar pstn12[] = 0;
scalar pstn13[] = 0;
scalar pstn14[] = 0;

void read_piston_data(){
  int count = 0;
  FILE *file;
  file = fopen(piston_file, "r");
  if(!file)
    {
        perror("Error opening piston position file");
    }
  int _running=1;
  while(_running && count< piston_timesteps ){
    _running = fscanf(file, "%lf", &(piston_positions[count]));
    piston_positions[count] /=100.; //convert to meters
    count++;
  }
  fclose(file);
  double start_offset = piston_positions[0];
  for (int i=0;i<count;i++){
    piston_positions[i] -= start_offset;
  }
  piston_position_p = piston_positions[0];
  piston_position = piston_positions[0];
}

event piston_update(i++){
  piston_counter = floor(t*100);
  double ins_ = 1;
  double counter_remainder = 0;
  counter_remainder = t*100.*ins_-piston_counter;
  piston_position = piston_positions[piston_counter] + (piston_positions[piston_counter+1] -piston_positions[piston_counter])*counter_remainder; //update the piston position
  //printf("t:%f, file_timestep:%d, %%to next file timestep:%.0f%%, piston_position:%f\n", t, piston_counter, counter_remainder*100, piston_position);
  U_X = (piston_position-piston_position_p)/dt;
  piston_position_p = piston_position;
  u.n[left] = dirichlet(U_X*(pstn1[]+pstn2[]-pstn3[]+pstn12[])); 
  //u.n[left] = dirichlet(U_X);
}

event init (i = 0)
{
  fprintf(stderr, "current threads:%d\n", omp_get_num_threads());
  fprintf(stderr, "max number of openmp threads:%d\n", omp_get_max_threads());
  u.n[right] = dirichlet(0);

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

foreach(){
  pstn1[] =  piston_function(y, 0., 0.5, Delta);
  pstn2[] =  piston_function(y, 0.5, 1., Delta);
  pstn3[] =  piston_function(y, 1., 1.5, Delta);
  pstn4[] =  piston_function(y, 1.5, 2., Delta);
  pstn5[] =  piston_function(y, 2., 2.5, Delta);
  pstn6[] =  piston_function(y, 2.5, 3., Delta);
  pstn7[] =  piston_function(y, 3., 3.5, Delta);
  pstn8[] =  piston_function(y, 3.5, 4., Delta);
  pstn9[] =  piston_function(y, 4., 4.5, Delta);
  pstn10[] =  piston_function(y, 4.5, 5., Delta);
  pstn11[] =  piston_function(y, 5., 5.5, Delta);
  pstn12[] =  piston_function(y, 5.5, 6., Delta);
  pstn13[] =  piston_function(y, 6., 6.5, Delta);
  pstn14[] =  piston_function(y, 6.5, 7., Delta);
}

}

int main(int argc, char *argv[])
  {

    //set max_LEVEL and run number from command line args
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) // This is your parameter name
        {                 
            LEVEL = atoi(argv[j + 1]);    // The next value in the array is your value
        }
    if (strcmp(argv[j], "-r") == 0) // This is your parameter name
        {                 
            run_number = atoi(argv[j + 1]);    // The next value in the array is your value
        }  
    if (strcmp(argv[j], "-nl") == 0) // This is your parameter name
        {                 
            nl = atoi(argv[j + 1]);    // The next value in the array is your value
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

  //piston data
  sprintf(piston_file, "piston_files/%d/fil3.dat", run_number);
  printf("%s\n", piston_file);
  read_piston_data();
  //origin(0,-0.5);

  if (argc>1)
    LEVEL = atoi(argv[1]);
    if (argc>2)
      nl = atoi(argv[2]);
    else
      nl = nl_;
  
  //origin (0, h_);
  N = 1<<LEVEL;
  L0 = Lx;
  G = g;
  breaking = 0.1;
  CFL_H = .5;
  TOLERANCE = 10e-5;
  //theta_H = 0.51;
  #if _OPENMP
  int num_omp = omp_get_max_threads();
  fprintf(stderr, "max number of openmp threads:%d\n", num_omp);
  omp_set_num_threads(set_n_threads);   
  fprintf(stderr, "set openmp threads:%d\n", 4);
  #elif _MPI
  fprintf(stderr, "mpiiiiii\n");
  fprintf(stderr, "npe:%d\n", npe());
  #endif

  run();
}


event output_field (t <= Tend; t += 1)
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
  fprintf (fp,pstn1[] =  piston_function(y, 0., 0.5, Delta);
  
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
  {"X_0",  8.009},
  {"X_1",  10.048},
  {"X_2",  10.745},
  {"X_3",  11.498},
  {NULL}
  };


//event output (t += 0.01; t <= Tend)
//  output_gauges (gauges, {eta});
