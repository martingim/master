/*
 * sets the left boundary condition to the speed of the piston movement
 * reads piston position data from file. Set up for file with 100hz sample rate 
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "adapt_wavelet_leave_interface.h"
#include "heights.h"
#include "output.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
//#include "tension.h"
#include "embed.h"
#include "reduced.h"
#include "profiling.h"
#include "output_vtu_foreach.h"

int set_n_threads = 1; //0 to use all available threads for OPENMP
int LEVEL = 6;
int max_LEVEL = 12; //Default level if none is given as command line argument
int padding = 2;

#define _h 0.6//water depth
double l = 25.6/2; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.01;
double uemax = 0.01;
double pemax = .01;
double Tend = 10.;

double probe_positions[144];
int n_probes = 144;
int n_extra_probe = 4;
double extra_probe_positions[]={8.00, 10.04, 10.75, 11.50};

vector h[]; //scalar field of the distance from the surface, using heights.h
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files

//piston file 
int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/fil3.dat";
char piston_file[40];
int file_samplerate = 125; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
int piston_counter;
double piston_ux[piston_timesteps];
double U_X = 0.; //the speed of the piston
//piston parameters 

void read_piston_data(){
  int count = 0;
  double piston_positions[piston_timesteps];
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
  for (int i=0;i<count-1;i++){
    piston_ux[i] = (piston_positions[i+1]-piston_positions[i])*file_samplerate;
  }
}

event setup_probe_positions(i=0){
  int j = 0;
  for(j=0;j<n_extra_probe;j++){
    probe_positions[j] = extra_probe_positions[j];
  }
  probe_positions[j] = 0.1;
  j++;
  for(j=j;j<n_probes;j++){
    probe_positions[j] = probe_positions[j-1]+0.1;
  }
}

void mask_domain(){
  //mask away the top of the domain
  mask(y > domain_height - _h ? top : none);
  u.n[top]  = neumann(0.);
  p[top]    = dirichlet(0.);
  pf[top]   = dirichlet(0.);
}

int main(int argc, char *argv[]) {
  
  //set max_LEVEL and run number from command line args
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) // This is your parameter name
        {                 
            max_LEVEL = atoi(argv[j + 1]);    // The next value in the array is your value
        }
    if (strcmp(argv[j], "-r") == 0) // This is your parameter name
        {                 
            run_number = atoi(argv[j + 1]);    // The next value in the array is your value
        }  
  }

  //make folders for saving the results
  sprintf(results_folder, "results/run%d/LEVEL%d", run_number, max_LEVEL);
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
  
  
  //copy the script to the results folder for later incpection if needed
  char copy_script[100];
  sprintf(copy_script, "cp wave_wall.c %s/wave_wall.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }

  //piston data
  sprintf(piston_file, "piston_files/%d/fil3.dat", run_number);
  printf("%s\n", piston_file);
  read_piston_data();
  
  L0 = l;
  //f.sigma = 0.078;
  rho1 = 997;
  rho2 = 1.204;
  mu1 = 8.9e-4;
  mu2 = 17.4e-6;
  G.y = - 9.81;
  N = 1 << LEVEL;
  DT = 0.01;
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
  u.n[left] = dirichlet(0.);
  u.n[left] = neumann(0.);
  u.n[right] = dirichlet(0.);
  u.n[top] = neumann(0.);
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

event init (i = 0) {
  origin(0, -_h);
  mask_domain();
  //refine((fabs(y)<l/N*0.49)&&(level<=max_LEVEL));
  //mvtu(42);
  fraction (f, - y); //set the water depth _h
  while (adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax,uemax,pemax, femax},max_LEVEL, LEVEL,padding).nf){  //for adapting more around the piston interface
    fraction (f, - y); //set the water level on the refined mesh
  }
  unrefine ((x < -0.1)&&(level>6));
  foreach(){
    pf[] = 0;
    p[] = pf[];
  }
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax, uemax, pemax, femax}, max_LEVEL, LEVEL,padding);
  unrefine ((y>0.1));//unrefine the air
}


/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
*/
event piston (i++, first) {
  piston_counter = floor(t*file_samplerate);
  //printf("t:%f, file_timestep:%d, %%to next file timestep:%.0f%%, piston_position:%f\n", t, piston_counter, counter_remainder*100, piston_position);
  U_X = piston_ux[piston_counter];
  u.n[left] = dirichlet(U_X);
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
  heights(f, h);
  double min_height; //vertical distance from the center of the cell to the air water interface
  double surface_elevation=0;
  for (int probe_n=0;probe_n<n_probes;probe_n++){
    min_height=1;
    foreach(serial){
      if((fabs(x-probe_positions[probe_n])<Delta/2.0)&&(fabs(y)<0.1)){
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

//save unordered mesh
event vtu(t+=0.1, last){
  printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);
}


// void mvtu(int s){
//   printf("Saving vtu file\n");
//   char filename[40];
//   sprintf(filename, "%s/vtu/TIME-%d", results_folder, s);
//   output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
// }

// event vtu(i++){
//   printf("Saving vtu file\n");
//   char filename[40];
//   sprintf(filename, "%s/vtu/step-%05d", results_folder, i);
//   output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
// }
// event move_results(t=0)
// {
//   char results_folder[40];
//   char make_folder[100];
//   sprintf(results_folder, "./results/run%d/LEVEL%d", run_number, max_LEVEL); 
//   sprintf(remove_folder, "rm -r %s", results_folder); 
//   sprintf(make_folder, "mkdir -p %s", results_folder);
  


//   if (system(remove_folder)==0){ 
//     printf("removed previous run results folder\n");
//   }
//   if (system(make_folder)==0){
//     printf("made folder for results:%s\n", results_folder);
//   }

  
  //mkdir(results_folder, 0755);
  //system(sprintf("rm -r %s", results_folder));//remove the results folder if it already exists
  //mkdir(results_folder,0755);
// }

/*
simulation stopped at Tend
*/
event stop (t = Tend);

event show_progress(i++)
{
  printf("t=%02.3f, i=%04d, dt=%.3g\n", t, i, dt);
  //float progress = 0;
  //progress = t /Tend;
  //printf("%.2f%%\r", progress*100);
}


event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
