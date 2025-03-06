/*
 * sets the left boundary condition to the speed of the piston movement
 * reads piston position data from file.
 * Thank you to Oystein Lande for how to implement the boundary conditions
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

int set_n_threads = 2; //0 to use all available threads for OPENMP
int LEVEL = 6;      //length_of domain/2**LEVEL the largest cells used
int max_LEVEL = 10; //Default level if none is given as command line argument length_of_domain/2**max_LEVEL the size of the smallest cells
int padding = 2;    //How many cells around the air water interface that should keep the highest refinement level

#define _h 0.6//water depth
double l = 14; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.01;
double uemax = 0.01;
double pemax = .01;
double Tend = 25.;

char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files

//surface probes
double probe_positions[144]; //surface probes that meausre the elevation of the free surface spaced 0.1m apart
int n_probes = 144;
int n_extra_probe = 4;
double extra_probe_positions[]={1.50, 10.04, 10.75, 11.50}; //extra surface probes at these positions(will be the first probes in the surfacs_probes file)
vector h[]; //scalar field for calculating the distance from the surface, using heights.h

//rampoff 
double ramp_start = 0.1;    //where the rampoff of the piston function starts. Meters above the still water level
double ramp_distance = 0.1; //the length of the rampoff zone

//piston file 
int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/piston_speed.dat";
int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
double piston_ux[piston_timesteps];

//pressure calculated from the piston movement for the boundary condition
# define neumann_pressure_variable(i) ((paddle_rampoff(y)*(Wave_VeloX(0,0,0,t+dt)-Wave_VeloX(0,0,0,t))/dt - a.n[i])*fm.n[i]/alpha.n[i])

void read_piston_data(){
  char piston_file[40];

  sprintf(piston_file, "piston_files/%d/piston_speed.dat", run_number);
  printf("piston file:%s\n", piston_file);

  int count = 0;
  FILE *file;
  file = fopen(piston_file, "r");
  if(!file)
    {
        perror("Error opening piston position file");
    }
  int _running=1;
  while(_running && count< piston_timesteps ){
    _running = fscanf(file, "%lf", &(piston_ux[count]));    
    count++;
  }
  fclose(file);
  for (int i=count;i<piston_timesteps;i++){
    piston_ux[i] = 0;
  }
}

double paddle_rampoff(double y){
  double ramp = 1;
  if (y>ramp_start){
    ramp = (ramp_start+ramp_distance-y)/ramp_distance;
    if (y>ramp_start+ramp_distance){
      ramp = 0;
    }
  }
  return ramp;
}

double Wave_VeloX(double x, double y, double z, double t){
  return piston_ux[(int)floor(t*file_samplerate)];
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

int main(int argc, char *argv[]) {
  
  //set max_LEVEL and run number from command line args
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) //check for command line flag
        {                 
            max_LEVEL = atoi(argv[j + 1]);
        }
    if (strcmp(argv[j], "-r") == 0)
        {                 
            run_number = atoi(argv[j + 1]);
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
  sprintf(copy_script, "cp piston.c %s/piston.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }
  
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


void mask_domain(){
  //mask away the top of the domain
  mask(y > domain_height - _h ? top : none);
}

event init (i = 0) {
  origin(0, -_h);
  mask_domain();
  // boundary conditions
  u.n[top]  = neumann(0.);
  p[top]    = dirichlet(0.);
  u.t[top]  = dirichlet(0.);
  
  u.n[left] = dirichlet(Wave_VeloX(0,0,0,t)*paddle_rampoff(y));
  p[left]   = neumann(neumann_pressure_variable(0));
  
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[right]  = neumann(0.);
  
  fraction (f, - y); //set the water depth _h
  while (adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax,uemax,pemax, femax},max_LEVEL, LEVEL,padding).nf){  //for adapting more around the piston interface
    fraction (f, - y); //set the water level on the refined mesh
  }
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y, p},{f},(double[]){uemax, uemax, pemax, femax}, max_LEVEL, LEVEL,padding);
  unrefine ((y>0.1)); //unrefine the air above the waves. make sure the waves are not so high that they get unrefined
                      //this is necessary to avoid problems when air vortexes reach the top boundary
}


event surface_probes(t+=0.01){
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
