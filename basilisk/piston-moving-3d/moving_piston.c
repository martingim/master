/*
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 * can read piston position data from file. Set up for file with 100hz sample rate 
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "grid/octree.h"
#include "utils.h"
#include "adapt_wavelet_leave_interface_two_levels.h"
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

int set_n_threads = 6; //0 to use all available threads for OPENMP
int LEVEL = 4;
int max_LEVEL = 7; //Default level if none is given as command line argument
int padding = 2;
int EXTRA_PISTON_LEVEL = 0; //extra refinement around the piston to make it leak less 

#define _h 0.6//water depth
double l = 7; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double domain_height = 1.0; //the height of the simulation domain
double domain_length = 6;
double domain_width = 7;
double femax = 0.2;
double uemax = 0.2;
double pemax = .2;
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
int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
#define n_pistons 14
double piston_positions[piston_timesteps][n_pistons];

//piston parameters 
#define piston_back_wall_offset 1. //the distance from the left of the domain to the front of the piston
double piston_height = 0.2;   //height of the piston above the still water level
double piston_bottom_clearance = 0.00; //piston distance above bottom
double piston_position[n_pistons];
double piston_position_p[n_pistons];
#define piston_depth .2 //depth of the piston
#define piston_width 0.5
//#define PISTON (piston_width/2.-fabs(z-(piston_number+0.5)*piston_width)) -(piston_position[piston_number]<x) - (x<(piston_position[piston_number]-piston_depth))-(y>piston_height) //subtract conditions outside the piston
#define PISTON intersection(intersection(piston_width/2.-fabs(z-(piston_number+0.5)*piston_width), piston_position[piston_number]-x), x-(piston_position[piston_number]-piston_depth))
scalar pstn[],pistons[];


void read_piston_data(){
  //This only take one piston position file as is
  int count = 0;
  FILE *file;
  file = fopen(piston_file, "r");
  if(!file)
    {
        perror("Error opening piston position file");
    }
  int _running=1;
  while(_running && count< piston_timesteps ){
    _running = fscanf(file, "%lf", &(piston_positions[count][0]));
    piston_positions[count][0] *=0.044; //convert to meters
    count++;
  }
  fclose(file);
  double start_offset = piston_positions[0][0];
  for (int i=0;i<count;i++){
    piston_positions[i][0] -= start_offset;
  }
  piston_position_p[0] = piston_positions[0][0];
  piston_position[0] = piston_positions[0][0];
}

event setup_probe_positions(i=0){
  //should add second dimension here
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
  mask(x > domain_length ? right: none);
  mask(z > domain_width ? back : none);

  u.n[top]  = neumann(0.);
  p[top]    = dirichlet(0.);
  pf[top]   = dirichlet(0.);
  u.n[back] = dirichlet(0.);
  u.n[front] = dirichlet(0.);
  u.n[right] = neumann(0.);
  u.n[left] = neumann(0.);
}

int main(int argc, char *argv[]) {
  
  //set max_LEVEL and run number from command line args
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0) // This is your parameter name
        {                 
            max_LEVEL = atoi(argv[j + 1]);    // The next value in the array is your value
        }
    if (strcmp(argv[j], "-P") == 0) // This is your parameter name
        {                 
            EXTRA_PISTON_LEVEL = atoi(argv[j + 1]);    // The next value in the array is your value
        }
    if (strcmp(argv[j], "-r") == 0) // This is your parameter name
        {                 
            run_number = atoi(argv[j + 1]);    // The next value in the array is your value
        }  
  }
  
  pstn.refine = pstn.prolongation = fraction_refine;
  pistons.refine = pistons.prolongation = fraction_refine;
  
  //make folders for saving the results
  sprintf(results_folder, "results/run%d/LEVEL%d_%d", run_number, max_LEVEL, EXTRA_PISTON_LEVEL);
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
  sprintf(piston_file, "piston_files/%d/padle_ut.dat", run_number);
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
  origin(-piston_back_wall_offset, -_h);
  mask_domain();
  //refine((fabs(y)<l/N*0.49)&&(level<=max_LEVEL));
  //mvtu(42);
  fraction (f, - y); //set the water depth _h
  
  for (int piston_number=0;piston_number<n_pistons;piston_number++){
    fraction (pstn, PISTON); //set the piston fraction
    foreach(){
      pistons[]+=pstn[];
    }
    
  }
  while (adapt_wavelet_leave_interface({u.x, u.y, p},{f, pistons},(double[]){uemax,uemax,femax,pemax, femax}, max_LEVEL+EXTRA_PISTON_LEVEL, LEVEL,padding, (int[]){max_LEVEL, max_LEVEL+EXTRA_PISTON_LEVEL}).nf){  //for adapting more around the piston interface
    fraction (f, - y); //set the water level on the refined mesh
    
    foreach(){
      pistons[]= 0;
    }
    for (int piston_number=0;piston_number<n_pistons;piston_number++){
      fraction (pstn, PISTON); //set the piston fraction
      foreach(){
        pistons[]+=pstn[];
      }
    }
  }
  //unrefine ((x < -0.1)&&(level>6));
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
  adapt_wavelet_leave_interface({u.x, u.y, p},{f, pistons},(double[]){uemax, uemax, pemax, femax, pemax}, max_LEVEL+EXTRA_PISTON_LEVEL, LEVEL,padding, (int[]){max_LEVEL, max_LEVEL+EXTRA_PISTON_LEVEL});
  //unrefine ((x < piston_position-piston_w*0.6)); //unrefine the area to the left of the piston
  //unrefine ((y<-0.4)&&(x>(piston_position+0.02))); //unrefine the bottom
  //unrefine ((y>0.1));
  
  // //set piston fraction after refining mesh
  // foreach(){
  //   pistons[]= 0;
  // }
  // for (int piston_number=0;piston_number<n_pistons;piston_number++){
  //   fraction (pstn, PISTON); //set the piston fraction
  //   foreach(){
  //     pistons[]+=pstn[];
  //   }
  // }

}


/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
*/
int piston_counter;
event piston (i++, first) {
  double U_X = 0;
  piston_counter = floor(t*file_samplerate);
  double counter_remainder = 0;
  counter_remainder = t*file_samplerate-piston_counter;
  
  for (int piston_number=0;piston_number<n_pistons;piston_number++){
    fraction (pstn, PISTON); //set the piston fraction
    piston_position[piston_number] = piston_positions[piston_counter][piston_number] + (piston_positions[piston_counter+1][piston_number] -piston_positions[piston_counter][piston_number])*counter_remainder; //update the piston position
    //printf("t:%f, file_timestep:%d, %%to next file timestep:%.0f%%, piston_position:%f\n", t, piston_counter, counter_remainder*100, piston_position);
    U_X = (piston_position[piston_number]-piston_position_p[piston_number])/dt;
    piston_position_p[piston_number] = piston_position[piston_number];
    foreach(){
      u.y[] = u.y[]*(1 - pstn[]);
      u.x[] = pstn[]*U_X + u.x[]*(1 - pstn[]);
    }
  }
  //printf("U_X:%f, piston_position:%f\n", U_X, piston_position);
}

event surface_probes(t+=0.1){
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
  double surface_elevation;
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
event vtu(t+=.1, last){
  fprintf(stderr, "openmp threads:%d\n", omp_get_max_threads());
  printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p,pistons}, (vector *) {u}, filename);
  fprintf(stderr, "openmp threads:%d\n", omp_get_max_threads());
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


event profiling (i += 1) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
