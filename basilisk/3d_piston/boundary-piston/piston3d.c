/*
 * sets the left boundary condition to the speed of the piston movement
 * reads piston position data from file. Set up for file with 100hz sample rate 
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "grid/octree.h"
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
int LEVEL = 4;
int max_LEVEL = 8; //Default level if none is given as command line argument
int padding = 0;

int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/fil3.dat";
#define symmetric_ 0 //if the simulation is symmetric around z=domain_width/2

#define n_pistons 14//number of pistons
#define _h 0.5//water depth
double l = 7; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double domain_length = 7;
double domain_width = 7; //the width of the simulation domain
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.2;
double uemax = 0.2;
double pemax = 0.2;
double Tend = 10.;

double probe_positions[144];
int n_probes = 144;
int n_extra_probe = 4;
double extra_probe_positions[]={8.00, 10.04, 10.75, 11.50};

vector h[]; //scalar field of the distance from the surface, using heights.h
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files

//piston file 

int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the max number of timesteps in the piston file

char piston_folder[40];
double piston_width = 0.5;
double piston_ux[n_pistons][piston_timesteps];
double ramp_start = 0.1; //ramp for the pressure boundary condition
double ramp_distance = 0.1;
# define neumann_pressure_variable(i) ((paddle_rampoff(y)*(Wave_VeloX(0,0,z,t+dt)-Wave_VeloX(0,0,z,t))/dt - a.n[i])*fm.n[i]/alpha.n[i])

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

double paddle_rampoff(double y){
  double ramp = 1;
  if (y<ramp_start){
    ramp = 1;
  }
  if (y>ramp_start){
    ramp = (ramp_start+ramp_distance-y)/ramp_distance;
    if (y>ramp_start+ramp_distance){
      ramp = 0;
    }
  }
  return ramp;
}

double Wave_VeloX(double x, double y, double z, double t){
  int piston_number = 0;
  piston_number = (int) floor(z/piston_width);
  return piston_ux[piston_number][(int) floor(t*file_samplerate)];
}

void mask_domain(){
  //mask away the top and side of the domain 
  mask(y > domain_height - _h ? top : none);
  mask(x > domain_length ? right: none);
  #if symmetric_
  mask(z > domain_width/2. ? back : none);
  #else
  mask(z > domain_width ? back : none);
  u.n[front]= dirichlet(0.);
  #endif
  u.n[back] = dirichlet(0.);
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
  sprintf(copy_script, "cp piston3d.c %s/piston3d.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }

  //piston data
  sprintf(piston_folder, "piston_files/%d/", run_number);
  printf("%s\n", piston_folder);
  read_piston_data();

  L0 = l;
  //f.sigma = 0.078;
  rho1 = 997;
  rho2 = 1.204;
  mu1 = 8.9e-4;
  mu2 = 17.4e-6;
  G.y = - 9.81;
  N = 1 << LEVEL;
  DT = 0.05;
  
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
  origin(0., -_h, 0.);
  mask_domain();
  u.n[top]  = neumann(0.);
  p[top]    = dirichlet(0.);
  u.t[top]   = dirichlet(0.);

  u.n[left] = dirichlet(Wave_VeloX(x,y,z,t)*paddle_rampoff(y));
  p[left] = neumann(neumann_pressure_variable(0));

  //u.n[bottom] = dirichlet(0.);
  //u.t[bottom] = dirichlet(0.);
  u.n[right] = neumann(0.);

  fraction (f, - y); //set the water depth _h
  while (adapt_wavelet_leave_interface({u.x, u.y, u.z, p},{f},(double[]){uemax,uemax,uemax, pemax, femax},max_LEVEL, LEVEL,padding).nf){  //for adapting more around the piston interface
    fraction (f, - y); //set the water level on the refined mesh
    // int piston_i=0;
    // for (scalar pstn in pistons){
    //   fraction(pstn, piston_width/2.-fabs(z-(piston_i+0.5)*piston_width));
    //   piston_i +=1;
    // }
  }

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
  adapt_wavelet_leave_interface({u.x, u.y, u.z, p},{f},(double[]){uemax, uemax, uemax, pemax, femax}, max_LEVEL, LEVEL,padding);
  //unrefine ((y>0.1));//unrefine the air
}


/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
*/
// double piston_function(double t_){
//   return 0.013552*(cos(8.95354*t_)/cosh(t_)/cosh(t_) - 8.95354*sin(8.95354*t_)*tanh(t_));
// }
// event piston (i++, first) {
//   foreach(){
//     pistn[]=1;
//   }
//   // double t_ = 0;
//   // for (scalar pstn in pistons){
//   //   foreach(){
//   //     pistn[]+=pstn[]*piston_function(t+t_);
//   //   }
//   //   t_+=0.05;

//   // }
//   //u.n[left] = dirichlet(pstn1_u[(int)floor(t*file_samplerate)]*pstn1[]+pstn1_u[pui]*pstn2[]);
//   // u.n[left] = dirichlet(pistn[]);
//   u.n[left] = dirichlet(pistn[]*pstn1_u[(int)floor(t*file_samplerate)]);
// }

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
  double surface_elevation=0;
  for (int probe_n=0;probe_n<n_probes;probe_n++){
    min_height=1;
    foreach(serial){
      if((fabs(x-probe_positions[probe_n])<Delta/2.0)&&(fabs(y)<0.1)&&(fabs(z-3.49)<Delta/2.)){
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
// event vtu(t+=0.1, last){
event vtu(i++){
  printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*1000));
  output_vtu((scalar *) {f,p}, (vector *) {u}, filename);

  //output all the piston fields
  //output_vtu((scalar *) {f,p,pstn1,pstn2,pstn3,pstn4,pstn5,pstn6,pstn7,pstn8,pstn9,pstn10,pstn11,pstn12,pstn13,pstn14}, (vector *) {u}, filename);
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
