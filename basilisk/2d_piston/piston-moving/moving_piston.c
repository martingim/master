/*
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 * can read piston position data from file. Set up for file with 100hz sample rate 
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "adapt_wavelet_leave_interface_two_levels.h"
#include "heights.h"
#include "output.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "profiling.h"
#include "output_vtu_foreach.h"

int set_n_threads = 8; //0 to use all available threads for OPENMP
int LEVEL = 4;
int max_LEVEL = 11; //Default level if none is given as command line argument
int padding = 2;
int EXTRA_PISTON_LEVEL = 0; //extra refinement around the piston to make it leak less 

#define _h 0.6//water depth
double l = 24.6; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double g_ = 9.81;
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.2;
double uemax = 0.2;
double pemax = 0.2;
double Tend = 50.;

double probe_positions[144];
int n_probes = 144;
int n_extra_probe = 4;
double extra_probe_positions[]={1.50, 10.04, 10.75, 11.50};

vector h[]; //scalar field of the distance from the surface, using heights.h
char results_folder[40]; //the location to save the results
char vtu_folder[50]; //the locaton to save the vtu files

//piston file 
int run_number = 1; //default run number if none is given in the command line piston files in "piston_files/%run_number/fil3.dat";
char piston_file[40];
char piston_speed_file[40];
int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
int piston_counter;
double piston_positions[piston_timesteps];
double piston_ux[piston_timesteps];
//piston parameters 
#define piston_back_wall_offset .1 //the distance from the left of the domain to the front of the piston
double piston_height = 0.2;   //height of the piston above the still water level
double piston_bottom_clearance = 0.00; //piston distance above bottom
#define piston_w .2 //width of the piston
#define PISTON Piston_Pos_x(x, y, z, t)-x
scalar pstn[];

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
    count++;
  }
  fclose(file);
  
  for (int i=count;i<piston_timesteps;i++){
    piston_positions[i] = 0;
  }
  count = 0;
  file = fopen(piston_speed_file, "r");
  if(!file)
    {
        perror("Error opening piston speed file");
    }
  _running=1;
  while(_running && count< piston_timesteps ){
    _running = fscanf(file, "%lf", &(piston_ux[count]));
    count++;
  }
  fclose(file);
  
  for (int i=count;i<piston_timesteps;i++){
    piston_ux[i] = 0;
  }
}

double Piston_Velo_x(double x , double y, double z, double t){
  int t0_i = (int)floor(t*file_samplerate);
  //linear interpolation of the piston speed f(t) = f(t0) +(f(t1)-f(t0))*(t-t0)/(t1-t0)
  return piston_ux[t0_i] + (piston_ux[t0_i+1]-piston_ux[t0_i])*(t*file_samplerate-t0_i);
}

double Piston_Pos_x(double x, double y, double z, double t){
  int t0_i = (int)floor(t*file_samplerate);
  return piston_positions[t0_i] + (piston_positions[t0_i+1]-piston_positions[t0_i])*(t*file_samplerate-t0_i);
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
  //set max_LEVEL, extra piston refinement
  //and run number from command line args
  for(int j=0;j<argc;j++){
    if (strcmp(argv[j], "-L") == 0)
        {
            max_LEVEL = atoi(argv[j + 1]);
        }
    if (strcmp(argv[j], "-P") == 0)
        {
            EXTRA_PISTON_LEVEL = atoi(argv[j + 1]);
        }
    if (strcmp(argv[j], "-r") == 0)
        {
            run_number = atoi(argv[j + 1]);
        }
  }

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

  //copy the script to the results folder for later inspection if needed
  char copy_script[100];
  sprintf(copy_script, "cp moving_piston.c %s/moving_piston.c", results_folder);
  if (system(copy_script)==0){
    printf("copied script to results folder\n");
  }

  sprintf(piston_file, "piston_files/%d/piston_position.dat", run_number);
  sprintf(piston_speed_file, "piston_files/%d/piston_speed.dat", run_number);

  printf("%s\n", piston_file);
  read_piston_data();

  L0 = l;
  rho1 = 997;
  rho2 = 1.204;
  mu1 = 8.9e-4;
  mu2 = 17.4e-6;
  G.y = - g_;
  N = 1 << LEVEL;
  CFL = 0.1;
  DT = l/(1<<(max_LEVEL));
  printf("DT=%.5f\n", DT);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
  u.n[left] = dirichlet(0.);
  u.n[right] = neumann(0.);
  u.n[top] = neumann(0.);
  pstn.refine = pstn.prolongation = fraction_refine;

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
  fraction (f, - y); //set the water depth _h
  fraction (pstn, PISTON); //set the piston fraction
  while (adapt_wavelet_leave_interface({u.x, u.y, p},{f, pstn},(double[]){uemax,uemax,femax,pemax, femax}, max_LEVEL+EXTRA_PISTON_LEVEL, LEVEL,padding, (int[]){max_LEVEL, max_LEVEL+EXTRA_PISTON_LEVEL}).nf){  //for adapting more around the piston interface
    fraction (f, - y); //set the water level on the refined mesh
    fraction (pstn, PISTON); //set the piston fraction on the refined mesh
  }
  unrefine ((x < -0.1)&&(level>6));
  foreach(){
    pf[] = 0;
    p[] = pf[];
  }
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field. Extra grid refinement around 
the piston is possible.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y, p},{f, pstn},(double[]){uemax, uemax, pemax, femax, pemax}, max_LEVEL+EXTRA_PISTON_LEVEL, LEVEL,padding, (int[]){max_LEVEL, max_LEVEL+EXTRA_PISTON_LEVEL});
  unrefine ((x>(Piston_Pos_x(x,y,z,t)+0.05))&&(level>=max_LEVEL));
  unrefine ((x < Piston_Pos_x(x,y,z,t)-piston_w*0.6)); //unrefine the area to the left of the piston
  unrefine ((y<-0.4)&&(x>(Piston_Pos_x(x,y,z,t)+0.02))); //unrefine the bottom
  unrefine ((y>0.1));
  fraction(pstn, PISTON);
}


/**
The moving piston is implemented via Stephane's trick. Which results in a leaking piston.
*/
event piston (i++, first) {
  fraction (pstn, PISTON);
  foreach(){
    u.y[] = u.y[]*(1 - pstn[]);
    u.x[] = pstn[]*Piston_Velo_x(x, y ,z, t) + u.x[]*(1 - pstn[]);
  }
  u.n[left]=dirichlet(Piston_Velo_x(x,y,z,t));
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
event vtu(t+=1, last){
  printf("Saving vtu file\n");
  char filename[100];
  sprintf(filename, "%s/TIME-%05.0f", vtu_folder, (t*100));
  output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
}

event save_energy(t+=0.01)
{
  printf("saving energy\n");
  char filename[200];
  sprintf(filename, "%s/energy.csv", results_folder);
  static FILE * fp = fopen(filename, "w");
  double ke = 0., gpe = 0., volume=0.;
  foreach (reduction(+:ke) reduction(+:gpe) reduction(+:volume)) {
    double norm2 = 0.;
    if (x>Piston_Pos_x(x,y,z,t)){
      foreach_dimension()
        norm2 += sq(u.x[]);
      ke += norm2*f[]*dv();
      gpe += (y+_h/2)*f[]*dv();
      volume += f[]*dv();
    }
  }
  if (i == 0)
    fprintf (fp, "t, ke, gpe, f\n");
  fprintf(fp, "%f, %f, %f, %f\n", t, rho1*ke/2., rho1*g_*gpe, volume);
}

event stop (t = Tend);

event show_progress(i++)
{
  printf("t=%02.3f, i=%04d, dt=%.3g, run:%d, LEVEL=%d, Extra_piston_level=%d\n", t, i, dt,run_number, max_LEVEL, EXTRA_PISTON_LEVEL);
}


event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
