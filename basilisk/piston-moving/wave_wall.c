/*
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 * reads piston position data from file.
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

#include "output_vtu_foreach.h"

int set_n_threads = 0; //0 use all available threads
int LEVEL = 6;
int MAXLEVEL = 12;
int padding = 2;
#define _h 0.6//water depth
double l = 25.6; //the size of the domain, preferable if l=(water_depth*2**LEVEL)/n where n is an integer
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.1;  //TODO change these to be based on LEVEL
double uemax = 0.01;
double pemax = 15.;
double Tend = 30.;
double probe_positions[]={8.00, 10.04, 10.75, 11.50};
int n_probes = 4;
vector h[]; //scalar field of the distance from the surface, using heights.h
char save_location[] = "./"; //the location to save the vtu files

//piston file 
char piston_file[] = "fil3.dat";
int file_samplerate = 100; //the samplerate of the piston position file
#define piston_timesteps 10000//the number of timesteps in the piston file
int piston_counter;
double piston_positions[piston_timesteps];
double piston_position = 0; //the starting position of the piston
double piston_position_p; //piston position at the previous timestep
double U_X = 0.; //the speed of the piston
//piston parameters 
#define piston_back_wall_offset 1. //the distance from the left of the domain to the front of the piston
double piston_height = 0.2;   //height of the piston above the still water level
double piston_bottom_clearance = 0.00; //piston distance above bottom
#define piston_w .2 //width of the piston
#define PISTON (1 -(piston_position<x) - (x<(piston_position-piston_w))-(y>piston_height)) //subtract conditions outside the piston
scalar pstn[];

void mask_domain(){
  //mask away the top of the domain
  mask(y > domain_height - _h ? top : none);
  //add masking the horizontal length
}

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
    piston_positions[count] /=-100.; //convert to meters
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

int main() {
  mkdir("./vtu",0755);
  read_piston_data();
  init_grid (1 << (LEVEL));
  mask_domain();
  origin(-piston_back_wall_offset, -_h);
  

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

  u.n[left] = neumann(0.);
  u.n[right] = dirichlet(0.);

  u.n[top] = neumann(0.);
  p[top]    = dirichlet(0.);
  pf[top]   = dirichlet(0.);
  

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
  
  fraction (f, - y); //set the water depth _h
  fraction (pstn, PISTON); //set the piston fraction
  int keep_refining=1;
  while (keep_refining){
    keep_refining = adapt_wavelet_leave_interface({u.x, u.y},{pstn,p,f},(double[]){uemax,uemax,femax,femax, 1.}, MAXLEVEL, LEVEL,padding).nf;
    fraction (f, - y); //set the water level on the refined mesh
    fraction (pstn, PISTON); //set the piston fraction on the refined mesh
  }
  unrefine ((x < -0.1)&&(level>6));
  foreach(){
    pf[] = 0;
    p[] = pf[];
  }
  fraction (f, - y); //set the water level on the refined mesh
  fraction (pstn, PISTON); //set the piston fraction on the refined mesh
}

/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
*/
event piston (i++, first) {
  piston_counter = floor(t*file_samplerate);
  double counter_remainder = 0;
  counter_remainder = t*100.-piston_counter;
  piston_position = piston_positions[piston_counter] + (piston_positions[piston_counter+1] -piston_positions[piston_counter])*counter_remainder; //update the piston position
  //printf("t:%f, file_timestep:%d, %%to next file timestep:%.0f%%, piston_position:%f\n", t, piston_counter, counter_remainder*100, piston_position);
  U_X = (piston_position-piston_position_p)/dt;
  piston_position_p = piston_position;
  fraction (pstn, PISTON);
  foreach() {
    u.y[] = u.y[]*(1 - pstn[]);
    u.x[] = pstn[]*U_X + u.x[]*(1 - pstn[]);
  }
  //printf("U_X:%f, piston_position:%f\n", U_X, piston_position);
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y},{pstn,p,f},(double[]){uemax, uemax, femax, pemax, femax}, MAXLEVEL, LEVEL,padding);
  unrefine ((x < piston_position-piston_w*0.6)); //unrefine the area to the left of the piston
  //unrefine ((x > piston_position+0.1)&&(y<-0.4)); //unrefine the bottom
  unrefine ((x>piston_position)&&(f[]<0.01)); //unrefine the air
}

event surface_probes(t+=0.01){
  char filename[40];
  sprintf(filename, "surface_probes.csv");
  FILE *fp = fopen(filename, "a");
  fprintf(fp, "%f, ", t);
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
    fprintf(fp, "%f, ", surface_elevation);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

//save unordered mesh
event vtu(t+=1){
  printf("Saving vtu file\n");
  char filename[40];
  sprintf(filename, "%svtu/TIME-%04.0f", save_location, (t*100));
  output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
}


// event vtu(i++){
//   printf("Saving vtu file\n");
//   char filename[40];
//   sprintf(filename, "%svtu/step-%04d", save_location, i);
//   output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
// }


/*
simulation stopped at Tend
*/
event stop (t = Tend);

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%02.3f, i=%04d, dt=%.3g\r", t, i, dt);
  //printf("%.2f%%\r", progress*100);
}
