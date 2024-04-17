/**
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 * can read piston position data from file. Set up for file with 100hz sample rate 
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "adapt_wavelet_leave_interface.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

#include "output_vtu_foreach.h"

int LEVEL = 7;
int MAXLEVEL = 10;
double l = 4; //the length of the wave tank
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.1;
double uemax = 0.1;
double Tend = 30.;

char save_location[] = "./"; //the location to save the vtu files
#define file_input 1 //whether the piston positions are read from file or given by function
char piston_file[] = "fil3.dat";

#define _h 0.6//water depth

//piston parameters
#define piston_starting_position 0.4
double piston_position = piston_starting_position; //the starting position of the piston (at rest leftmost position)
double h = 1.9;   //height of the piston 
double hb = 0.0; //piston height above bottom
#define w .1 //width of the piston
#define PISTON (w - (fabs(piston_position-x)) - (y > h) - (y < hb))
scalar pstn[];

//differences based on reading piston data from file or from function
#if file_input    //reading data from file
int piston_counter;
#define piston_timesteps 10000
double piston_positions[piston_timesteps];
double piston_time[piston_timesteps];
double piston_position_p; //previous piston position
double U_X = 0.; //the speed of the piston
#else   //piston data from function
double omega = 8.95;
double piston_amplitude = 0.0129;
#define Piston_Position (piston_starting_position+piston_amplitude*(1-cos(t*omega)))
#define U_X (piston_amplitude*omega*sin(t*omega))  //piston velocity from pos -cos(t*omega)
#endif




void mask_domain(){
  //mask away the top of the domain
  mask(y > domain_height ? top : none);
}

#if file_input
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
    piston_positions[count] /=100.;
    piston_positions[count] += piston_starting_position;
    count++;
  }
  fclose(file);
  piston_position_p = piston_positions[0];
  piston_position = piston_positions[0];
}
#endif

int main() {
  #if file_input
  read_piston_data();
  #endif
  mkdir("./vtu",0755);
  L0 = l;
  f.sigma = 0.01;
  rho1 = 1025;
  rho2 = 1.225;
  mu1 = 8.9e-4; 
  mu2 = 17.4e-6;
  G.y = -9.81;
  //Z.y = _h; //set the reduced gravity reference level to the water level at rest
  N = 1 << LEVEL;
  DT = 0.01;
  //u.n[top]= neumann(0.);
  //p[left] = neumann(0.);
  u.n[left]  = neumann (0.);
  u.t[bottom] = dirichlet(0.);
  run();
}


event init (i = 0) {
  init_grid (1 << (MAXLEVEL));
  //pstn.prolongation = pstn.refine = fraction_refine;
  fraction (f, _h - y); //Water depth _h 
  fraction (pstn, PISTON);
  //adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax,uemax,femax,1.0}, MAXLEVEL, LEVEL,1);
  //fraction (f, _h - y); //Water depth _h 
  //fraction (pstn, PISTON);
  
  //while (adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax,uemax,femax,1.0}, MAXLEVEL, LEVEL,1).nc){
    //fraction (f, _h - y); //Water depth _h
    //fraction (pstn, PISTON);
  //}
  /* 
  char filename2[40];
  sprintf(filename2, "%svtu/adapt-%06g", save_location, (t*1000));
  output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename2);
 */
  mask_domain();
}


/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
 */
event piston (i++) {
#if file_input
  piston_counter = floor(t*100);
  double counter_remainder = 0;
  counter_remainder = t*100.-piston_counter;
  piston_position = piston_positions[piston_counter] + (piston_positions[piston_counter+1] -piston_positions[piston_counter])*counter_remainder; //update the piston position
  printf("t:%f, file_timestep:%d, %%to next file timestep:%.0f%%, piston_position:%f\n", t, piston_counter, counter_remainder*100, piston_position);
  U_X = (piston_position-piston_position_p)/dt;
  piston_position_p = piston_position;
#else
  piston_position = Piston_Position;
#endif

  fraction (pstn, PISTON);
  //u.n[left]  = dirichlet (-U_X);
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
    //p[] *= (1 - pstn[]);
  }
  printf("U_X:%f, piston_position:%f\n", U_X, piston_position);
  
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax, uemax, femax, 1.0}, MAXLEVEL, LEVEL,0);
}

//save unordered mesh
event vtu(t+=.1){
  char filename[40];
  sprintf(filename, "%svtu/TIME-%d", save_location, (i));
  output_vtu((scalar *) {f,p,pstn}, (vector *) {u}, filename);
}

/**
## Movie

   For systems using Bview, `movie.mp4` may be generated.

![The water, piston and grid structure (mirrored)](wave_wall/movie.mp4)

We can see vortex shedding from the piston, wave overtopping and no
significant jetting at the walls. Meaning that the setup could be much
improved.
 */

#if 0
#include "view.h"
event bview_movie (t += 0.1) {
  view (tx = -0.5);
  draw_vof ("pstn", filled = 1, fc = {0.2,0.2,0.2});
  draw_vof ("f", filled = 1, fc = {0.1,0.1,0.9});
  begin_mirror ({0,-1});
  cells();
  end_mirror();
  save ("movie.mp4");
}
/**
Else, `f.mp4` will have to reaveal the dynamics.
 */
#else
/*
event movie (t += 0.1) {
  foreach() {
    if (pstn[] > 0.5)
      pstn[] = -1;
  }
  output_ppm (f, file = "f.mp4", mask = pstn,
	      n = 512, box = {{0,0},{l,domain_height}});
}
*/
#endif
/**
At t = 10, the simulation is stopped.
 */
event stop (t = Tend);

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%f, ", t, i, dt);
  printf("%.2f%%\n", progress*100);
}
