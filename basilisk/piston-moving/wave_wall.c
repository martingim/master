/**
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 */
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "adapt_wavelet_leave_interface.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "output_vtu_foreach.h"

int LEVEL = 4;
int MAXLEVEL = 13;
double l = 16; //the length of the wave tank
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.1;
double uemax = 0.01;
double Tend = 30.;
double piston_amplitude = 0.0129;
char save_location[] = "./";
double omega = 8.95;

int piston_counter;
#define piston_timesteps 10000
double piston_positions[piston_timesteps];
double piston_time[piston_timesteps];

#define piston_starting_position 0.1
double piston_position = piston_starting_position; //the starting position of the piston (at rest leftmost position)


#define file_input 1
#if file_input
double piston_position_p=piston_starting_position; //previous piston position
double U_X = 0.; //the speed of the piston
#else
#define Piston_Position (0.1+piston_amplitude-piston_amplitude*cos(t*omega))
#define U_X (piston_amplitude*omega*sin(t*omega))  //piston velocity from pos -cos(t*omega)
#endif

double h = .9;   //Max. height
double hb = 0.0; //Height above bottom floor

#define w .1 //width of the piston
#define _h 0.6//water depth
#define PISTON (w - (x > piston_position) - (y > h) - (y < hb))
scalar pstn[];

//mask away the top of the domain
void mask_domain(){
  mask(y > domain_height ? top : none);
}

void read_piston_data(){
  int count = 0;
  FILE *file;
  file = fopen("fil3.dat", "r");
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

int main() {
  read_piston_data();
  mkdir("./vtu",0755);
  L0 = l;
  f.sigma = 0.01;
  rho1 = 1025;
  rho2 = 1.225;
  mu1 = 8.9e-4; 
  mu2 = 17.4e-6;
  N = 1 << LEVEL;
  DT = 0.01;
  run();
}
/**
initialization
 */
event init (i = 0) {
  init_grid (1 << (MAXLEVEL));
  //refine ( abs(y-_h)<0.05 && level < MAXLEVEL);
  printf("pisoton_pis:%f\n", piston_position);
  pstn.prolongation = pstn.refine = fraction_refine;
  fraction (f, _h - y); //Water depth _h
  //f.prolongation = f.refine = fraction_refine;
  //adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax,uemax,femax,1.0}, MAXLEVEL, LEVEL,0);
  mask_domain();
}
/**
The force of gravity is included with $g = 9.81$
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -9.81;
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
  printf("t:%f, piston_counter:%d, counter_remainder:%f, piston_position:%f, piston_position_p:%f, asdf:%f\n", t, piston_counter, counter_remainder, piston_position, piston_position_p, piston_positions[piston_counter+1] -piston_positions[piston_counter]);
  U_X = (piston_position-piston_position_p)/dt;
  piston_position_p = piston_position;
#else
  piston_position = Piston_Position;
#endif

  fraction (pstn, PISTON);
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
    /*//set constant waterlevel inside piston
    if (pstn[] > 0.5 && y < _h){
      f[] = 1.;
      }
    */
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
event vtu(t+=1 ;t<Tend){
  char filename[40];
  sprintf(filename, "%svtu/TIME-%06g", save_location, (t*1000));
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
