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
int MAXLEVEL = 11;
double l = 4; //the length of the wave tank
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.1;
double uemax = 0.01;
double Tend = 30.;
double piston_amplitude = 0.0129;
char save_location[] = "./";
double omega = 8.95;

int piston_counter = 0;
#define piston_timesteps 10000
double piston_positions[piston_timesteps];


double xp = 0.1; //the starting position of the piston (at rest leftmost position)
double h = .9;   //Max. height
double hb = 0.0; //Height above bottom floor

#define w .1 //width of the piston
#define _h 0.6//water depth
//#define XP (0.1+piston_amplitude-piston_amplitude*cos(t*omega))
//#define U_X (piston_amplitude*omega*sin(t*omega))  //piston velocity from pos -cos(t*omega)
#define PISTON (w - (x > xp) - (y > h) - (y < hb))
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

  while(!feof(file) && count< piston_timesteps ){
    fscanf(file, "%lf", &(piston_positions[count++]));
  }
  fclose(file);
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
  run();
}
/**
initialization
 */
event init (i = 0) {
  init_grid (1 << (LEVEL));
  refine ( abs(y-_h)<0.05 && level < MAXLEVEL);
  //pstn.prolongation = pstn.refine = fraction_refine;
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
event piston (t+=0.01; t < Tend) {
  piston_counter += 1;
  xp = piston_positions[piston_counter]; //update the piston position based on the piston speed
  //xp = XP;      //update the piston position based on the given formula
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
  printf("U_X:%f, xp:%f\n", U_X, xp);
  
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax, uemax, femax, 1.0}, MAXLEVEL, LEVEL,0);
}

//save unordered mesh
event vtu(t+=1;t<Tend){
  char filename[40];
  sprintf(filename, "%svtu/TIME-%06g", save_location, (t*100));
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
event movie (t += 0.1) {
  foreach() {
    if (pstn[] > 0.5)
      pstn[] = -1;
  }
  output_ppm (f, file = "f.mp4", mask = pstn,
	      n = 512, box = {{0,0},{l,domain_height}});
}
#endif
/**
At t = 10, the simulation is stopped.
 */
event stop (t = Tend);

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%f", t, i, dt);
  printf("%.2f%%\n", progress*100);
}
