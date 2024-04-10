/**
 * based on http://basilisk.fr/sandbox/Antoonvh/wave_wall.c
 */
#include <sys/stat.h>
#include "utils.h"
#include "adapt_wavelet_leave_interface__.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "output_vtu_foreach.h"

int LEVEL = 4;
int MAXLEVEL = 10;
double femax = 0.1;
double uemax = 0.01;
double Tend = 10.;
double piston_amplitude = 0.02;
char save_location[] = "./";
double omega = 8.95;

double xp = 0.1; //the starting position of the piston (at rest leftmost position)
double h = .9;   //Max. height
double hb = 0.0; //Height above bottom floor

#define w .1 //width of the piston
#define _h 0.6//water depth
//#define XP (-piston_amplitude*cos(t*omega))
#define U_X (piston_amplitude*omega*sin(t*omega))  //piston velocity from pos -cos(t*omega)
#define PISTON (w - (x > xp) - (y > h) - (y < hb))
scalar pstn[];

/**
Furthermore, domain size and fluid properties are also chosen adhoc.
 */

void mask_domain(){
  mask(y > 1 ? top : none);
}


int main() {
  mkdir("./vtu",0755);
  L0 = 8;
  f.sigma = 0.01;
  rho1 = 1025;
  rho2 = 1.225;
  mu1 = 8.9e-4; 
  mu2 = 17.4e-6;
  N = 1 << LEVEL;
  run();
}
/**
## initialization

Ideally, a proper wave should be initialized here to avoid the
overhead of the piston wave making process. Now virtually all effort
is targeted towards the piston and not towards the actual wave - basin
wall collision.
 */
event init (i = 0) {
  init_grid (1 << (MAXLEVEL));  
  pstn.prolongation = pstn.refine = fraction_refine;
  fraction (f, _h - y); //Water depth _h
  //while(adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax,uemax,femax,femax}, MAXLEVEL, LEVEL,1).nf);{
  //  fraction(f, _h-y);
  //  fraction(pstn, PISTON);
  //}
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
event piston (i++; t < Tend) {
  xp += dt*U_X; //update the piston position based on the piston speed
  //xp = XP;      //update the piston position based on the given formula
  fraction (pstn, PISTON);
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
  }
  printf("\nU_X:%f\n", U_X);
}
/**
The grid is adapted to accuracy represent the piston fraction field,
the water fraction fild and the velocity components. 
 */

event adapt (i++){
  printf("\nadapting mesh\n");
  adapt_wavelet_leave_interface({u.x, u.y},{pstn,f},(double[]){uemax, uemax, femax, 1.0}, MAXLEVEL, LEVEL,0);
}
//save unordered mesh
event vtu(t+=.01;t<Tend){
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
	      n = 512, box = {{0,0},{10,5}});
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
  printf("\r t=%.3f, i=%d, ", t, i);
  printf("%.2f%%\r", progress*100);
}
