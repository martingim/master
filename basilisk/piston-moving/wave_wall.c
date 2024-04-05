/**
 * 
 */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "output_vtu_foreach.h"

int LEVEL = 10;
int MAXLEVEL = 12;
double femax = 0.001;
double uemax = 0.01;
double Tend = 10.;
double piston_amplitude = 0.1;
char save_location[] = "./";
double frequency = 1.4244;
/**
No-slip side walls and set and the moving piston parameters are chosen
adhoc.
*/
//u.t[right] = dirichlet (0.);
//u.t[left] = dirichlet (0.);

double U_max = 1.5, //Maximum velocity
  tau = 0.25,       //Time scale for moving
  xp = 1.,         //Start location
  w = 0.25,         //Width
  h = .9,           //Max. height
  hb = 0.0;           //Height above bottom floor

#define _h 0.6//water depth
#define U_X (piston_amplitude*sin(t*frequency*3.14*2))
#define PISTON (w - (x > xp) - (y > h) - (y < hb))
scalar pstn[];

/**
Furthermore, domain size and fluid properties are also chosen adhoc.
 */


void mask_domain(){
  mask(y > 2.5 ? top : none);
}


int main() {
  L0 = 10;
  mu1 = 0.01;
  mu2 = 0.04;
  f.sigma = 0.01;
  rho1 = 1000.;
  rho2 = 1.;
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
  
  pstn.prolongation = pstn.refine = fraction_refine;
  fraction (f, _h - y); //Water depth $h = 2$
 
  while (adapt_wavelet((scalar *){pstn, f, u}, (double[]){femax, femax, uemax, uemax} , MAXLEVEL).nf){
    fraction (f, _h - y);
    adapt_wavelet((scalar *){pstn, f, u}, (double[]){femax, femax, uemax, uemax} , MAXLEVEL);
  }
   mask_domain();
}
/**
The force of gravity is included with $g = 10$
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -10;
}
/**
The moving piston is implemented via Stephane's trick. Note that this
piston is leaky.
 */
event piston (i++; t < 10) {
  xp += dt*U_X;
  fraction (pstn, PISTON);
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
  }
}
/**
The grid is adapted to accuracy represent the piston fraction field,
the water fraction fild and the velocity components. 
 */

event adapt (i++)
  adapt_wavelet ((scalar *){pstn, f, u}, (double[]){femax, femax, uemax, uemax} , MAXLEVEL);

//save unordered mesh
event vtu(t+=.1;t<10){
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
