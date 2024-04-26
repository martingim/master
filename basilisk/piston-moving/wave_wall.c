/*
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
//#include "reduced.h"

#include "output_vtu_foreach.h"

int LEVEL = 6;
int MAXLEVEL = 12;
int padding = 0;
double tank_length = 24.6; //length of the wave tank from the piston
double l = 32; //the size of the domain later masked to match the length of the wave tank from the piston
double domain_height = 1.0; //the height of the simulation domain
double femax = 0.1;
double uemax = 0.0001;
double Tend = 30.;
double height_probes[]={8.00, 10.04, 10.75, 11.50};
char save_location[] = "./"; //the location to save the vtu files
#define file_input 1 //whether the piston positions are read from file or given by function
char piston_file[] = "fil3.dat";

#define _h 0.6//water depth

//piston parameters
#define piston_starting_position 0.4
double piston_position = piston_starting_position; //the starting position of the piston (at rest leftmost position)
double piston_height = 0.2;   //height of the piston above the still water level
double piston_bottom_clearance = 0.0; //piston distance above bottom
#define w .1 //width of the piston
#define PISTON (0.1-(x>piston_position)-(x<piston_position-w)-(y<-_h+piston_bottom_clearance)-(y>piston_height)) //subtract conditions outside the piston
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
#define Piston_Position (piston_amplitude*(1-cos(t*omega)))
#define U_X (piston_amplitude*omega*sin(t*omega))  //piston velocity from pos -cos(t*omega)
#endif


void mask_domain(){
  //mask away the top of the domain
  mask(y > domain_height - _h ? top : none);
  mask(x > tank_length ? right : none);
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
    piston_positions[count] /=100.; //convert to meters
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
#endif

int main() {
  origin(-piston_starting_position, -_h);
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
  N = 1 << LEVEL;
  DT = 0.01;
  u.t[bottom] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
#if _OPENMP
  int num_omp = omp_get_max_threads();
  fprintf(stderr, "number of openmp threads:%d\n", num_omp);
#endif
  run();
  
}


event init (i = 0) {
  init_grid (1 << (LEVEL));
  mask_domain();
  fraction (f, - y); //set the water depth _h 
  fraction (pstn, PISTON); //set the piston fraction
  int keep_refining=1;
  while (keep_refining){
    keep_refining = adapt_wavelet_leave_interface({u.x, u.y},{pstn,p,f},(double[]){uemax,uemax,femax,femax, 1.}, MAXLEVEL, LEVEL,padding).nf;
    foreach(){
      pf[] = f[]*(-y)/10./1e-5;
      p[] = pf[];
    }
    fraction (f, - y); //set the water level on the refined mesh
    fraction (pstn, PISTON); //set the piston fraction on the refined mesh
  }
  unrefine ((x < -0.05)&&(level>6));
  foreach(){
    pf[] = f[]*(-y)/10./1e-5;
    p[] = pf[];
  }
  fraction (f, - y); //set the water level on the refined mesh
  fraction (pstn, PISTON); //set the piston fraction on the refined mesh
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
  foreach() {
    u.y[] *= (1 - pstn[]);
    u.x[] = u.x[]*(1 - pstn[]) + pstn[]*U_X;
  }
  printf("U_X:%f, piston_position:%f\n", U_X, piston_position);
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -9.81;
}

/*
The grid is adapted to keep max refinement at the air water interface.
And to minimise the error in the velocity field.
 */
event adapt (i++){
  adapt_wavelet_leave_interface({u.x, u.y},{pstn,p,f},(double[]){uemax, uemax, femax, femax,1.0}, MAXLEVEL, LEVEL,padding);
  unrefine ((x < -0.05)&&(level>6)); //unrefine the area to the left of the piston
  unrefine (y>0.05);
}

//save unordered mesh
event vtu(t+=.1){
//event vtu(i++){
  printf("Saving vtu file\n");
  char filename[40];
  sprintf(filename, "%svtu/TIME-%04.0f", save_location, (t*100));
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
/*
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
/*
simulation stopped at Tend
*/
event stop (t = Tend);

event show_progress(i++)
{
  float progress = 0;
  progress = t /Tend;
  printf("t=%.3f, i=%d, dt=%g, ", t, i, dt);
  printf("%.2f%%\n", progress*100);
}
