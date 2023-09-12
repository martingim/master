//based on http://basilisk.fr/sandbox/M1EMN/Exemples/dambNS.c

#include "navier-stokes/centered.h"
#include "vof.h"
#define RHOF 1.e-3
#define MU 1./1000
#define LEVEL 8

scalar f[], *interfaces = {f};
face vector alphav[];
face vector muv[];

int main() {
    L0 = 10.;
    DT = 1e-2;
    init_grid (1<< LEVEL);
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = neumann(0);

    run();
}

event init (t=0) {
    mask (y > L0/4 ? top :
        none);
    const face vector g[] = {0,-1};
    a = g;
    alpha = alphav;
    scalar phi[];
    foreach_vertex(){
        phi[] = min(1-y,L0/2-x);
    }
    fractions (phi, f);

    foreach() {
        u.x[] = f[] * (1e-8);
        u.y[] = 0;
    }
}

#define rho(f) ((f) + RHOF*(1. - (f)))
#define muc(f) ((f)*MU + MU/100.*(1. - (f)))

event properties (i++) {
   mu = muv;
  	
    foreach_face() {
        double fm = (f[] + f[-1])/2.;
        alphav.x[] = 1./rho(fm);
        muv.x[] = muc(fm);
    }
    boundary ((scalar *){muv,alphav});
}

void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}

event logfile (i++) {
    stats s = statsf (f);
    fprintf (stderr, "%g %d dt=%g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print (mgp);
    mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}

event interface (t += 1;t <=3){
    output_facets (f);
    fprintf(stdout,"\n");
    fprintf(stderr, "------~~~~~~-----  \n" );
}
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = -1 + f[]* (.25+    sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l}); 
    output_ppm (l, min = -1,  max= 3.25,
      linear = true,
                n = 1024, box = {{0,-1},{10,2.5}},
                file = "dambNS.mp4"
                );                   
    foreach()
    l[] = level;
    boundary ({l}); 
    output_ppm (l, n = 1024,  min = 0.,   max= 12.,
      box = {{0,-1},{10,2.5}},
      file = "level.mp4")  ;         
}
#if 0
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = -1 + f[]* (.25+    sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    static FILE * fp2 = popen ("ppm2mpeg > dambNS.mpg", "w");
    output_ppm (l, fp2 , min = -1,  max= 3.25,
      linear = true,
                n = 2048, box = {{0,-1},{10,2.5}}
                );          
    
    foreach()
    l[] = level;
    boundary ({l});
    static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
    output_ppm (l, fp1,   box = {{0,0},{10,2.5}})  ;          
}
#endif
event pictures (t=0.05) {
    scalar l[];
    foreach()
    l[] = f[]* (   sqrt(sq(u.x[]) + sq(u.y[])));; 
    boundary ({l});
     output_ppm (l, file = "dambNS.png", min = 0,   max= 2,
     // linear = true,
                //n = 1024 ,
                 box = {{0,0},{10,2.5}}
                );
    foreach()
    l[] = level; 
    boundary ({l});
     output_ppm (l, file = "level.png", min = 0,   max= 12,
     // linear = true,
                //n = 1024 ,
                 box = {{0,0},{10,2.5}}
                );
             
}

#if QUADTREE
event adapt(i++){
 scalar g[];
 foreach()
   g[]=f[]*(1+noise())/2;
 boundary({g});
 //adapt_wavelet ({g,f,u.x,u.y}, (double[]){0.001,.01,0.01,0.01}, minlevel = 5,maxlevel = LEVEL);
 adapt_wavelet({g,f},(double[]){0.001,0.01},minlevel = 5,maxlevel = LEVEL);
} 
#endif