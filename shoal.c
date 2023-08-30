#include "grid/multigrid.h"
#include "output_vtu_foreach.h"

#if ML
# include "layered/hydro.h"
# include "layered/nh.h"
# include "layered/remap.h"
# include "layered/perfs.h"
#else //!ML
# include "green-naghdi.h"
#endif

int main()
{
    X0 = -10.;
    L0 = 25.;
    Y0 = -L0/2.;
    G = 9.81;

#if ML
    N = 512;
    nl = 2;
    CFL_H = 1;
#else // Green-Naghdi
    N = 1024;
    gradient = NULL;

#endif
    run();
}

scalar maxa[];
event init(i = 0)
{
#if ML
    u.n[left]  = - radiation (0.06*sin(2.*pi*t/1.)); // 0.049
#else
  u.n[left]  = - radiation (0.042*sin(2.*pi*t/1.));
#endif
  u.n[right] = + radiation (0);
  double h0 = 0.45;
  double cosa = cos (20.*pi/180.), sina = sin(20.*pi/180);
  foreach() {
      double xr = x*cosa - y*sina, yr = x*sina + y*cosa;
      double z0 = xr >= -5.82 ? (5.82 + xr)/50. :0.;
      double zs = sq(xr/3.) + sq(yr/4.) <= 1. ?
      -0.3 + 0.5*sqrt(1. - sq(xr/3.75) - sq(yr/5.)) : 0.;
      zb[] = z0 + zs -h0;
#if ML
    foreach_layer()
        h[] = max(-zb[], 0.)*beta[point.l];
#else
    h[] = max(-zb[], 0.);
#endif  
    maxa[] = 0.;
  }
}

event friction (i++) {
    foreach()
        if (x > 12.) {
            double a = h[] < dry ? HUGE : 1. + 2.*(x - 12.)*dt*norm(u)/h[];
#if ML
        foreach_layer()
#endif 
            foreach_dimension()
                u.x[] /= a;        
    }
}

#if 1
event movie(t+=1) {
    char name[350];
    int t_;
    t_ = t;
    sprintf (name, "out_vtu%d", t_);
    output_ppm (eta, min = -0.04, max = 0.04, n=512);
    output_vtu ((scalar *) {h, eta}, (vector *) {u.x, u.y}, name);
}
#endif
event maximum (t = 40; i++){
    foreach()
        if (fabs(eta[])> maxa[])
            maxa[] = fabs(eta[]);
}

event end (t = 50) {
    FILE *fp = fopen("end", "w");
    output_field ({eta,zb,maxa}, fp, linear=true);
    FILE * fp2 = fopen ("section2", "w");
    FILE * fp5 = fopen ("section5", "w");
    FILE * fp7 = fopen ("section7", "w");
    for (double y = -10.; y<=10.; y += 0.02) {
        fprintf (fp2, "%g %g\n", y, interpolate(maxa, 3., y));
        fprintf (stderr, "%g %g\n", y, interpolate (maxa, 5., y));
        fprintf (fp5, "%g %g\n", y, interpolate (maxa, 9., y));
    }
    for (double x = -10.; x<=10.; x+=0.02)
        fprintf (fp7, "%g %g\n", x, interpolate (maxa, x, 0.));
}