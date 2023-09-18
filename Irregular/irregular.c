/**
# Irregular wave case using periodic domain (multilayer solver)

outputs .vts files which can be viewed in paraview
The solution is obtained using the layered model and demonstrates its
robustness and a degree of realism even for this complex case.


Made by: Oystein lande 2022  (modified version of basilisk/src/test/breaking.c)

*/

//#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "view.h"


#include "CFDwavemaker.h"

/**
The initial conditions is set in the waveinput.dat file and read through
the CFDwavemaker lib. This particular example is a spatial periodic
solution with a domain size L 1738m. We run the simulation for approximately
5 mean wave periods.*/

#define l_  1738.
#define k_  (2.*pi)/l_
#define h_  76.4
#define g_  9.81
#define T0  15.
#define Tmax 5*T0

/**
The domain is periodic in $x$ and resolved using 128$^2$
points and 20 layers. */

int main()
{
        size(l_);
        //omp_set_num_threads(1);
        origin(-L0 / 2., -L0 / 2.);
        periodic(right);
        periodic(top);
        N = 128;
        nl = 20;
        G = g_;
        //nu = 1. / RE;
        nu = 0;

        /* Initialize CFDwavemaker. waveinput.dat file is read when calling this function*/
        wave_Initialize();

        /**
        Some implicit damping is necessary to damp fast modes. This may be
        related to the slow/fast decoupling of the $\theta$-scheme mentioned
        by [Durran & Blossey, 2012](#durran2012). */

        //theta_H = 0.51;

        run();

        /*Release memory allocated to CFDwavemaker after simulation end*/
        wave_Cleanup();
}

/**
The initial conditions for the free-surface and velocity are given by
the third-order Stokes solution. */

event init(i = 0)
{

        /**
        We can use a larger CFL, in particular because we are not dealing
        with shallow-water/wetting/drying. */

        CFL = 0.8;

        /**
        The layer thicknesses follow a geometric progression, starting from
        a top layer with a thickness of 1/3 that of the regular
        distribution. */

        geometric_beta((1./3) * h_ / nl, true);


        // We set the seabed reference (zb), layer heights (h) and velocities (u.x u.y and w)
        foreach() {
                zb[] = -h_;
                double H = wave_SurfElev(x, y, 0) - zb[];
                double z = zb[];
                foreach_layer() {
                        h[] = H * beta[point.l];
                        z += h[] / 2.;
                        u.x[] = wave_VeloX(x, y, z, 0);
                        u.y[] = wave_VeloY(x, y, z, 0);
                        w[] = wave_VeloZ(x, y, z, 0);
                        z += h[] / 2.;
                }
        }
        boundary(all);
}

/**
We add (an approximation of) horizontal viscosity. */

event viscous_term(i++)
horizontal_diffusion((scalar*) {u}, nu, dt);

/**
We log the evolution of the kinetic and potential energies.

~~~gnuplot Evolution of the kinetic, potential and total energy
set xlabel 't/T0'
plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:3 w l t 'potential', \
        '' u 1:(($2+$3)/2.) w l t 'total/2'
~~~
*/

event logfile(i++; t <= Tmax)
{
        double ke = 0., gpe = 0.;
        foreach(reduction(+:ke) reduction(+:gpe)) {
                foreach_layer() {
                        double norm2 = sq(w[]);
                        foreach_dimension()
                                norm2 += sq(u.x[]);
                        ke += norm2 * h[] * dv();
                }
                gpe += sq(eta[]) * dv();
        }
        fprintf(stderr, "%g %g %g\n", t / T0, ke / 2., g_ * gpe / 2.);
}


/**
And generate the movie of the free surface, coloring with the horizontal
particle velocity u.x */


event movie(t += 0.1)
{
        view(fov = 17.3106, quat = { 0.475152,0.161235,0.235565,0.832313 },
                tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768);
        char s[80];
        sprintf(s, "t = %.2f T0", t / T0);
        draw_string(s, size = 80);
        //for (double x = -1; x <= 1; x++)
                //translate(x);
        squares("u.x", linear = true, z = "eta", min = -10., max = 12.);
        save("movie.mp4");
}