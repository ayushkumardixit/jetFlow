#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char nameOut[80], dumpFile[80];

const double radius = 25e-3;
const double Q = 0.125;
double U = Q / (pi * radius * radius);
double length = 60e-3;
int MAXlevel = 9;
double Ldomain = 0.5;

scalar f0[], f[];
face vector muv[];
scalar rhov[];

#define tsnap (1e-2) // 0.001 only for some cases.
#define tmax (1e2)
// Error tolerancs
#define VelErr (1e-2) 

// Boundary conditions

u.n[left] = dirichlet(f0[] * U);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

// u.t[top] = neumann(0.);
// u.n[top] = neumann(0.);
// p[top] = dirichlet(0.);

int main(int argc, char const *argv[])
{
    L0 = Ldomain;
    init_grid(1 << 7);
    char comm[80];
    sprintf(comm, "mkdir -p intermediate");
    system(comm);
    // Name of the restart file. See writingFiles event.
    sprintf(dumpFile, "dump");

    rho = rhov;
    mu = muv;
    run();
}
event properties(i++)
{
    foreach ()
        rhov[] = 1e-3;
    foreach_face()
        muv.x[] = 1.81e-5;
}
event init(t = 0)
{
    if (!restore(file = dumpFile))
    {
        refine(y < 5. * radius && level < MAXlevel);
        fraction(f0, radius - y);
        f0.refine = f0.prolongation = fraction_refine;
        restriction({f0}); // for boundary conditions on levels
        foreach ()
        {
            f[] = f0[] * (x < length);
            u.x[] = U * f[];
        }
    }
}
// event adapt(i++)
// {
//     adapt_wavelet((scalar *){u.x, u.y}, (double[]){VelErr, VelErr}, MAXlevel, MAXlevel - 3);
// }
event writingFiles(t = 0; t += tsnap; t <= tmax)
{
    dump(file = dumpFile);
    sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
    dump(file = nameOut);
}
event logWriting(i++)
{
    double ke = 0.;
    foreach (reduction(+ : ke))
    {
        ke += (2 * pi * y) * (0.5 * rho[] * (sq(u.x[]) + sq(u.y[]))) * sq(Delta);
    }
    FILE * fp = fopen ("log", "a");
    fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
}
