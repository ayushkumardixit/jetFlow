#include "axi.h"
#include "navier-stokes/centered.h"

char filename[80];
int n = 40;
double array[2][40] = {0};

int main(int a, char const *arguments[])
{
    sprintf(filename, "%s", arguments[1]);
    restore(file = filename);
    boundary((scalar *){u.x, u.y});
    for (int i = 0; i < n; i++)
    {
        int nocells = 0;
        array[0][i] = i + 0.5;
        foreach ()
        {
            if (x < (i + 1) && x > i)
            {
                array[1][i] +=  2 * pi * y * Delta * u.x[];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        fprintf(ferr, "%d %g %g\n", nocells, array[0][i], array[1][i]);
    }
}