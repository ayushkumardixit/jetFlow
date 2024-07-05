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
        for (int j = 0; j < n; j++)
        {
            array[1][i] += 2 * pi * (j + 0.5) * interpolate(u.x, i + 0.5, j + 0.5);
        }
        array[0][i] = i + 0.5;

        fprintf(ferr, "%.2e %.2e\n", array[0][i], array[1][i]);
    }
}
