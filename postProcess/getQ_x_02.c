#include "axi.h"
#include "navier-stokes/centered.h"

char filename[80];
const int n = 80;
double delta_get = 0.5 * 40 / n;
double array[2][n] = {0};

int main(int a, char const *arguments[])
{
    sprintf(filename, "%s", arguments[1]);
    restore(file = filename);
    boundary((scalar *){u.x, u.y});
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            array[1][i] += 2 * pi * (j + delta_get) * interpolate(u.x, i + delta_get, j + delta_get);
        }
        array[0][i] = i + delta_get;

        fprintf(ferr, "%.2e %.2e\n", array[0][i], array[1][i]);
    }
}
