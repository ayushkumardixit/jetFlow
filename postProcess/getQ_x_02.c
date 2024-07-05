#include "axi.h"
#include "navier-stokes/centered.h"

#define n 40

char filename[80];
//const int n = 80;
double delta_get = 0.5 * 40 / n;
double array[2][n] = {0};
double x_get, y_get;

int main(int a, char const *arguments[])
{
    sprintf(filename, "%s", arguments[1]);
    restore(file = filename);
    boundary((scalar *){u.x, u.y});
    for (int i = 0; i < n; i++)
    {
	x_get = 40*i/n + delta_get;
	array[0][i] = x_get;	
        for (int j = 0; j < n; j++)
        {
	    y_get = 40*j/n + delta_get;	
            array[1][i] += 2 * pi * y_get * delta_get * interpolate(u.x, x_get, y_get);
        }	
    	fprintf(ferr, "%d %.3e %.2e\n", i, array[0][i], array[1][i]);
    }
}
