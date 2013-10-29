#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* bitcount : count 1 bits in x */
int bitcount_(double *x)
{
int b;
unsigned long *y;
double z;

z = *x;
y = (unsigned long *) &z;
for (b = 0; *y !=0; *y >>= 1)
  if (*y & 01)
    b++;
return b;
}

/* wrapper for IBM system */
int bitcount(double *x) 
{
  return bitcount_(x);
}
