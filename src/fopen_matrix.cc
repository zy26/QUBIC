#include "fopen_matrix.h"

int FopenMatrix::internal::fscanf2(FILE *file, short *x)
{
  return fscanf(file, "%hd", x);
}

int FopenMatrix::internal::fscanf2(FILE *file, double *x)
{
  return fscanf(file, "%lf", x);
}

int FopenMatrix::internal::fscanf2(FILE *file, float *x)
{
  return fscanf(file, "%f", x);
}

int FopenMatrix::internal::fscanf2(FILE *file, int *x)
{
  return fscanf(file, "%d", x);
}
