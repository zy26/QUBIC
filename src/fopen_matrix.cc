#include "fopen_matrix.h"

void FopenMatrix::internal::fscanf2(FILE *file, short *x)
{
  fscanf(file, "%hd", x);
}

void FopenMatrix::internal::fscanf2(FILE *file, double *x)
{
  fscanf(file, "%lf", x);
}

void FopenMatrix::internal::fscanf2(FILE *file, float *x)
{
  fscanf(file, "%f", x);
}

void FopenMatrix::internal::fscanf2(FILE *file, int *x)
{
  fscanf(file, "%d", x);
}
