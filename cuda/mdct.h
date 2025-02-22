#ifndef MDCT_H
#define MDCT_H
extern "C" {
#include <stdlib.h>
__global__ void easy_IMDCT_kernel(double *X, double *y, double *window, int n);
void IMDCT_gpu(double *X, double *y, double *window, int n);
}
#endif
