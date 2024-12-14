#include "mdct.h"
#include <math.h>
#include <vector>

__global__ void easy_IMDCT_kernel (double *X, double *y, double *window, int n) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	while (idx < n) {
		y[idx] = 0.0;
		for(int m = 0; m < n / 2; ++m) {
			y[idx] += X[m] * cos(M_PI / 2 / n * (2 * idx + 1 + n / 2) * (2 * m + 1));
		}
		y[idx] *= window[idx];
		idx += stride;
	}
}
void IMDCT_gpu (double *X, double *y, double *window, int n) {
	static double *X_d, *y_d, *window_d;
	static int ok = 0;
	if (!ok) {
		cudaMalloc(&X_d, sizeof(double) * 8192);
		cudaMalloc(&y_d, sizeof(double) * 8192);
		cudaMalloc(&window_d, sizeof(double) * 8192);
		ok++;
	}
	cudaMemcpy(X_d, X, sizeof(double) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y, sizeof(double) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(window_d, window, sizeof(double) * n, cudaMemcpyHostToDevice);
	int sm_count = n / 256;
	easy_IMDCT_kernel<<<sm_count, 256>>>(X_d, y_d, window_d, n);
	cudaDeviceSynchronize();
	cudaMemcpy(y, y_d, sizeof(double) * n, cudaMemcpyDeviceToHost);
}
