// System includes
#include <stdio.h>
#include <assert.h>

// CUDA runtime
#include <cuda_runtime.h>

// helper functions and utilities to work with CUDA
#include <helper_functions.h>
#include <helper_cuda.h>
#include <cuda.h>

	
#include "cublas.h"
	
int main(int argc, char **argv)
{
// 	int devID;
// 	cudaDeviceProp props;
// 
// 	// This will pick the best possible CUDA capable device
// 	devID = findCudaDevice(argc, (const char **)argv);
// 
// 	//Get GPU information
// 	checkCudaErrors(cudaGetDevice(&devID));
// 	checkCudaErrors(cudaGetDeviceProperties(&props, devID));
// 	printf("Device %d: \"%s\" with Compute %d.%d capability\n",
// 		devID, props.name, props.major, props.minor);
	cublasStatus status;
	status = cublasInit();
	if (status != CUBLAS_STATUS_SUCCESS)
	{
	 	printf("!!!! CUBLAS initialization error\n");
	 	return;
	}
	cublasShutdown();
	int imageSize = 512;
	unsigned int*p_devRgbU; 
	cudaMalloc((void**)&p_devRgbU,/* sizeof(unsigned int)**/imageSize);
	cudaFree(p_devRgbU);

	return 1;
}