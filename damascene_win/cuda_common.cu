#include "cuda_common.cuh"
#include <cutil.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <helper_cuda.h>

float* getImage(uint width, uint height, float* devImage) 
{
	int imageSize = width * height * sizeof(float);
	float* result = (float*)malloc(imageSize);
	CUDA_SAFE_CALL(cudaMemcpy(result, devImage, imageSize, cudaMemcpyDeviceToHost));
	return result;
}

int* getImage(uint width, uint height, int* devImage) 
{
	int imageSize = width * height * sizeof(int);
	int* result = (int*)malloc(imageSize);
	CUDA_SAFE_CALL(cudaMemcpy(result, devImage, imageSize, cudaMemcpyDeviceToHost));
	return result;
}

