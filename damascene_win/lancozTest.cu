#include "stdafx.h"


#define MAXITER 6000
#define CHECKITER 500
#define LUMPTOL 1e-5
#define TOLERANCE 1e-3
#define SPURTOLERANCE 1e-10
#define MAXEIGNUM 10



#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "getopt.h"
#include <time.h>
#include <math.h>
#include "cublas.h"
#include <assert.h>
#include "cutil.h"
#include <windows.h>
#include "csr.h"
#include "unistd.h"
#include "fcntl.h"
#include <vector>
#include <iostream>
#include <time.h>
#include<fstream>
#include "lanczos.h"
#include "stencilMVM.h"
#include "intervening.h"

using namespace std;


int main()
{
	  int width, height;
    int N;

    int fp;
    fp = open("E:\\ct-contour-reconstruction\\Release\\matrix.dat", O_RDONLY);


	read(fp, &width, sizeof(int)); 
	read(fp, &height, sizeof(int));

    N = width*height;

    vector<int> rowptr(N+1);
    vector<int> cols;
    vector<float> vals;
    
    vector<float> temp(81*sizeof(float)*N);
	read(fp, &temp[0], 81*sizeof(float)*N);
    close(fp);

	printf("\n¾ØÕóÎ¬ÊýÎª:%d  and  %d£¡",width,height);

    int xoffsets[81];
    int yoffsets[81];
    int p = 0;
    for(int y=-5;y<=5;y++)
    {
        for(int x=-5;x<=5;x++)
        {
            if(x*x+y*y <= 25)
            {
                xoffsets[p] = x;
                yoffsets[p] = y;
                p++;
            }
        }
    }

    rowptr[0] = 0;
    int elt=0;

    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            int id = y*width+x;
            for(int o=0;o<81;o++)
            {
               if(y+yoffsets[o] >= 0 && y+yoffsets[o]<height && x+xoffsets[o]>=0 && x+xoffsets[o]<width)
               {
                   cols.push_back( id + yoffsets[o]*width+xoffsets[o]);
                   vals.push_back(temp[ o*width*height + id]);
                   elt++;
               }
            }
            rowptr[id+1] = elt;
        }
    }


    /* Actual code starts here */
	int nMatrixDimension = width * height;

	int radius = 5;
    CSRMatrix mat(N, elt, &rowptr[0], &cols[0], &vals[0]);

	dim3 blockDim(XBLOCK, 1);
	dim3 gridDim((width * height - 1)/XBLOCK + 1, 1);

	int matrixPitchInFloats = findPitchInFloats(nMatrixDimension);
	Stencil myStencil(radius, width, height, matrixPitchInFloats);

	float* devMatrix;

	printf("Reading matrix from file...\n");
	float* hostMatrix = &vals[0];
	printf("Copying matrix to GPU\n");


	uint nDimension = myStencil.getStencilArea();


	cudaMalloc((void**)&devMatrix, nDimension * nMatrixDimension * sizeof(float));

	CUDA_SAFE_CALL(cudaMemcpy(devMatrix, hostMatrix, nMatrixDimension * nDimension * sizeof(float), cudaMemcpyHostToDevice));
	//intervene(myStencil, devMatrix, &devDstMatrix);

	float* eigenValues;
	float* devEigenVectors = 0;
	float fTolerance = 0.001;
	generalizedEigensolve(myStencil, devMatrix, matrixPitchInFloats, MAXEIGNUM, &eigenValues, &devEigenVectors, fTolerance);

	return 1;
}