// System includes
#include <stdio.h>
#include <assert.h>

// CUDA runtime
#include <cuda_runtime.h>

// helper functions and utilities to work with CUDA
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cutil.h>
#include <fcntl.h>
#include <float.h>
#include <unistd.h>
#include "texton.h"
#include "convert.h"
#include "intervening.h"
#include "lanczos.h"
#include "stencilMVM.h"

#include "localcues.h"
#include "combine.h"
#include "nonmax.h"
#include "spectralPb.h"
#include "globalPb.h"
#include "skeleton.h"

#include "common_func.h"
#include "CircleTemplateTrace.h"
#include "LineSegTrace.h"

#include "pgm.cuh"

int main(int argc, char **argv)
{
	for (int i =0; i<= 74; i++)
	{
		int fileIdx = i;
		char filename[MAX_PATH];
		sprintf(filename, "%s%d", argv[1], fileIdx);
		char inputfile[MAX_PATH];
		char outputColorfile[MAX_PATH];

		char* period = strrchr(filename, '.');
		if (period == 0) {
			period = strrchr(filename, 0);
		}
		strncpy(inputfile, filename, period - filename);
		sprintf(&inputfile[0] + (period - filename) , "bin.pgm");

		strncpy(outputColorfile, filename, period - filename);
		sprintf(&outputColorfile[0] + (period - filename) , "_traced.pgm");

		/**null 从cutil内部申请内存*/
		//float* data = NULL;	
		int width, height;

		//cutLoadPGMf(inputfile, (float**)&data, &width, &height);
		float* data = loadPGM(inputfile, &width, &height);
		assert(width > 0 && height > 0);

		CLineSegTrace trace;
		int nAmount = trace.initTracePoints((float*)data, width, height);
		printf("Outputting %s \n", inputfile);
		printf("Valid Pt number:%d \n", nAmount);
		trace.traceLineSegs();
		trace.debugPrintOutput(outputColorfile, 70, 100000);

/*		cutFree(data);*/
	}

	//system("pause");
	return 1;
}