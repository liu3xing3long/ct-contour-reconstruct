#include "stdafx.h"
#include "junction_cpp.h"

void junctionComputation( float* fControus, int width, int height )
{
	float* fWeight = new float[height*width];
	memcpy(fWeight, fControus, width*height*sizeof(float));




	delete fWeight;
	fWeight = NULL;
}
