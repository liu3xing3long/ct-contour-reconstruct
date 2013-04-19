#ifndef GRADIENT_H
#define GRADIENT_H
typedef unsigned int uint;

int initializeGradients(uint widthIn, uint heightIn, uint borderIn, uint maxbins, uint norientsIn, uint nscaleIn);

float* gradients(float* devImage, uint nbins, bool blur, float sigma, uint* radii);
float* gradients(int* devImage, uint nbins, bool blur, float sigma, uint* radii);


void finalizeGradients();

#endif
