#include <stdio.h>
#include <vector>
#include "../Stencil.h"

class Point {
 public:
 Point(int rowIn = 0, int colIn = 0, double valueIn = 0) {
   row = rowIn;
   col = colIn;
   value = valueIn;
 }
 int row;
 int col;
 double value;
};


int main(int argc, char** argv) {
  FILE* fp = fopen("tiny.sma", "w");
  int width = 4;
  int height = 4;
  int radius = 1;
  int nPixels = width * height;
  fwrite(&nPixels, sizeof(int), 1, fp);
  Stencil myStencil(radius, width, height, width);
  
  int nDiags = myStencil.getStencilArea();
  int* diags = (int*)malloc(sizeof(int)*nDiags);
  myStencil.copyOffsets(diags);
  for(int i = 1; i < nDiags; i++) {
    diags[i] = diags[i-1] + diags[i];
  }
  
  std::vector<Point> nonzeros;
  for(int diag = 0; diag < nDiags; diag++) {
    int differential = diags[diag];
    int xOffset = differential % width;
    int yOffset = differential / width;
    if (xOffset > radius) {
      xOffset = xOffset - width;
      yOffset++;
    }
    if (xOffset < -radius) {
      xOffset = xOffset + width;
      yOffset--;
    }
    
    double diagBase = (double)(diag + 1)/10;
    double diagInc = 0.0;
    for(int row = 0; row < height; row++) {
      for(int col = 0; col < width; col++) {
        int currentX = col + xOffset;
        int currentY = row + yOffset;
        if ((currentX >= 0) && (currentX < width) &&
            (currentY >= 0) && (currentY < height)) {
          int currentRow = row * width + col;
          int currentCol = currentY * width + currentX;
          if (currentCol >= currentRow) {
            nonzeros.push_back(Point(currentRow, currentCol, diagBase + diagInc));
            if (currentCol != currentRow) {
              nonzeros.push_back(Point(currentCol, currentRow, diagBase + diagInc));
            }
          }
          diagInc += 0.01;
          if (diagInc > 0.095) {
            diagInc = 0.0;
          }
        }
      }
    }    
  }

  
  
  int nnz = nonzeros.size();
  fwrite(&nnz, sizeof(int), 1, fp);
  int* n = new int[width * height];
  int* cols = new int[nnz];
  double* vals = new double[nnz];

  memset(n, 0, sizeof(int) * width * height);
  
  int* colPtr = cols;
  double* valPtr = vals;
  for(int index = 0; index < width * height; index++) {
    //printf("index: %i\n", index);
    for(int i = 0; i < nnz; i++) {
      Point currentPoint = nonzeros[i];
      if (index == currentPoint.row) {
        //printf("Row: %i, Col: %i, Val: %.2f\n", index, currentPoint.col, currentPoint.value); 
        n[index]++;
        *colPtr = currentPoint.col;
        *valPtr = currentPoint.value;
        colPtr++;
        valPtr++;
      }
    }
  }
  printf("Nonzeros found\n");

  
  int current = 0;
  for (int row = 0; row < nPixels; row++) {
    int nz = n[row];
    fwrite(&nz, sizeof(int), 1, fp);
    fwrite(&vals[current], sizeof(double), nz, fp);
    fwrite(&cols[current], sizeof(int), nz, fp);
    current = current + nz;
  }
  fclose(fp);
  printf("%i nonzeros found\n", nnz);
}
