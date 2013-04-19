#include <stdio.h>
#include <vector>

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
  int width = 17;
  int height = 17;
  int radius = 1;
  int nPixels = width * height;
  fwrite(&nPixels, sizeof(int), 1, fp);
  
  int nDiags = 5;
  int diagX[] = {0, -1, 0, 1, 0};
  int diagY[] = {-1, 0, 0, 0, 1};
  std::vector<Point> nonzeros;
  for(int diag = 0; diag < nDiags; diag++) {
    int xOffset = diagX[diag];
    int yOffset = diagY[diag];
    double diagBase = (double)(diag + 1)/10;
    double diagInc = 0.0;
    for(int row = 0; row < height; row++) {
      for(int col = 0; col < width; col++) {
        int currentX = col + xOffset;
        int currentY = row + yOffset;
        if ((currentX >= 0) && (currentX < width) &&
            (currentY >= 0) && (currentY < width)) {
          int currentRow = row * width + col;
          int currentCol = currentY * width + currentX;
          nonzeros.push_back(Point(currentRow, currentCol, diagBase + diagInc));
          diagInc += 0.01;
          if (diagInc >= 0.09) {
            diagInc = 0.0;
          }
        }
      }
    }    
  }
  

  
  /* int nnz = 64; */
/*   int n[] = {3, 4, 4, 3, 4, 5, 5, 4, 4, 5, 5, 4, 3, 4, 4, 3}; */
/*   int cols[] = {      0, 1, 4, */
/*                    0, 1, 2, 5, */
/*                    1, 2, 3, 6, */
/*                    2, 3,    7, */
/*                 0,    4, 5, 8, */
/*                 1, 4, 5, 6, 9, */
/*                 2, 5, 6, 7,10, */
/*                 3, 6, 7,   11, */
/*                 4,    8, 9,12, */
/*                 5, 8, 9,10,13, */
/*                 6, 9,10,11,14, */
/*                 7,10,11,   15, */
/*                 8,   12,13, */
/*                 9,12,13,14, */
/*                10,13,14,15, */
/*                11,14,15}; */
/*   double vals[] = {            0.31, 0.41, 0.51, */
/*                          0.21, 0.32, 0.42, 0.52, */
/*                          0.22, 0.33, 0.43, 0.53, */
/*                          0.23, 0.34,       0.54, */
/*                    0.11,       0.35, 0.44, 0.55, */
/*                    0.12, 0.24, 0.36, 0.45, 0.56, */
/*                    0.13, 0.25, 0.37, 0.46, 0.57, */
/*                    0.14, 0.26, 0.38,       0.58, */
/*                    0.15,       0.39, 0.47, 0.59, */
/*                    0.16, 0.27, 0.31, 0.48, 0.51, */
/*                    0.17, 0.28, 0.32, 0.49, 0.52, */
/*                    0.18, 0.29, 0.33,       0.53, */
/*                    0.19,       0.34, 0.41,  */
/*                    0.11, 0.21, 0.35, 0.42,  */
/*                    0.12, 0.22, 0.36, 0.43,  */
/*                    0.13, 0.23, 0.37}; */
  int nnz = nonzeros.size();
  fwrite(&nnz, sizeof(int), 1, fp);
  int* n = new int[width * height];
  int* cols = new int[nnz];
  double* vals = new double[nnz];

  memset(n, 0, sizeof(int) * width * height);
  
  int* colPtr = cols;
  double* valPtr = vals;
  for(int index = 0; index < width * height; index++) {
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
