#ifndef LCUDA_PGM_H_
#define LCUDA_PGM_H_

float* loadPGM(const char *filename, int *width, int *height);

///[BUGFIX] 这个函数貌似会往文件里多写点东西
void savePGM(const char *filename, float* imageData, int width, int height);

#endif
