#ifndef windows_acml_h__
#define windows_acml_h__

#define LAPACK_ORIG

#ifdef LAPACK_ORIG
#include <complex>
/// Custom Lapacke complex type
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float    float 
#define lapack_complex_double   double
#include "..\lapacke\include\lapacke_mangling.h"
#include "..\lapacke\include\lapacke.h"

#else


extern "C"{
#include "..\clapack\include\f2c.h"
#include "..\clapack\include\clapack.h"
	
	int sgetrf_(integer *m, integer *n, real *a, integer *lda, 
		integer *ipiv, integer *info);
	int sgetri_(integer *n, real *a, integer *lda, integer *ipiv,
		real *work, integer *lwork, integer *info);
	int dstein_(integer *n, doublereal *d__, doublereal *e, 
		integer *m, doublereal *w, integer *iblock, integer *isplit, 
		doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
		integer *ifail, integer *info);
	int dstebz_(char *range, char *order, integer *n, doublereal 
		*vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
		doublereal *d__, doublereal *e, integer *m, integer *nsplit, 
		doublereal *w, integer *iblock, integer *isplit, doublereal *work, 
		integer *iwork, integer *info);
};

#define  LAPACK_sgetrf sgetrf_
#define  LAPACK_sgetri sgetri_
#define  LAPACK_dstein dstein_
#define  LAPACK_dstebz dstebz_

#endif



#endif // windows_acml_h__
