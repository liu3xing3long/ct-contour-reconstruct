#ifndef CSR
#define CSR

class CSRMatrix
{
    public:
        int * rowptr;
        int * cols;
        float * vals;

        int size;
        int nnz;

    public:
        
        CSRMatrix();
        CSRMatrix(int mysize, int nnzeros);
        CSRMatrix(int mysize, int nnzeros, int *rowpointers, int *colindices, float *values);

        void rowSum(float** p_rs);
        void SpMV(float* vecin, float* vecout);

        ~CSRMatrix();
};

#endif
