#include "stdafx.h"
#include "csr.h"
#include "cuda.h"
#include "cutil.h"


CSRMatrix::CSRMatrix()
{
}

CSRMatrix::CSRMatrix(int mysize, int nonzeros)
{
    size = mysize;
    nnz = nonzeros;

    CUDA_SAFE_CALL(cudaMalloc(&rowptr, (size+1)*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc(&cols, (nnz)*sizeof(float)));
    CUDA_SAFE_CALL(cudaMalloc(&vals, (nnz)*sizeof(float)));
}

CSRMatrix::CSRMatrix(int mysize, int nonzeros, int* rowpointers, int* colindices, float *values)
{
    size = mysize;
    nnz = nonzeros;

    CUDA_SAFE_CALL(cudaMalloc(&rowptr, (size+1)*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc(&cols, (nnz)*sizeof(float)));
    CUDA_SAFE_CALL(cudaMalloc(&vals, (nnz)*sizeof(float)));
    
    CUDA_SAFE_CALL(cudaMemcpy(rowptr, rowpointers, (size+1)*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(cols, colindices, (nnz)*sizeof(float), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(vals, values, (nnz)*sizeof(float), cudaMemcpyHostToDevice));

}

__global__ void spmv(int* rowptr, int* cols, float* vals, float* vecin, float* vecout, int size)
{
    int x = blockIdx.x*blockDim.x+threadIdx.x;
    //int y = blockIdx.y*blockDim.y+threadIdx.y;
    
    int tid = x;

    if(x < size)
    {
        int start = rowptr[tid];
        int end = rowptr[tid+1]; 
        float sum = 0.0;

        for(int j=start; j<end; j++)
        {
            sum += vals[j]*vecin[cols[j]];
        }

        vecout[tid] = sum;
    }
}

void CSRMatrix::SpMV(float* gvecIn, float* gvecOut)
{        
    dim3 block = dim3(256);
    dim3 grid = dim3( (size+block.x-1)/block.x);
    spmv<<<grid, block>>>(rowptr, cols, vals, gvecIn, gvecOut, size);
}

__global__ void setValue_kernel(float* a, int length, float val)
{
    int gid = blockDim.x*blockIdx.x+threadIdx.x;

    if(gid<length) 
    {
        a[gid] = val;
    }   
}

void CSRMatrix::rowSum(float** p_rs)
{

    dim3 block = dim3(256);
    dim3 grid  = dim3( (size+block.x-1)/block.x);

    float* rs = *p_rs;
    if(!rs) 
    {
        //rs = new float[nPixels];
        CUDA_SAFE_CALL(cudaMalloc((void**)&rs, size*sizeof(float)));
        CUDA_SAFE_CALL(cudaMemset(rs, 0, size*sizeof(float)));
    }
    float* ones;
    CUDA_SAFE_CALL(cudaMalloc((void**)&ones, size*sizeof(float)));
    setValue_kernel<<<grid, block>>>(ones, size, 1.0);
        
    SpMV(rs, ones);
    CUDA_SAFE_CALL(cudaFree(ones));

    *p_rs = rs;
}

CSRMatrix::~CSRMatrix()
{
    CUDA_SAFE_CALL(cudaFree(rowptr));
    CUDA_SAFE_CALL(cudaFree(cols));
    CUDA_SAFE_CALL(cudaFree(vals));
}
