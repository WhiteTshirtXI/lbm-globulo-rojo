#include "helper.h"
#include "math.h"

__device__ float norm_d(float a[3])
{
	return sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]));
}

__device__ void cross_d(float a[3], float b[3], float temp[3])
{
		temp[0] = a[1]*b[2] - a[2]*b[1];
		temp[1] = a[2]*b[0] - a[0]*b[2];
		temp[2] = a[0]*b[1] - a[1]*b[0];
}

__device__ float darAreaElemento(int i, int *faces_d, float *vertex_d, int nNodos, int nCeldas) {

	float a[3], b[3], c[3], v1[3], v2[3], temp[3];
	int A, B, C;
	A=FACES_D(i, 0);
	B=FACES_D(i, 1);
	C=FACES_D(i, 2);

	a[0]=VERTEX_D(A, 0);
	a[1]=VERTEX_D(A, 1);
	a[2]=VERTEX_D(A, 2);

	b[0]=VERTEX_D(B, 0);
	b[1]=VERTEX_D(B, 1);
	b[2]=VERTEX_D(B, 2);

	c[0]=VERTEX_D(C, 0);
	c[1]=VERTEX_D(C, 1);
	c[2]=VERTEX_D(C, 2);

	v1[0] = b[0] - a[0];
	v1[1] = b[1] - a[1];
	v1[2] = b[2] - a[2];

	v2[0] = c[0] - a[0];
	v2[1] = c[1] - a[1];
	v2[2] = c[2] - a[2];

	cross_d(v1, v2, temp);
	return norm_d(temp)/2.0;
}

__global__ void calcular_cambio_area(int nCeldas, int nNodos, int *faces_d, float *vertex_d, int *faces_ref_d, float *vertex_ref_d, float *area_d) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < nCeldas) {

		float a = darAreaElemento(i, faces_d, vertex_d, nNodos, nCeldas);
		float b = darAreaElemento(i, faces_ref_d, vertex_ref_d, nNodos, nCeldas);

		area_d[i] = a - b;


	}
}


void calcular_cambio_area_wrapper(int nCeldas, int nNodos, int *faces_d, float *vertex_d, int *faces_ref_d, float *vertex_ref_d, float *area_d) {


	//X*Y*Z = 9261;
	//Maximum number of threads per block:           1024

	dim3 grid_size;
	grid_size.x = nCeldas/1024 + 1;


	dim3 block_size;
	// 1000 threads per blocks
	block_size.x = 1024;

	//Launch kernel
	calcular_cambio_area<<<grid_size, block_size>>>(nCeldas, nNodos, faces_d, vertex_d, faces_ref_d, vertex_ref_d, area_d);


}
