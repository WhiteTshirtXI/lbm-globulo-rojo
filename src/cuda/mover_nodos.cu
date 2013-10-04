#include "helper.h"

__global__ void mover_nodos(int dt, int nNodos, float *vertex_d, float *velocidad_d, float *velocidad2_d) {

	int u = blockIdx.x*blockDim.x + threadIdx.x;

	if (u < nNodos) {

		VERTEX_D(u, 0) += ((3/2)*(VELOCIDAD_D(u, 0)) - (1./2.0)*(VELOCIDAD2_D(u, 0)))*dt;
		VERTEX_D(u, 1) += ((3/2)*(VELOCIDAD_D(u, 1)) - (1./2.0)*(VELOCIDAD2_D(u, 1)))*dt;
		VERTEX_D(u, 2) += ((3/2)*(VELOCIDAD_D(u, 2)) - (1./2.0)*(VELOCIDAD2_D(u, 2)))*dt;

	}
}


void mover_nodos_wrapper(int dt, int nNodos, float *vertex_d, float *velocidad_d, float *velocidad2_d) {


	//X*Y*Z = 9261;
	//Maximum number of threads per block:           1024

	dim3 grid_size;
	grid_size.x = nNodos/1024 + 1;


	dim3 block_size;
	// 1000 threads per blocks
	block_size.x = 1024;

	//Launch kernel
	mover_nodos<<<grid_size, block_size>>>(dt, nNodos, vertex_d, velocidad_d, velocidad2_d);


}
