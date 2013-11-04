#include "helper.h"

__global__ void calcular_fuerzas_helfrich(int nNodos, float *fuerza_mesh_d) {

	int u = blockIdx.x*blockDim.x + threadIdx.x;

	if (u < nNodos) {

		FUERZA_MESH_D(u, 0) = 0.0;
		FUERZA_MESH_D(u, 1) = 0.0;
		FUERZA_MESH_D(u, 2) = 0.0;
		/*float c0 = 1.0e-9;
				float kh = darKhPorNodo(u);
				float kg = darKgPorNodo(u);
				float lkh = darLaplaceKh(u);
				float mag = 0.0;
				mag = (kb*(((2*kh) + c0)*(2*(kh*kh)-(2*kg)-(kh*c0)))) + (2*kb*lkh);
				FUERZA_MESH(u, 0) = FUERZA_MESH(u, 0)-mag*normalesPorNodo[u][0];
				FUERZA_MESH(u, 1) = FUERZA_MESH(u, 1)-mag*normalesPorNodo[u][1];
				FUERZA_MESH(u, 2) = FUERZA_MESH(u, 2)-mag*normalesPorNodo[u][2];*/

	}
}


void calcular_fuerzas_helfrich_wrapper(int nNodos, float *fuerza_mesh_d) {


	//X*Y*Z = 9261;
	//Maximum number of threads per block:           1024

	dim3 grid_size;
	grid_size.x = nNodos/100 + 1;


	dim3 block_size;
	// 1000 threads per blocks
	block_size.x = 100;

	//Launch kernel
	calcular_fuerzas_helfrich<<<grid_size, block_size>>>(nNodos, fuerza_mesh_d);


}
