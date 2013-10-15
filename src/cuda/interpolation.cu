#include "helper.h"

/**
*	Función phi_4 para calcular función delta de Dirac con soporte cuatro
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación para distancia
*/
__device__ float phi_4_d(float r)
{
	float phi = 0.0;
	if((0.0 <= fabs(r)) && (fabs(r) <= 1.0))
        phi = ((1./8.)*(3.0-(2.0*fabs(r))+sqrt(1.0+(4.0*fabs(r))-(4.0*r*r))));
    if((1.0 <= fabs(r)) && (fabs(r) <= 2.0))
        phi = ((1./8.)*(5.0-(2.0*fabs(r))-sqrt(-7.0+(12.0*fabs(r))-(4.0*r*r))));
    if(fabs(r) >= 2.0)
        phi = 0.0;
    return phi;
}

/**
*	Función dirac_4 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
__device__ float dirac_4_d(float *x)
{
    float d = phi_4_d(x[0])*phi_4_d(x[1])*phi_4_d(x[2]);
    return d;
}

__global__ void interpolation(float *vel_d, int nNodos, float *vertex_d, float *velocidad_d, float *velocidad2_d, int X, int Y, int Z) {

	int u = blockIdx.x*blockDim.x + threadIdx.x;

	if (u < nNodos) {

		//Recorrer todos los nodos de la malla
		float pos[3], distancia[3], delta, a, A, b, B, c, C, ux=0.0, uy=0.0, uz=0.0;

		pos[0] = VERTEX_D(u, 0);
		pos[1] = VERTEX_D(u, 1);
		pos[2] = VERTEX_D(u, 2);

		a = pos[0]-3.0;
		A = pos[0]+3.0;
		b = pos[1]-3.0;
		B = pos[1]+3.0;
		c = pos[2]-3.0;
		C = pos[2]+3.0;


		for(int i = (int) a; i < A; i++)
			for(int j = (int) b; j < B; j++)
				for(int k=(int) c; k < C; k++)
				{
					distancia[0]=pos[0]-i;
					distancia[1]=pos[1]-j;
					distancia[2]=pos[2]-k;
					delta = dirac_4_d(distancia);
					ux+=delta*VEL_D(i, j, k, 0);
					uy+=delta*VEL_D(i, j, k, 1);
					uz+=delta*VEL_D(i, j, k, 2);
				}
		VELOCIDAD2_D(u, 0) = VELOCIDAD_D(u, 0);
		VELOCIDAD2_D(u, 1) = VELOCIDAD_D(u, 1);
		VELOCIDAD2_D(u, 2) = VELOCIDAD_D(u, 2);

		VELOCIDAD_D(u, 0) = ux;
		VELOCIDAD_D(u, 1) = uy;
		VELOCIDAD_D(u, 2) = uz;


	}
}


void interpolation_wrapper(float *vel_d, int nNodos, float *vertex_d, float *velocidad_d, float *velocidad2_d, int X, int Y, int Z) {


	//X*Y*Z = 9261;
	//Maximum number of threads per block:           1024

	dim3 grid_size;
	grid_size.x = nNodos/100 + 1;


	dim3 block_size;
	// 1000 threads per blocks
	block_size.x = 100;

	//Launch kernel
	interpolation<<<grid_size, block_size>>>(vel_d, nNodos, vertex_d, velocidad_d, velocidad2_d, X, Y, Z);


}
