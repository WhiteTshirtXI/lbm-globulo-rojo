#include "helper.h"

/**
 *	Función phi_4 para calcular función delta de Dirac con soporte cuatro
 *	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
 *   @return d, ponderación para distancia
 */
__device__ float phi_4_e(float r)
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
__device__ float dirac_4_e(float *x)
{
	float d = phi_4_e(x[0])*phi_4_e(x[1])*phi_4_e(x[2]);
	return d;
}

__global__ void spread(int nNodos, float *vertex_d, float *fuerza_mesh_d, float *fuerza_d, int X, int Y, int Z) {

	int u = blockIdx.x*blockDim.x + threadIdx.x;

	if (u < nNodos) {

		float pos[3]={0.0,0.0,0.0}, distancia[3]={0.0,0.0,0.0}, delta, df[3]={0.0,0.0,0.0}, fNodo[3]={0.0,0.0,0.0};
		float a, A, b, B, c, C;

		pos[0] = VERTEX_D(u, 0);
		pos[1] = VERTEX_D(u, 1);
		pos[2] = VERTEX_D(u, 2);

		fNodo[0] = FUERZA_MESH_D(u, 0);
		fNodo[1] = FUERZA_MESH_D(u, 1);
		fNodo[2] = FUERZA_MESH_D(u, 2);

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
					delta = dirac_4_e(distancia);
					df[0]=fNodo[0]*delta;
					df[1]=fNodo[1]*delta;
					df[2]=fNodo[2]*delta;
					FUERZA_D(i, j, k, 0) += df[0];
					FUERZA_D(i, j, k, 1) += df[1];
					FUERZA_D(i, j, k, 2) += df[2];
				}
		df[0] = 0.0;
		df[1] = 0.0;
		df[2] = 0.0;
	}

}



void spread_wrapper(int nNodos, float *vertex_d, float *fuerza_mesh_d, float *fuerza_d, int X, int Y, int Z) {

	//X*Y*Z = 9261;
	//Maximum number of threads per block:           1024

	dim3 grid_size;
	grid_size.x = nNodos/1024 + 1;


	dim3 block_size;
	// 1000 threads per blocks
	block_size.x = 1024;

	//Launch kernel
	spread<<<grid_size, block_size>>>(nNodos, vertex_d, fuerza_mesh_d, fuerza_d, X, Y, Z);


}
