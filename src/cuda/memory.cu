#include <stdio.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

void alloc_memory_GPU(float **data_d, size_t size)
{
	gpuErrchk( cudaMalloc(data_d, size) );

}

void free_memory_GPU(float *cells_d, float *flags_d, float *vel_d, float *rho_d, float *fuerza_d, float *vertex_d, float *velocidad_d, float *velocidad2_d) {

	gpuErrchk( cudaFree(cells_d) );
	gpuErrchk( cudaFree(flags_d) );
	gpuErrchk( cudaFree(vel_d) );
	gpuErrchk( cudaFree(rho_d) );
	gpuErrchk( cudaFree(fuerza_d) );

	gpuErrchk( cudaFree(vertex_d) );
	gpuErrchk( cudaFree(velocidad_d) );
	gpuErrchk( cudaFree(velocidad2_d) );
}

void send_data_to_GPU(int X, int Y, int Z, float *cells, float *cells_d, float *flags, float *flags_d, float *vel, float *vel_d, float *rho, float *rho_d, float *fuerza, float *fuerza_d, int nNodos, float *vertex, float *vertex_d, float *velocidad, float *velocidad_d, float *velocidad2, float *velocidad2_d) {

	gpuErrchk( cudaMemcpy(cells_d, cells, 2*X*Y*Z*19*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(flags_d, flags, X*Y*Z*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(vel_d, vel, X*Y*Z*3*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(rho_d, rho, X*Y*Z*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(fuerza_d, fuerza, X*Y*Z*3*sizeof(float), cudaMemcpyHostToDevice) );

	gpuErrchk( cudaMemcpy(vertex_d, vertex, nNodos*3*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(velocidad_d, velocidad, nNodos*3*sizeof(float), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(velocidad2_d, velocidad2, nNodos*3*sizeof(float), cudaMemcpyHostToDevice) );

}

void retrieve_data_from_GPU(int X, int Y, int Z, float *cells, float *cells_d, float *flags, float *flags_d, float *vel, float *vel_d, float *rho, float *rho_d, float *fuerza, float *fuerza_d, int nNodos, float *vertex, float *vertex_d, float *velocidad, float *velocidad_d, float *velocidad2, float *velocidad2_d) {

	gpuErrchk( cudaMemcpy(cells, cells_d, 2*X*Y*Z*19*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(flags, flags_d, X*Y*Z*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(vel, vel_d, X*Y*Z*3*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(rho, rho_d, X*Y*Z*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(fuerza, fuerza_d, X*Y*Z*3*sizeof(float), cudaMemcpyDeviceToHost) );

	gpuErrchk( cudaMemcpy(vertex, vertex_d, nNodos*3*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(velocidad, velocidad_d, nNodos*3*sizeof(float), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(velocidad2, velocidad2_d, nNodos*3*sizeof(float), cudaMemcpyDeviceToHost) );

}
