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

void alloc_memory_GPU(void **data_d, size_t size)
{
	gpuErrchk( cudaMalloc(data_d, size) );

}

void free_memory_GPU(void *data_d) {

	gpuErrchk( cudaFree(data_d) );
}

void send_data_to_GPU(void *data, void *data_d, size_t size) {

	gpuErrchk( cudaMemcpy(data_d, data, size, cudaMemcpyHostToDevice) );

}

void retrieve_data_from_GPU(void *data, void *data_d, size_t size) {

	gpuErrchk( cudaMemcpy(data, data_d, size, cudaMemcpyDeviceToHost) );

}
