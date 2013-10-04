#ifndef _MEMORY_H_
#define _MEMORY_H_

void alloc_memory_GPU(float **data_d, size_t size);
void free_memory_GPU(float *cells_d, float *flags_d, float *vel_d, float *rho_d, float *fuerza_d, float *vertex_d, float *velocidad_d, float *velocidad2_d);

void send_data_to_GPU(int X, int Y, int Z, float *cells, float *cells_d, float *flags, float *flags_d, float *vel, float *vel_d, float *rho, float *rho_d, float *fuerza, float *fuerza_d, int nNodos, float *vertex, float *vertex_d, float *velocidad, float *velocidad_d, float *velocidad2, float *velocidad2_d);
void retrieve_data_from_GPU(int X, int Y, int Z, float *cells, float *cells_d, float *flags, float *flags_d, float *vel, float *vel_d, float *rho, float *rho_d, float *fuerza, float *fuerza_d, int nNodos, float *vertex, float *vertex_d, float *velocidad, float *velocidad_d, float *velocidad2, float *velocidad2_d);


#endif
