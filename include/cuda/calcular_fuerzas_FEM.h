#ifndef _CALCULAR_FUERZAS_FEM_H
#define _CALCULAR_FUERZAS_FEM_H


void calcular_fuerzas_FEM_wrapper(int nNodos, int nCeldas, int *faces_d, float *vertex_d, float *vertex_ref_d, float ks, float *fuerza_mesh_d);


#endif
