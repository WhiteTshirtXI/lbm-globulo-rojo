#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fluid.h"
#include "mesh.h"
#include "ibm.h"
#include "helper.h"
#include "interpolation.h"

using namespace std;

/**
*	Función phi_2 para calcular función delta de Dirac con soporte de dos nodos
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación para distancia
*/
float phi_2(float r)
{
	float phi = 0.0;
	if ((0.0 <= fabs(r)) &&  (fabs(r)<= 1.0))
        phi = (1.0-fabs(r));
    if ( 1.0 <= fabs(r))
        phi = 0.0;
    return phi;
}

/**
*	Función phi_3 para calcular función delta de Dirac con soporte cuatro
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return float d, ponderación para distancia 700,  10100
*/
float phi_3(float r)
{
	float phi = 0.0;
	if((0.0 <= fabs(r)) && ( fabs(r) <= (1./2.)))
        phi = ((1./3.)*(1.0+sqrt(1.0-3.0*r*r)));
    if(((1./2.) <= fabs(r))&&(fabs(r) <= (3./2.)))
        phi = ((1./6.)*(5.0-3.0+fabs(r)-sqrt(-2.0+6.0*fabs(r)-3.0*r*r)));
    if((3./2.) <= fabs(r))
        phi = 0.0;
    return phi;
}

/**
*	Función phi_4 para calcular función delta de Dirac con soporte cuatro
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación para distancia
*/
float phi_4(float r)
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
float dirac_2(float *x)
{
    float d = phi_2(x[0])*phi_2(x[1])*phi_2(x[2]);
    return d;
}

/**
*	Función dirac_3 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
float dirac_3(float *x)
{
    float d = phi_3(x[0])*phi_3(x[1])*phi_3(x[2]);
    return d;
}

/**
*	Función dirac_4 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
float dirac_4(float *x)
{
    float d = phi_4(x[0])*phi_4(x[1])*phi_4(x[2]);
    return d;
}

/**
* TODO: Implementar las funcion spread
*/

void spread(fluid &fluido, mesh &membrana, int x, int y, int z)
{
	// Recorrer todos los nodos del fluido
	float pos[3]={0.0,0.0,0.0}, distancia[3]={0.0,0.0,0.0}, delta, df[3]={0.0,0.0,0.0}, fNodo[3]={0.0,0.0,0.0};
	float a, A, b, B, c, C;
	int nodos = membrana.darNumeroNodos();
	for(int u=0;u<nodos;u++)
	{
		membrana.darPosNodo(u, pos);
		membrana.darFuerzaNodo(u, fNodo);
		a = pos[0]-3.0;
		A = pos[0]+3.0;
		b = pos[1]-3.0;
		B = pos[1]+3.0;
		c = pos[2]-3.0;
		C = pos[2]+3.0;

		for(int i = (int) a;i<A;i++)
			for(int j = (int) b;j<B;j++)
				for(int k=(int) c;k<C;k++)
				{
					distancia[0]=pos[0]-i;
					distancia[1]=pos[1]-j;
					distancia[2]=pos[2]-k;
					delta = dirac_4(distancia);
					df[0]=fNodo[0]*delta;
					df[1]=fNodo[1]*delta;
					df[2]=fNodo[2]*delta;
					fluido.addFuerza(i,j,k,df);
				}
		df[0] = 0.0;
		df[1] = 0.0;
		df[2] = 0.0;
	}
}

/**
* TODO: Implementar las funcion interpolation
*/
void interpolation(float *vel_d, float *vertex_d, float *velocidad_d, float *velocidad2_d, int nNodos, int X, int Y, int Z)
{

	//vel fluido
	//membrana needs nNodos, vertex, velocidad, velocidad2

	interpolation_wrapper(vel_d, nNodos, vertex_d, velocidad_d, velocidad2_d, X, Y, Z);

}
