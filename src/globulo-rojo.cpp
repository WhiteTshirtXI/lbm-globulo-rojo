#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "ibm.h"
#include "fluid.h"
#include "mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
	mesh membrana, referencia;
	fluid fluido;
	double dt = 1.0;
	double dx = 1.0;
	double X = 21;
	double Y = 21;
	double Z = 21;
	int VTK = 50;

	// Parametros adimensionales
	double rho = 1.0;
	double nu = 1./6.;
	double Re = 0.5; 
	double G = 0.5;
	double R = Z/5;
	double gamma_dot = (Re*nu)/(rho*pow(R,2));
	double ks = (gamma_dot*nu*R)/(G);
	double kb = ks*1.0e-6;
	double kp = (gamma_dot)/(G);
	double STEPS = 12.0/kp;
	printf("A completar %f iteraciones\n", STEPS);

	// Membrana
	membrana.setID(1);
	membrana.mesh_refine_tri4();
	membrana.mesh_refine_tri4();
	membrana.mesh_refine_tri4();
	membrana.proyectarEsfera(R);
	membrana.proyectarRBC(R);
	membrana.moverCentro((X-1)/2.0, (Y-1)/2.0, (Z-1)/2.0);
	membrana.iniciarGeometria();

	referencia.setID(1);
	referencia.mesh_refine_tri4();
	referencia.mesh_refine_tri4();
	referencia.mesh_refine_tri4();
	referencia.proyectarEsfera(R);
	referencia.proyectarRBC(R);
	referencia.moverCentro((X-1)/2.0, (Y-1)/2.0, (Z-1)/2.0);
	referencia.iniciarGeometria();

	// Fluido
	fluido.inicializar(X,Y,Z);
	fluido.setVelocidad(gamma_dot);

	for(int ts = 0 ; ts < STEPS ; ts++)
		{
		// -----------------------------------------------------------------------//
		// 1. Interpolation
		// -----------------------------------------------------------------------//
		interpolation(fluido, membrana, X, Y, Z);

		// -----------------------------------------------------------------------//
		// 2. Encontrar nuevas posiciones de la membrana
		// -----------------------------------------------------------------------//
		membrana.moverNodos(dt, dx);

		// -----------------------------------------------------------------------//
		// 3. Calcular fuerzas en los nodos de la membrana
		// -----------------------------------------------------------------------//
		membrana.calcularFuerzasHelfrich(kb);
		membrana.calcularFuerzasFEM(referencia, ks);

		// -----------------------------------------------------------------------//
		// 4. Propagar la densidad de fuerza hacia el fluido
		// -----------------------------------------------------------------------//
		spread(fluido, membrana, X, Y, Z);

		// -----------------------------------------------------------------------//
		// 5. Solucionar la dinámica del fluido
		// -----------------------------------------------------------------------//
		fluido.collide();
		fluido.stream();

		// -----------------------------------------------------------------------//
		// 6. Calcular propiedades macro del fluido
		// -----------------------------------------------------------------------//
		fluido.calcularMacro();

		// -----------------------------------------------------------------------//
		// 7. Calcular propiedades macro de la membrana
		// -----------------------------------------------------------------------//
		membrana.calcularCambioArea(referencia);

		// -----------------------------------------------------------------------//
		// 9. Visualización
		// -----------------------------------------------------------------------//
		if(ts%VTK==0)
		{
			fluido.guardar(ts);
			membrana.guardarVTU(ts);
			printf("%d\n",ts);
		}
	}//Ciclo principal

	return 0;
}
