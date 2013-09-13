#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "ibm.h"
#include "fluid.h"
#include "mesh.h"
#include "memory.h"

using namespace std;

int main(int argc, char *argv[])
{
	double X = 21;
	double Y = 21;
	double Z = 21;

	mesh membrana, referencia;
	fluid fluido(X, Y, Z);

	double dt = 1.0;
	double dx = 1.0;

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

	float *cells_f = fluido.get_cells();
	float *flags_f = fluido.get_flags();
	float *rho_f = fluido.get_rho();
	float *fuerza_f = fluido.get_fuerza();
	float *vel_f = fluido.get_vel();

	float *cells_d = NULL;
	float *flags_d = NULL;
	float *rho_d = NULL;
	float *fuerza_d = NULL;
	float *vel_d = NULL;

	alloc_memory_GPU(X, Y, Z, &cells_d, &flags_d, &vel_d, &rho_d, &fuerza_d);

	send_data_to_GPU(X, Y, Z, cells_f, cells_d, flags_f, flags_d, vel_f, vel_d, rho_f, rho_d, fuerza_f, fuerza_d);


	// Fluido
	//From fluid constructor
	fluido.calcularMacro(cells_d, rho_d, vel_d, fuerza_d);
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
		fluido.collide(cells_d, fuerza_d);
		fluido.stream(cells_d, flags_d);

		// -----------------------------------------------------------------------//
		// 6. Calcular propiedades macro del fluido
		// -----------------------------------------------------------------------//
		fluido.calcularMacro(cells_d, rho_d, vel_d, fuerza_d);

		// -----------------------------------------------------------------------//
		// 7. Calcular propiedades macro de la membrana
		// -----------------------------------------------------------------------//
		membrana.calcularCambioArea(referencia);

		// -----------------------------------------------------------------------//
		// 9. Visualización
		// -----------------------------------------------------------------------//
		if(ts%VTK==0)
		{
			retrieve_data_from_GPU(X, Y, Z, cells_f, cells_d, flags_f, flags_d, vel_f, vel_d, rho_f, rho_d, fuerza_f, fuerza_d);

			fluido.guardar(ts);
			membrana.guardarVTU(ts);
			printf("%d\n",ts);
		}
	}//Ciclo principal

	free_memory_GPU(cells_d, flags_d, vel_d, rho_d, fuerza_d);


	return 0;
}
