#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mesh.h"
#include "elemento.h"
#include <omp.h>
#include "debug.h"
#include "helper.h"

float* mesh::get_vertex(void) {

	return vertex;

}

void mesh::print() {

	printf("id: %d\n", id);
	printf("nNodos: %d\n", nNodos);
	printf("nCeldas: %d\n", nCeldas);
	printf("cX: %f, cY: %f, cZ: %f\n", cX, cY, cZ);

	printf("VERTEX: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", VERTEX(i, j));
	printf("\n");

	/*
	printf("VELOCIDAD: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", VELOCIDAD(i, j));
	printf("\n");

	printf("VELOCIDAD2: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", VELOCIDAD2(i, j));
	printf("\n");

	printf("FUERZA: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", FUERZA_MESH(i, j));
	printf("\n");

	printf("AREA: ");
	for(int i = 0; i < nCeldas; i++)
		printf("%f ", area[i]);
	printf("\n");
*/
	printf("FACES: ");
	for(int i = 0; i < nCeldas; i++)
		for(int j = 0; j < 3; j++)
			printf("%d ", FACES(i, j));
	printf("\n");
/*
	printf("NORMALESPORNODO: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", NORMALESPORNODO(i, j));
	printf("\n");

	printf("NORMALESPORCARA: ");
	for(int i = 0; i < nCeldas; i++)
		for(int j = 0; j < 3; j++)
			printf("%f ", NORMALESPORCARA(i, j));
	printf("\n");


	printf("CARASPORNODO: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 7; j++)
			printf("%f ", CARASPORNODO(i, j));
	printf("\n");

	printf("ANGULOSPORNODO: ");
	for(int i = 0; i < nNodos; i++)
		for(int j = 0; j < 7; j++)
			printf("%f ", ANGULOSPORNODO(i, j));
	printf("\n");

	printf("LAPLACEKG: ");
	for(int i = 0; i < nNodos; i++)
		printf("%f ", laplaceKg[i]);
	printf("\n");

	printf("LAPLACEKH: ");
	for(int i = 0; i < nNodos; i++)
		printf("%f ", laplaceKh[i]);
	printf("\n");

	printf("NODOSPROBLEMAS: ");
	for(int i = 0; i < 12; i++)
		printf("%f ", nodosProblema[i]);
	printf("\n");
*/
	printf("areaS: %f\n", areaS);
	printf("volumenE: %f\n", volumenE);

}

// Constructor
mesh::mesh()
{

	cX = cY = cZ = 0;

	float t = (1.+sqrt(5.))/2.;
	float tau = t/sqrt(1.+t*t);
	float one = 1./sqrt(1.+t*t); // Unit sphere

	// Twelve vertices of icosahedron on unit sphere
	// Creacion dinamica de arreglos al apuntador

	nNodos = 12;

	vertex = (float*)malloc(12*3*sizeof(float));
	if (vertex == NULL) {

		_DEBUG("Error allocating vertex");
		exit(-1);
	}

	float nodos[12][3] ={{  tau,  one,    0},
			{ -tau,  one,    0 },
			{ -tau, -one,    0 },
			{  tau, -one,    0 },
			{  one,   0 ,  tau },
			{  one,   0 , -tau },
			{ -one,   0 , -tau },
			{ -one,   0 ,  tau },
			{   0 ,  tau,  one },
			{   0 , -tau,  one },
			{   0 , -tau, -one },
			{   0 ,  tau, -one }};

	for(int i = 0; i < nNodos; i++)
		for(int j = 0 ; j < 3; j++)
			VERTEX(i, j) = nodos[i][j];



	// Structure for unit icosahedron
	nCeldas = 20;

	faces = (int*)malloc(20*3*sizeof(int));
	if (faces == NULL) {

		_DEBUG("Error allocating faces");
		exit(-1);
	}

	int f[20][3]= {{4,  7,  8},
			{4,  9,  7},
			{5, 11,  6},
			{5,  6, 10},
			{0,  3,  4},
			{0,  5,  3},
			{2,  1,  7},
			{2,  6,  1},
			{8, 11,  0},
			{8,  1, 11},
			{9,  3, 10},
			{9, 10,  2},
			{8,  0,  4},
			{11,  5, 0},
			{4,  3, 9},
			{5, 10,  3},
			{7,  1,  8},
			{6, 11,  1},
			{7, 9,  2},
			{6,  2, 10}};

	for(int i = 0; i < nCeldas; i++)
		for(int j = 0 ; j < 3; j++)
			FACES(i, j) = f[i][j];

}

//Destructor: Destruye los elementos de la malla nodos y celdas
mesh::~mesh()
{

	if (vertex != NULL) free(vertex);
	if (faces != NULL) free(faces);

	if (normalesPorCara != NULL) free(normalesPorCara);
	if (carasPorNodo != NULL) free(carasPorNodo);
	if (normalesPorNodo != NULL) free(normalesPorNodo);
	if (angulosPorNodo != NULL) free(angulosPorNodo);
	if (fuerza != NULL) free(fuerza);
	if (velocidad != NULL) free(velocidad);
	if (velocidad2 != NULL) free(velocidad2);
}

// Trasladar en el espacio el centro de la esfera
void mesh::moverCentro(float x, float y, float z)
{
	cX = x;
	cY = y;
	cZ = z;

	for(int i = 0; i < nNodos; i++)
	{
		VERTEX(i, 0) = VERTEX(i, 0) + x;
		VERTEX(i, 1) = VERTEX(i, 1) + y;
		VERTEX(i, 2) = VERTEX(i, 2) + z;
	}
}

// Rotar la estructura completa dados tres angulos
// theta, phi, psi
void mesh::rotarEstructura(float alpha, float phi, float theta)
{
	float Rx[3][3] = {{1,          0,          0},
			{0, cos(alpha), -sin(alpha)},
			{0, sin(alpha),  cos(alpha)}};

	float Ry[3][3] = {{ cos(phi), 0, sin(phi)},
			{        0, 1,        0},
			{-sin(phi), 0, cos(phi)}};

	float Rz[3][3] = {{cos(theta), -sin(theta),0},
			{sin(theta), cos(theta), 0},
			{0, 0, 1}};
	float Rxy[3][3];
	float Rxyz[3][3];

	// Realizar composición de rotaciones
	/* Realiza el producto de matrices y guarda el resultado en una tercera matriz*/
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			Rxy[i][j]=0;
			for(int k=0;k<3;k++){
				Rxy[i][j]=Rxy[i][j]+(Rx[i][k]*Ry[k][j]);
			}
		}
	}

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			Rxyz[i][j]=0;
			for(int k=0;k<3;k++){
				Rxyz[i][j]=Rxyz[i][j]+(Rxy[i][k]*Rz[k][j]);
			}
		}
	}

	// Rotar cada uno de los nodos por la matriz
	for(int i = 0; i<nNodos ; i++)
	{
		float x = Rxyz[0][0]*VERTEX(i, 0) + Rxyz[0][1]*VERTEX(i, 1) + Rxyz[0][2]*VERTEX(i, 2);
		float y = Rxyz[1][0]*VERTEX(i, 0) + Rxyz[1][1]*VERTEX(i, 1) + Rxyz[1][2]*VERTEX(i, 2);
		float z = Rxyz[2][0]*VERTEX(i, 0) + Rxyz[2][1]*VERTEX(i, 1) + Rxyz[2][2]*VERTEX(i, 2);
		VERTEX(i, 0) = x;
		VERTEX(i, 1) = y;
		VERTEX(i, 2) = z;
	}

}


// Proyectar los nodos a una superficie RCB
void mesh::proyectarRBC(float r)
{

	float c0 = 0.207;
	float c1 = 2.000;
	float c2 = -1.123;
	float R  = r;

	for(int i = 0;i<nNodos;i++)
	{
		float X = VERTEX(i, 0);
		float Y = VERTEX(i, 1);
		float a = ((X*X)+(Y*Y))/(R*R);
		float b = (c0 + c1*a + c2*(a*a));
		float Z = (0.5*R*(sqrt(fabs(1.0-a))))*b;
		VERTEX(i, 0)=X;
		VERTEX(i, 1)=Y;

		// Calcular magnitud del vector
		float x, y, z;
		x = VERTEX(i, 0)*VERTEX(i, 0);
		y = VERTEX(i, 1)*VERTEX(i, 1);
		z = VERTEX(i, 2)*VERTEX(i, 2);
		float mag = sqrt(x+y+z);
		if(VERTEX(i, 2)/mag<0){
			VERTEX(i, 2)=-Z;// + (X*X)/16;
		}else{
			VERTEX(i, 2)=Z;// + (X*X)/16;
		}
	}
}


// Proyectar los nodos a una superficie esferica
void mesh::proyectarEsfera(float r)
{
	for(int i=0; i<nNodos;i++)
	{
		float x, y, z;
		x = VERTEX(i, 0)*VERTEX(i, 0);
		y = VERTEX(i, 1)*VERTEX(i, 1);
		z = VERTEX(i, 2)*VERTEX(i, 2);
		float mag= sqrt(x+y+z);
		VERTEX(i, 0) = VERTEX(i, 0)*r/mag;
		VERTEX(i, 1) = VERTEX(i, 1)*r/mag;
		VERTEX(i, 2) = VERTEX(i, 2)*r/mag;
	}
}

// Proyectar los nodos a una superficie
void mesh::proyectarElipsoide(float a, float b, float c)
{
	for(int i=0; i<nNodos;i++)
	{
		float x, y, z;
		x = VERTEX(i, 0)*VERTEX(i, 0);
		y = VERTEX(i, 1)*VERTEX(i, 1);
		z = VERTEX(i, 2)*VERTEX(i, 2);
		float mag= sqrt(x+y+z);
		VERTEX(i, 0) = VERTEX(i, 0)*a/mag;
		VERTEX(i, 1) = VERTEX(i, 1)*b/mag;
		VERTEX(i, 2) = VERTEX(i, 2)*c/mag;
	}
}


// Funcion que agrega un nodo a la lista de nodos de la malla
// Retorna la posicion del nuevo nodo
int mesh::agregarNodo(float x, float y, float z)
{

	float *nuevaLista = (float*)malloc((nNodos+1)*3*sizeof(float));
	if (nuevaLista == NULL) {

		_DEBUG("Error allocating nuevaLista");
		exit(-1);
	}
	memset(nuevaLista, 0, (nNodos+1)*3*sizeof(float));

	for(int i = 0; i < nNodos; i++)
	{
		ACCESS2(nuevaLista, nNodos+1, 3, i, 0) = VERTEX(i, 0);
		ACCESS2(nuevaLista, nNodos+1, 3, i, 1) = VERTEX(i, 1);
		ACCESS2(nuevaLista, nNodos+1, 3, i, 2) = VERTEX(i, 2);
	}

	ACCESS2(nuevaLista, nNodos+1, 3, nNodos, 0) = x;
	ACCESS2(nuevaLista, nNodos+1, 3, nNodos, 1) = y;
	ACCESS2(nuevaLista, nNodos+1, 3, nNodos, 2) = z;


	free(vertex);
	vertex = nuevaLista;

	nNodos=nNodos+1;

	return nNodos-1;
}

// Funcion que determina si un nodo existe en la lista
// retorna -1 si no existe
// retorna  N si existe
int mesh::existeNodo(float x, float y, float z)
{
	int existe = -1;
	for(int i=0;i<nNodos;i++)
	{
		bool a = (VERTEX(i, 0)==x);
		bool b = (VERTEX(i, 1)==y);
		bool c = (VERTEX(i, 2)==z);
		if(a &&  b && c)
		{
			existe = i;
		}
	}
	return existe;
}


//Función para refinamiento de malla
void mesh::mesh_refine_tri4()
{
	/*
	// Dividir cada cara en 4 triangulos
	// Make new midpoints
	// a = (A+B)/2
	// b = (B+C)/2
	// c = (C+A)/2
	//        B
	//       /\
	//      /  \
	//    a/____\b       Construct new triangles
	//    /\    /\       [A,a,c]
	//   /  \  /  \      [a,B,b]
	//  /____\/____\     [c,b,C]
	// A      c     C    [a,b,c]
	// -------------------------------------------
	 *
	 */

	printf("Refinando...\n");
	// Crear la estructura para las nuevas celdas

	int *nuevaLista = (int*)malloc((4*nCeldas)*3*sizeof(int));
	if (nuevaLista == NULL) {

		_DEBUG("Error allocating nuevaLista");
		exit(-1);
	}
	memset(nuevaLista, 0, (4*nCeldas)*3*sizeof(int));

	// Dividir todas las caras
	for(int f = 0; f < nCeldas; f++)
	{
		int NA, NB, NC;
		// Encontrar los vertices correspondientes a cada cara
		NA = FACES(f, 0);
		NB = FACES(f, 1);
		NC = FACES(f, 2);

		// Encontrar las coordenadas de cada vertice
		float A[3], B[3], C[3];
		A[0]=VERTEX(NA, 0);
		A[1]=VERTEX(NA, 1);
		A[2]=VERTEX(NA, 2);

		B[0]=VERTEX(NB, 0);
		B[1]=VERTEX(NB, 1);
		B[2]=VERTEX(NB, 2);

		C[0]=VERTEX(NC, 0);
		C[1]=VERTEX(NC, 1);
		C[2]=VERTEX(NC, 2);

		// Hallar los puntos medios de cada coordenada
		float a[3], b[3], c[3];
		a[0] = (A[0]+B[0])/2.;
		a[1] = (A[1]+B[1])/2.;
		a[2] = (A[2]+B[2])/2.;

		b[0] = (B[0]+C[0])/2.;
		b[1] = (B[1]+C[1])/2.;
		b[2] = (B[2]+C[2])/2.;

		c[0] = (C[0]+A[0])/2.;
		c[1] = (C[1]+A[1])/2.;
		c[2] = (C[2]+A[2])/2.;

		// Hallar el indice de cada nodo de las nuevas celdas
		int Na, Nb, Nc;
		Na = posicionNodo(a[0], a[1], a[2]);
		Nb = posicionNodo(b[0], b[1], b[2]);
		Nc = posicionNodo(c[0], c[1], c[2]);

		// Crear las nuevas caras
		int C1[3]={NA, Na, Nc}, C2[3]={Na, NB, Nb};
		int C3[3]={Nc, Nb, NC}, C4[3]={Na, Nb, Nc};

		// Agregarlas a la nueva estructura de celdas
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-4, 0) = C1[0];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-4, 1) = C1[1];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-4, 2) = C1[2];

		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-3, 0) = C2[0];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-3, 1) = C2[1];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-3, 2) = C2[2];

		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-2, 0) = C3[0];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-2, 1) = C3[1];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-2, 2) = C3[2];

		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-1, 0) = C4[0];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-1, 1) = C4[1];
		ACCESS2(nuevaLista, 4*nCeldas, 3, (f+1)*4-1, 2) = C4[2];
	}

	free(faces);
	faces = nuevaLista;

	nCeldas = 4*nCeldas;
}

// Funcion auxiliar para obtener la posicion de los vertices
int mesh::posicionNodo(float x, float y, float z)
{
	//Decidir si existe o no existe en la lista
	int pos = existeNodo(x,y,z);
	if(pos>-1)
	{
		return pos;
	}
	return agregarNodo(x,y,z);
}

/**

 */
int mesh::guardarVTU(int t)
{
	FILE *archivo;/*El manejador de archivo*/
	char ruta[80];
	char numero[10];
	char nIden[10];
	char iden[10];

	sprintf(nIden,"%i",id);
	strcpy(iden, "cel");
	strcat(iden,nIden);
	strcat(iden,"-");

	strcpy(ruta, "temp/esfera-");
	sprintf(numero,"%d",t);
	strcat(ruta,iden);
	strcat(ruta,numero);

	//Agrega el numero de la celula
	strcat(ruta,".vtk");
	archivo=fopen(ruta, "w");
	if(archivo==NULL){/*Si no lo logramos abrir, salimos*/
		printf("No se puede guardar archivo");
		return 1;}
	else{

		// Escribir datos al archivo
		// 1. Escribir cabecera.
		fprintf(archivo, "# vtk DataFile Version 3.0\n");
		fprintf(archivo, "vtk output\n");
		fprintf(archivo, "ASCII\n");
		fprintf(archivo, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(archivo, "POINTS %i float\n",nNodos);

		// Escribir coordenadas de cada nodo
		for(int i = 0; i<nNodos ; i++){
			for(int j = 0; j<3 ; j++){
				fprintf(archivo, "%f ", VERTEX(i, j));
			}
			fprintf(archivo, "\n");}

		// Escribir celdas
		fprintf(archivo, "CELLS %d %d\n",nCeldas,nCeldas*4);
		for(int i = 0; i<nCeldas ; i++){
			fprintf(archivo, "3 ");
			for(int j = 0; j<3 ; j++){
				fprintf(archivo, "%i ", FACES(i, j));
			}
			fprintf(archivo, "\n");}
		fprintf(archivo,"CELL_TYPES %d\n",nCeldas);
		for(int i = 0; i<nCeldas ; i++){
			fprintf(archivo, "5 \n");}

		// Escribir velocidad
		fprintf(archivo,"POINT_DATA %d\n", nNodos);
		fprintf(archivo,"VECTORS Velocidad float\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f %f %f\n", VELOCIDAD(i, 0), VELOCIDAD(i,1), VELOCIDAD(i, 2));
		}

		// Escribir vectores normales por cada nodo
		fprintf(archivo,"VECTORS Normales float\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f %f %f\n", NORMALESPORNODO(i, 0), NORMALESPORNODO(i, 1), NORMALESPORNODO(i, 2));
		}

		// Escribir vectores de fuerza por cada nodo
		fprintf(archivo,"VECTORS Fuerzas float\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f %f %f\n", FUERZA_MESH(i, 0), FUERZA_MESH(i, 1), FUERZA_MESH(i, 2));
		}

		// Escribir valores de curvatura media en cada nodo
		fprintf(archivo,"SCALARS Curvatura-media float\n");
		fprintf(archivo,"LOOKUP_TABLE default\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f\n", darKhPorNodo(i));
		}

		// Escribir valores de curvatura de Gauss en cada nodo
		fprintf(archivo,"SCALARS Curvatura-Gauss float\n");
		fprintf(archivo,"LOOKUP_TABLE default\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f\n", darKgPorNodo(i));
		}

		// Escribir valores de curvatura de Gauss en cada nodo
		fprintf(archivo,"SCALARS k1 float\n");
		fprintf(archivo,"LOOKUP_TABLE default\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f\n", darK1PorNodo(i));
		}

		// Escribir valores de curvatura de Gauss en cada nodo
		fprintf(archivo,"SCALARS k2 float\n");
		fprintf(archivo,"LOOKUP_TABLE default\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f\n", darK1PorNodo(i));
		}

		// Escribir los valores del operador Laplace en la curvatura media
		fprintf(archivo,"SCALARS OLB-Kh float\n");
		fprintf(archivo,"LOOKUP_TABLE default\n");
		for(int i=0;i<nNodos;i++)
		{
			fprintf(archivo, "%f\n", darLaplaceKh(i));
		}
		fclose(archivo);
		return 0;
	}
}

void mesh::iniciarGeometria()
{
	// 1. Iniciar el área superficial
	areaS = calcularAreaSuperficial();

	// 2. Iniciar el volumen encerrado por la membrana
	volumenE = calcularVolumen();

	// 3. Inicia la estructura de normales sobre cada cara
	normalesPorCara = (float*)malloc(nCeldas*3*sizeof(float));
	if (normalesPorCara == NULL) {

		_DEBUG("Error allocating normalesPorCara");
		exit(-1);
	}

	for(int i = 0 ; i<nCeldas ; i++)
	{
		float normales[3];

		darNormalCara(i, normales);

		NORMALESPORCARA(i, 0) = normales[0];
		NORMALESPORCARA(i, 1) = normales[1];
		NORMALESPORCARA(i, 2) = normales[2];
	}


	// 4. Construye la estructura de caras por Nodo
	carasPorNodo = (float*)malloc(nNodos*7*sizeof(float));
	if (carasPorNodo == NULL) {

		_DEBUG("Error allocating carasPorNodo");
		exit(-1);
	}

	// Construir la estructura que contiene las caras que comparten cada nodo
	// Recorrer nodos
	int contador = 0;
	int seguidor = 1;
	for(int i=0; i<nNodos ; i++ )
	{
		//Recorrer las caras buscando el nodo i
		for(int j = 0 ; (j < nCeldas); j++)
		{
			if(contieneNodo(i,j))
			{
				//printf("i: %d, j: %d, c: %d, s: %d\n", i, j, contador, seguidor);
				// Agregar la cara a la lista de caras por nodo
				CARASPORNODO(i, seguidor) = j;

				contador++;
				seguidor++;
			}
		}
		// Escribe el numero de caras que tiene cada nodo
		CARASPORNODO(i, 0) = contador;
		contador = 0;
		seguidor = 1;
	}

	// 5. Iniciar vectores normales por cada nodo como el promedio de caras
	normalesPorNodo = (float*)malloc(nNodos*3*sizeof(float));
	if (normalesPorNodo == NULL) {

		_DEBUG("Error allocating normalesPorNodo");
		exit(-1);
	}

	for(int i = 0 ; i<nNodos ; i++)
	{

		float normales[3];

		darNormalPromedio(i, normales);

		NORMALESPORNODO(i, 0) = normales[0];
		NORMALESPORNODO(i, 1) = normales[1];
		NORMALESPORNODO(i, 2) = normales[2];
	}

	// 6. Inicia la estructura de angulos sobre cada nodo
	angulosPorNodo = (float*)malloc(nNodos*7*sizeof(float));
	if (angulosPorNodo == NULL) {

		_DEBUG("Error allocating angulosPorNodo");
		exit(-1);
	}
	for(int i = 0 ; i<nNodos ; i++)
	{
		float angulos[7];

		calcularAngulosPorNodo(i, angulos);

		ANGULOSPORNODO(i, 0) = angulos[0];
		ANGULOSPORNODO(i, 1) = angulos[1];
		ANGULOSPORNODO(i, 2) = angulos[2];
		ANGULOSPORNODO(i, 3) = angulos[3];
		ANGULOSPORNODO(i, 4) = angulos[4];
		ANGULOSPORNODO(i, 5) = angulos[5];
		ANGULOSPORNODO(i, 6) = angulos[6];

	}

	// 9. 10. Inicia la estructura de valores para Laplace beltrami sobre cada nodo
	laplaceKg = new float[nNodos];
	laplaceKh = new float[nNodos];
	for(int i = 0; i < nNodos ; i++)
	{
		laplaceKg[i] = darLaplaceKg(i);
		laplaceKh[i] = darLaplaceKh(i);
	}

	// 11. Iniciar la fuerza en cada nodo
	fuerza = (float*)malloc(nNodos*3*sizeof(float));
	if (fuerza == NULL) {

		_DEBUG("Error allocating fuerza");
		exit(-1);
	}
	memset(fuerza, 0, nNodos*3*sizeof(float));

	// 12. Iniciar la velocidad en cada nodo
	velocidad = (float*)malloc(nNodos*3*sizeof(float));
	if (velocidad == NULL) {

		_DEBUG("Error allocating velocidad");
		exit(-1);
	}
	memset(velocidad, 0, nNodos*3*sizeof(float));


	velocidad2 = (float*)malloc(nNodos*3*sizeof(float));
	if (velocidad2 == NULL) {

		_DEBUG("Error allocating velocidad2");
		exit(-1);
	}
	memset(velocidad2, 0, nNodos*3*sizeof(float));


	// Iniciar la estructura para cambio de nodo
	area = new float[nCeldas];

	// Encontrar los nodos problema 5 caras
	encontrarNodosProblema();
}


// Retorna el valor de Laplace por cada nodo
float mesh::darLaplaceKg(int nodo)
{
	// Implementar operador laplaciano
	//1. Encontrar todos los vecinos del nodo
	int vecinos[7];
	float laplace=0.0;
	darNodosVecinos(nodo, vecinos);
	// 2. Por cada vecino encontrar las dos caras
	for(int i = 1; i<=vecinos[0];i++)
	{
		int pj = vecinos[i];
		int pj1, pj2;
		int caras[2];
		darCarasSegunDosNodos(nodo,pj,caras);

		// 3. Buscar los nodos Pj1 y Pj2
		// Busca a pj1
		int nodoA = FACES(caras[0], 0);
		int nodoB = FACES(caras[0], 1);
		int nodoC = FACES(caras[0], 2);
		if((nodoA != nodo) && (nodoA != pj))
		{
			pj1 = nodoA;
		}else if((nodoB != nodo) && (nodoB != pj))
		{
			pj1 = nodoB;
		}else if((nodoC != nodo) && (nodoC != pj))
		{
			pj1 = nodoC;
		}

		// Busca a pj2
		nodoA = FACES(caras[1], 0);
		nodoB = FACES(caras[1], 1);
		nodoC = FACES(caras[1], 2);
		if((nodoA != nodo) && (nodoA != pj))
		{
			pj2 = nodoA;
		}else if((nodoB != nodo) && (nodoB != pj))
		{
			pj2 = nodoB;
		}else if((nodoC != nodo) && (nodoC != pj))
		{
			pj2 = nodoC;
		}

		// Encuentra posiciones y curvaturas en los nodos nodo, pj, pj1 y pj2
		float posNodo[3], posPj[3], posPj1[3], posPj2[3];
		float kgNodo, kgPj;
		float pj1j[3], pj2j[3], pjnodo[3],fji;
		float d1, d2, d3;

		darPosNodo(nodo, posNodo);
		darPosNodo(pj, posPj);
		darPosNodo(pj1, posPj1);
		darPosNodo(pj2, posPj2);
		kgNodo = darKgPorNodo(nodo);
		kgPj = darKgPorNodo(pj);

		//Calcular vectores y magnitudes
		pj1j[0] = posPj1[0] - posPj[0];
		pj1j[1] = posPj1[1] - posPj[1];
		pj1j[2] = posPj1[2] - posPj[2];

		d1 = norm(pj1j);
		d2 = norm(pj2j);
		d3 = norm(pjnodo);
		fji = kgPj - kgNodo;
		laplace += ((d1+d2)*fji)/(2.0*d3);
	}
	return laplace/darAreaVoronoiPorNodo(nodo);
}

// Retorna el valor de Laplace por cada nodo
float mesh::darLaplaceKh(int nodo)
{
	// Implementar operador laplaciano
	//1. Encontrar todos los vecinos del nodo
	int vecinos[7];
	float laplace=0.0;
	darNodosVecinos(nodo, vecinos);
	// 2. Por cada vecino encontrar las dos caras
	if(~(CARASPORNODO(nodo, 0)<6))
	{
		for(int i = 1; i<=vecinos[0];i++)
		{
			int pj = vecinos[i];
			int pj1, pj2;
			int caras[2];
			darCarasSegunDosNodos(nodo,pj,caras);

			// 3. Buscar los nodos Pj1 y Pj2
			// Busca a pj1
			int nodoA = FACES(caras[0], 0);
			int nodoB = FACES(caras[0], 1);
			int nodoC = FACES(caras[0], 2);
			if((nodoA != nodo) && (nodoA != pj))
			{
				pj1 = nodoA;
			}else if((nodoB != nodo) && (nodoB != pj))
			{
				pj1 = nodoB;
			}else if((nodoC != nodo) && (nodoC != pj))
			{
				pj1 = nodoC;
			}

			// Busca a pj2
			nodoA = FACES(caras[1], 0);
			nodoB = FACES(caras[1], 1);
			nodoC = FACES(caras[1], 2);
			if((nodoA != nodo) && (nodoA != pj))
			{
				pj2 = nodoA;
			}else if((nodoB != nodo) && (nodoB != pj))
			{
				pj2 = nodoB;
			}else if((nodoC != nodo) && (nodoC != pj))
			{
				pj2 = nodoC;
			}

			// Encuentra posiciones y curvaturas en los nodos nodo, pj, pj1 y pj2
			float posNodo[3], posPj[3], posPj1[3], posPj2[3];
			float khNodo, khPj;
			float pj1j[3], pj2j[3], pjnodo[3],fji;
			float d1, d2, d3;

			darPosNodo(nodo, posNodo);
			darPosNodo(pj, posPj);
			darPosNodo(pj1, posPj1);
			darPosNodo(pj2, posPj2);
			khNodo = darKhPorNodo(nodo);
			khPj = darKhPorNodo(pj);

			//Calcular vectores y magnitudes
			pj1j[0] = posPj1[0] - posPj[0];
			pj1j[1] = posPj1[1] - posPj[1];
			pj1j[2] = posPj1[2] - posPj[2];

			pj2j[0] = posPj2[0] - posPj[0];
			pj2j[1] = posPj2[1] - posPj[1];
			pj2j[2] = posPj2[2] - posPj[2];

			pjnodo[0] = posNodo[0] - posPj[0];
			pjnodo[1] = posNodo[1] - posPj[1];
			pjnodo[2] = posNodo[2] - posPj[2];

			d1 = norm(pj1j);
			d2 = norm(pj2j);
			d3 = norm(pjnodo);
			fji = khPj - khNodo;
			laplace += ((d1+d2)*fji)/(2.0*d3);
		}
		float areaB = darAreaAlrededorPorNodo(nodo);
		laplace = laplace/areaB;
	}else
	{
		laplace = darLaplaceKhPromedioPorNodo(nodo);
	}
	return laplace;
}

// Dar Kg como el promedio de los nodos vecinos
float mesh::darKgPromedioPorNodo(int nodo)
{
	int vecinos[7];
	float kg = 0.0;
	darNodosVecinos(nodo, vecinos);
	for(int i=1;i<vecinos[0];i++)
	{
		kg = kg + darKgPorNodo(vecinos[i]);
	}
	return kg/vecinos[0];
}

// Dar Kh como el promedio de los nodos vecinos
float mesh::darKhPromedioPorNodo(int nodo)
{
	int vecinos[7];
	float kh = 0.0;
	darNodosVecinos(nodo, vecinos);
	for(int i=1;i<=vecinos[0];i++)
	{
		kh = kh + darKhPorNodo(vecinos[i]);
	}
	return kh/vecinos[0];
}

// Dar Kh como el promedio de los nodos vecinos
float mesh::darLaplaceKgPromedioPorNodo(int nodo)
{
	int vecinos[7];
	float lkg = 0.0;
	darNodosVecinos(nodo, vecinos);
	for(int i=1;i<=vecinos[0];i++)
	{
		lkg = lkg + darLaplaceKg(vecinos[i]);
	}
	return lkg/vecinos[0];
}

// Dar Kh como el promedio de los nodos vecinos
float mesh::darLaplaceKhPromedioPorNodo(int nodo)
{
	int vecinos[7];
	float lkh = 0.0;
	darNodosVecinos(nodo, vecinos);
	for(int i=1;i<=vecinos[0];i++)
	{
		lkh = lkh + darLaplaceKh(vecinos[i]);
	}
	return lkh/vecinos[0];
}


// Retorna la curvatura Gaussiana por cada nodo
float mesh::darKgPorNodo(int nodo)
{
	float angulos[7];
	float gauss = 0.0;
	float area = darAreaVoronoiPorNodo(nodo);
	float ang = calcularAngulosPorNodo(nodo, angulos);
	gauss = ang/area;
	return gauss;
}

// Retorna la curvatura media por cada nodo
float mesh::darKhPorNodo(int nodo)
{

	float kh=0;
	if(CARASPORNODO(nodo, 0)<6)
	{
		kh = darKhPromedioPorNodo(nodo);
	}else{
		float K[3] = {0.0,0.0,0.0}; // Mean curvature normal operator
		int vecinos[7];
		float xi[3];
		darPosNodo(nodo, xi);
		darNodosVecinos(nodo, vecinos);
		// Por cada nodo vecino debe encontrar la suma de la cotangente y las direcciones
		for(int i=1; i<=vecinos[0];i++)
		{
			float xj[3], d[3];
			darPosNodo(vecinos[i], xj);
			d[0] = xi[0] - xj[0];
			d[1] = xi[1] - xj[1];
			d[2] = xi[2] - xj[2];
			float magd = norm(d);
			float sumCot = darAreaVoronoiParcial(nodo, vecinos[i]);
			// Multiplicar por ocho para ajustar la suma completa
			K[0] += (sumCot*(d[0])*(8.0))/pow(magd,2);
			K[1] += (sumCot*(d[1])*(8.0))/pow(magd,2);
			K[2] += (sumCot*(d[2])*(8.0))/pow(magd,2);
		}

		float areaM = darAreaVoronoiPorNodo(nodo);

		// Calcular la magnitud del operador K
		kh = sqrt(pow(K[0],2) + pow(K[1],2) + pow(K[2],2));
		kh = kh/(2.*areaM);
	}
	// Build normal vectors
	//normalesPorNodo[nodo][0] = K[0]/kh;
	//normalesPorNodo[nodo][1] = K[1]/kh;
	//normalesPorNodo[nodo][2] = K[2]/kh;
	return kh;
}

// Retorna la suma de los angulos formados por las caras alrededor del nodo i
float mesh::calcularAngulosPorNodo(int nodo, float angulos[7])
{
	// Obtiene la cantidad de caras que comparten el nodo
	int numeroCaras = CARASPORNODO(nodo, 0);
	float sumTheta_i = 0.0;
	float pi[3];
	darPosNodo(nodo, pi);

	for(int i = 1; i<=numeroCaras; i++)
	{
		int nCara = CARASPORNODO(nodo, i); // iesima cara

		// Obtiene los tres vertices de cada cara
		int n1, n2, n3;
		n1 = FACES(nCara, 0);
		n2 = FACES(nCara, 1);
		n3 = FACES(nCara, 2);

		// Calcular el angulo de la i-esima cara
		float pA[3];
		float pB[3];

		// Casos dependiendo cual sea el vertice que se repite
		if(n1==nodo)
		{
			darPosNodo(n2, pA);
			darPosNodo(n3, pB);
		}else{
			if(n2==nodo)
			{
				darPosNodo(n1, pA);
				darPosNodo(n3, pB);
			}else{
				if(n3==nodo)
				{
					darPosNodo(n1, pA);
					darPosNodo(n2, pB);
				}
			}
		}

		// Construir los vectores desde i hasta los demas puntos de la cara
		float vec1[3], vec2[3];
		vec1[0] = pA[0] - pi[0];
		vec1[1] = pA[1] - pi[1];
		vec1[2] = pA[2] - pi[2];

		vec2[0] = pB[0] - pi[0];
		vec2[1] = pB[1] - pi[1];
		vec2[2] = pB[2] - pi[2];

		// Realizar producto interior entre vec1 y vec2
		float AdotB = (vec1[0]*vec2[0]) + (vec1[1]*vec2[1]) + (vec1[2]*vec2[2]);
		float magA = sqrt((vec1[0]*vec1[0]) + (vec1[1]*vec1[1]) + (vec1[2]*vec1[2]));
		float magB = sqrt((vec2[0]*vec2[0]) + (vec2[1]*vec2[1]) + (vec2[2]*vec2[2]));
		float theta_i = acos(AdotB/(magA*magB));
		angulos[i] = theta_i;
		sumTheta_i += theta_i;
	}
	// Indica cuantos ángulos estan formados alrededor del nodo i
	angulos[0] = numeroCaras;
	float retorno = ((2*3.141617) - sumTheta_i);
	return retorno;
}


// Encuentra el area de voronoi parcial utilizando dos nodos
// Retorna (1/8) de la suma de cotangentes
float mesh::darAreaVoronoiParcial(int nodoA, int nodoB)
{
	float area= 0.0;
	int caras[2];
	darCarasSegunDosNodos(nodoA, nodoB, caras);

	// Calcular el area de voronoi parcial utilizando las dos caras
	// el angulo alpha corresponde a la cara 2

	// Encontrar el tercer nodo en la cara 1
	int nodoC;
	for(int k = 0; k < 3; k++)
	{
		if((FACES(caras[1], k) != nodoA) && (FACES(caras[1], k) != nodoB))
		{
			nodoC = FACES(caras[1], k);
		}
	}

	// Trazar los vectores entre nodos A, B y C
	float vec1[3], vec2[3], vec3[3];

	vec1[0] = VERTEX(nodoC, 0) - VERTEX(nodoB, 0);
	vec1[1] = VERTEX(nodoC, 1) - VERTEX(nodoB, 1);
	vec1[2] = VERTEX(nodoC, 2) - VERTEX(nodoB, 2);

	vec2[0] = VERTEX(nodoC, 0) - VERTEX(nodoA, 0);
	vec2[1] = VERTEX(nodoC, 1) - VERTEX(nodoA, 1);
	vec2[2] = VERTEX(nodoC, 2) - VERTEX(nodoA, 2);

	float beta = 0.0;
	float alpha = 0.0;
	float mag1, mag2, mag3;


	// Producto punto entre los vectores
	mag1 = sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2));
	mag2 = sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2));
	float v1Dotv2 = (vec1[0]*vec2[0]) + (vec1[1]*vec2[1]) + (vec1[2]*vec2[2]);

	//Producto cross entre los vectores
	cross(vec1, vec2, vec3[0],vec3[1],vec3[2] );
	mag3 = sqrt(pow(vec3[0],2) + pow(vec3[1],2) + pow(vec3[2],2));

	beta = atan2(mag3, v1Dotv2);
	//beta = acos((v1Dotv2)/(mag1*mag2));

	// Encontrar el tercer nodo en la cara 2
	for(int k = 0; k < 3; k++)
	{
		if((FACES(caras[0], k) != nodoA) && (FACES(caras[0], k) != nodoB))
		{
			nodoC = FACES(caras[0], k);
		}
	}

	// Trazar los vectores entre nodos A, B y C
	vec1[0] = VERTEX(nodoC, 0) - VERTEX(nodoB, 0);
	vec1[1] = VERTEX(nodoC, 1) - VERTEX(nodoB, 1);
	vec1[2] = VERTEX(nodoC, 2) - VERTEX(nodoB, 2);

	vec2[0] = VERTEX(nodoC, 0) - VERTEX(nodoA, 0);
	vec2[1] = VERTEX(nodoC, 1) - VERTEX(nodoA, 1);
	vec2[2] = VERTEX(nodoC, 2) - VERTEX(nodoA, 2);

	mag1 = sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2));
	mag2 = sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2));
	v1Dotv2 = (vec1[0]*vec2[0]) + (vec1[1]*vec2[1]) + (vec1[2]*vec2[2]);

	//Producto cross entre los vectores
	cross(vec1, vec2, vec3[0],vec3[1],vec3[2] );
	mag3 = sqrt(pow(vec3[0],2) + pow(vec3[1],2) + pow(vec3[2],2));

	alpha = atan2(mag3, v1Dotv2);
	//alpha = acos((v1Dotv2)/(mag1*mag2));

	float d[3];
	d[0] = VERTEX(nodoA, 0) - VERTEX(nodoB, 0);
	d[1] = VERTEX(nodoA, 1) - VERTEX(nodoB, 1);
	d[2] = VERTEX(nodoA, 2) - VERTEX(nodoB, 2);

	float magd = 0.0;
	magd = sqrt((pow(d[0],2) + pow(d[1],2) + pow(d[2],2)));
	area = (1./8.)*(pow(magd,2))*((1./tan(alpha)) + (1./tan(beta)));
	return area;
}

float mesh::darAreaVoronoiPorNodo(int nodo)
{
	float area = 0.0;
	// Encontrar los J nodos que rodean al nodo i
	int vecinos[7];
	darNodosVecinos(nodo, vecinos);
	// Suma de cotangentes
	for(int i = 0; i < vecinos[0] ; i++)
	{
		// Encontrar los angulos alpha y beta
		area += darAreaVoronoiParcial(nodo, vecinos[i+1]);
	}
	return area;
}

void mesh::darCarasPorNodo(int nodo, int caras[7])
{
	for(int i = 0; i<7; i++)
	{
		caras[i] = CARASPORNODO(nodo, i);
	}
}

// Dar nodos que forman la cara
void mesh::darNodosPorElemento(int cara, int nodos[3])
{
	nodos[0] = FACES(cara, 0);
	nodos[1] = FACES(cara, 1);
	nodos[2] = FACES(cara, 2);
}

// Dar las dos caras que tienen dos nodos en comun
void mesh::darCarasSegunDosNodos(int nodoA, int nodoB, int caras[2])
{
	// Buscar en todas las caras si estan presentes los dos nodos
	int lugar = 0;
	for(int cara = 0; cara < nCeldas ; cara++)
	{
		if(contieneNodo(nodoA, cara))
		{
			if(contieneNodo(nodoB, cara))
			{
				caras[lugar] = cara;
				lugar++;
			}
		}
	}
}


// Retorna el area alrededor de un nodo utilizando el baricentro de cada cara
float mesh::darAreaBaricentricaPorNodo(int nodo)
{
	int caras = CARASPORNODO(nodo, 0);
	float area = 0.0;
	// Sumar el area de todas las caras y dividirlo por 3
	for(int i = 1; i <= caras; i++)
	{
		area += darAreaElemento(CARASPORNODO(nodo, i));
	}
	return area/3.;
}

// Encuentra todos los nodos vecinos en el primer anillo dado un nodo i
void mesh::darNodosVecinos(int nodo, int vecinos[7])
{
	//Inicializar vecinos en -1, la posición 0 indica el numero de vecinos
	for(int i = 0 ; i<7 ; i++)
	{
		vecinos[i] = -1;
	}
	// Encuentra la cantidad de caras que comparten el nodo
	int caras = CARASPORNODO(nodo, 0);
	// Crea el arreglo para guardar todos los vertices
	int *nodos = new int[(caras*3)-3];
	// El numero de caras es igual al numero de vecinos del nodo
	vecinos[0] = caras;
	// Recorre la estructura de nodos y selecciona aquellos que no estan en vertices
	int seguidor = 1;
	int lugar = 0;
	bool repetido = false;

	// Inicializar la estructura nodos con los nodos de cada cara
	for(int i=1; i<caras; i++)
	{
		// Por cada cara debe agregar tres nodos a la estructura
		int nodoA, nodoB, nodoC;
		int cara = CARASPORNODO(nodo, i);
		nodoA = FACES(cara, 0);
		nodoB = FACES(cara, 1);
		nodoC = FACES(cara, 2);
		nodos[lugar] = nodoA;
		lugar++;
		nodos[lugar] = nodoB;
		lugar++;
		nodos[lugar] = nodoC;
		lugar++;
	}

	for(int i = 0; i < (caras*3)-3 ; i++){
		// Recorre el arreglo de vecinos para buscar el nodo i
		if(nodos[i] == nodo){
			repetido = true;
		}
		else {
			for(int j=1;j<seguidor;j++) {
				//Compara cada vecino con nodo i
				if((vecinos[j] == nodos[i])) {
					repetido = true;
				}
			}
		}

		if(!repetido){
			vecinos[seguidor] = nodos[i];
			seguidor++;
		}
		repetido = false;
	}
	delete nodos;
}

// Entrega el vector normal dado un nodo i
void mesh::darNormalPromedio(int nodo, float normal[3])
{
	// Encuentra el numero de caras que comparten el nodo
	int numeroCaras = CARASPORNODO(nodo, 0);

	float pX = 0.0;
	float pY = 0.0;
	float pZ = 0.0;
	int nCara = 0;
	for(int c = 1 ; c<=numeroCaras ; c++)
	{
		nCara = CARASPORNODO(nodo, c);
		pX+=NORMALESPORCARA(nCara, 0);
		pY+=NORMALESPORCARA(nCara, 1);
		pZ+=NORMALESPORCARA(nCara, 2);
	}
	normal[0] = pX/numeroCaras;
	normal[1] = pY/numeroCaras;
	normal[2] = pZ/numeroCaras;
}

// Retorna true o false si el nodo esta en la cara
bool mesh::contieneNodo(int nodo, int cara)
{
	bool flag = false;
	for(int p = 0; p < 3; p++)
	{
		if(FACES(cara, p)==nodo){
			flag = true;}
	}
	return flag;
}

void mesh::darNormalCara(int j, float normal[3]){

	// Encontrar los tres nodos que componen la cara
	float nodoA[3], nodoB[3], nodoC[3], a[3], b[3];
	darPosNodo(FACES(j, 0), nodoA);
	darPosNodo(FACES(j, 1), nodoB);
	darPosNodo(FACES(j, 2), nodoC);

	// Construir vector 1
	a[0] = nodoB[0] - nodoA[0];
	a[1] = nodoB[1] - nodoA[1];
	a[2] = nodoB[2] - nodoA[2];

	// Construir vector 2
	b[0] = nodoC[0] - nodoA[0];
	b[1] = nodoC[1] - nodoA[1];
	b[2] = nodoC[2] - nodoA[2];

	// Utilizar el producto cross para encontrar el vector normal
	normal[0] = a[1]*b[2] - a[2]*b[1];
	normal[1] = b[0]*a[2] - a[0]*b[2];
	normal[2] = a[0]*b[1] - b[0]*a[1];

	// Entrega el vector normalizado
	float mag = sqrt(pow(normal[0],2) + pow(normal[1],2)+ pow(normal[2],2));
	normal[0] = (normal[0]/mag)*(-1);
	normal[1] = (normal[1]/mag)*(-1);
	normal[2] = (normal[2]/mag)*(-1);
}

// 1. Actualiza normales por cara
// 2. Actualiza normales por nodo
// 3. Actualiza curvaturas por nodo

void mesh::actualizarGeometria()
{
	// 1. Actualiza normales por cara
	for(int i = 0; i<nCeldas ; i++)
	{
		float normales[3];

		darNormalCara(i, normales);

		NORMALESPORCARA(i, 0) = normales[0];
		NORMALESPORCARA(i, 1) = normales[1];
		NORMALESPORCARA(i, 2) = normales[2];

	}
}

void mesh::setVelocidad(int n, float ux, float uy, float uz)
{
	VELOCIDAD2(n, 0) = VELOCIDAD(n, 0);
	VELOCIDAD2(n, 1) = VELOCIDAD(n, 1);
	VELOCIDAD2(n, 2) = VELOCIDAD(n, 2);

	VELOCIDAD(n, 0) = ux;
	VELOCIDAD(n, 1) = uy;
	VELOCIDAD(n, 2) = uz;
}


// Alamacena las tres componentes de laposicion en pos
void mesh::darPosNodo(int n, float pos[3])
{
	pos[0] = VERTEX(n, 0);
	pos[1] = VERTEX(n, 1);
	pos[2] = VERTEX(n, 2);
}

void mesh::darFuerzaNodo(int n, float f[3])
{
	f[0] = FUERZA_MESH(n, 0);
	f[1] = FUERZA_MESH(n, 1);
	f[2] = FUERZA_MESH(n, 2);
}


void mesh::moverNodos(float dt, float dx)
{
	for(int u=0;u<nNodos;u++)
	{
		VERTEX(u, 0) += ((3/2)*(VELOCIDAD(u, 0)) - (1./2.0)*(VELOCIDAD2(u, 0)))*dt;
		VERTEX(u, 1) += ((3/2)*(VELOCIDAD(u, 1)) - (1./2.0)*(VELOCIDAD2(u, 1)))*dt;
		VERTEX(u, 2) += ((3/2)*(VELOCIDAD(u, 2)) - (1./2.0)*(VELOCIDAD2(u, 2)))*dt;
	}
}

void mesh::calcularFuerzasFEM(mesh referencia, float ks)
{
	int a,b,c;
	float va[3], vb[3], vc[3], vA[3], vB[3], vC[3];
	// Crear estructuras para pasar al algoritmo de fuerzas
	float ref[3][3], def[3][3], fuerzas[3][3];

	// Obtener indices de los nodos de cada elemento
	for(int e=0;e<nCeldas;e++)
	{
		a = FACES(e, 0);
		b = FACES(e, 1);
		c = FACES(e, 2);

		// Obtener las posiciones de cada vertice no deformado
		referencia.darPosNodo(a,va);
		referencia.darPosNodo(b,vb);
		referencia.darPosNodo(c,vc);

		// Obtener las posiciones de cada vertice deformado
		darPosNodo(a,vA);
		darPosNodo(b,vB);
		darPosNodo(c,vC);

		ref[0][0] = va[0];
		ref[0][1] = va[1];
		ref[0][2] = va[2];

		ref[1][0] = vb[0];
		ref[1][1] = vb[1];
		ref[1][2] = vb[2];

		ref[2][0] = vc[0];
		ref[2][1] = vc[1];
		ref[2][2] = vc[2];

		def[0][0] = vA[0];
		def[0][1] = vA[1];
		def[0][2] = vA[2];

		def[1][0] = vB[0];
		def[1][1] = vB[1];
		def[1][2] = vB[2];

		def[2][0] = vC[0];
		def[2][1] = vC[1];
		def[2][2] = vC[2];

		rotacion(ref, def, fuerzas, ks);

		//Agregar cada fuerza a cada nodo
		FUERZA_MESH(a, 0) = FUERZA_MESH(a, 0)-fuerzas[0][0];
		FUERZA_MESH(a, 1) = FUERZA_MESH(a, 1)-fuerzas[0][1];
		FUERZA_MESH(a, 2) = FUERZA_MESH(a, 2)-fuerzas[0][2];

		FUERZA_MESH(b, 0) = FUERZA_MESH(b, 0)-fuerzas[1][0];
		FUERZA_MESH(b, 1) = FUERZA_MESH(b, 1)-fuerzas[1][1];
		FUERZA_MESH(b, 2) = FUERZA_MESH(b, 2)-fuerzas[1][2];

		FUERZA_MESH(c, 0) = FUERZA_MESH(c, 0)-fuerzas[2][0];
		FUERZA_MESH(c, 1) = FUERZA_MESH(c, 1)-fuerzas[2][1];
		FUERZA_MESH(c, 2) = FUERZA_MESH(c, 2)-fuerzas[2][2];
	}
}


void mesh::calcularFuerzasVolumen(float v0, float ka)
{
	float mag = ka*(v0 - calcularVolumen());

	float normal[3];
	for(int i = 0; i < nNodos ; i++)
	{
		normal[0] = NORMALESPORNODO(i, 0);
		normal[1] = NORMALESPORNODO(i, 1);
		normal[2] = NORMALESPORNODO(i, 2);
		FUERZA_MESH(i, 0) = FUERZA_MESH(i, 0)+normal[0]*mag;
		FUERZA_MESH(i, 1) = FUERZA_MESH(i, 1)+normal[1]*mag;
		FUERZA_MESH(i, 2) = FUERZA_MESH(i, 2)+normal[2]*mag;
	}
}
void mesh::calcularFuerzasHelfrich(float kb)

{
	for(int i = 0; i < nNodos; i++)
	{
		FUERZA_MESH(i, 0) = 0.0;
		FUERZA_MESH(i, 1) = 0.0;
		FUERZA_MESH(i, 2) = 0.0;
		/*float c0 = 1.0e-9;
		float kh = darKhPorNodo(i);
		float kg = darKgPorNodo(i);
		float lkh = darLaplaceKh(i);
		float mag = 0.0;
		mag = (kb*(((2*kh) + c0)*(2*(kh*kh)-(2*kg)-(kh*c0)))) + (2*kb*lkh);
		FUERZA_MESH(i, 0) = FUERZA_MESH(i, 0)-mag*normalesPorNodo[i][0];
		FUERZA_MESH(i, 1) = FUERZA_MESH(i, 1)-mag*normalesPorNodo[i][1];
		FUERZA_MESH(i, 2) = FUERZA_MESH(i, 2)-mag*normalesPorNodo[i][2];*/
	}
}

/**
+	Retorna el volumen del tetrahedro formado por los vertices de cada elemento
 *   y el centro de la esfera.
 *
 */
float mesh::darVolumenElemento(int i)
{
	float a[3], b[3], c[3], d[3], v1[3], v2[3], v3[3], temp[3], V;
	int A, B, C;

	A=FACES(i, 0);
	B=FACES(i, 1);
	C=FACES(i, 2);

	a[0]=VERTEX(A, 0);
	a[1]=VERTEX(A, 1);
	a[2]=VERTEX(A, 2);

	b[0]=VERTEX(B, 0);
	b[1]=VERTEX(B, 1);
	b[2]=VERTEX(B, 2);

	c[0]=VERTEX(C, 0);
	c[1]=VERTEX(C, 1);
	c[2]=VERTEX(C, 2);

	d[0] = cX;
	d[1] = cY;
	d[2] = cZ;

	v1[0] = a[0] - d[0];
	v1[1] = a[1] - d[1];
	v1[2] = a[2] - d[2];

	v2[0] = b[0] - d[0];
	v2[1] = b[1] - d[1];
	v2[2] = b[2] - d[2];

	v3[0] = c[0] - d[0];
	v3[1] = c[1] - d[1];
	v3[2] = c[2] - d[2];

	cross(v2, v3, temp[0], temp[1], temp[2]);
	V = (temp[0]*v1[0] + temp[1]*v1[1] + temp[2]*v1[2])/6.0;
	return fabs(V);
}

float mesh::calcularVolumen()

{
	volumenE = 0.0;
	for(int i = 0; i<nCeldas ; i++)
	{
		volumenE += darVolumenElemento(i);
	}
	return volumenE;
}

float mesh::darAreaElemento(int i)
{
	float a[3], b[3], c[3], v1[3], v2[3], temp[3];
	int A, B, C;
	A=FACES(i, 0);
	B=FACES(i, 1);
	C=FACES(i, 2);

	a[0]=VERTEX(A, 0);
	a[1]=VERTEX(A, 1);
	a[2]=VERTEX(A, 2);

	b[0]=VERTEX(B, 0);
	b[1]=VERTEX(B, 1);
	b[2]=VERTEX(B, 2);

	c[0]=VERTEX(C, 0);
	c[1]=VERTEX(C, 1);
	c[2]=VERTEX(C, 2);

	v1[0] = b[0] - a[0];
	v1[1] = b[1] - a[1];
	v1[2] = b[2] - a[2];

	v2[0] = c[0] - a[0];
	v2[1] = c[1] - a[1];
	v2[2] = c[2] - a[2];

	cross(v1, v2, temp[0], temp[1], temp[2]);
	return norm(temp)/2.0;
}


float mesh::calcularAreaSuperficial()
{
	float area = 0.0;
	for(int i = 0 ; i < nCeldas; i++)
	{
		area += darAreaElemento(i);
	}
	return area;
}



void mesh::calcularCambioArea(mesh ref)
{
	for(int i = 0 ; i < nCeldas; i++)
	{
		area[i] = darAreaElemento(i)-ref.darAreaElemento(i);
	}
}

// Calcula la fuerza neta sobre la membrana
void mesh::calcularFuerzaNeta(float fNeta[3])
{
	float fx=0., fy=0., fz=0., fNodo[3] ={0.,0.,0.};
	for(int u = 0.0; u < nNodos ; u++)
	{
		darFuerzaNodo(u, fNodo);
		fx += fNodo[0];
		fy += fNodo[1];
		fz += fNodo[2];
	}
	fNeta[0] = fx;
	fNeta[1] = fy;
	fNeta[2] = fz;
	fx = 0.0;
	fy = 0.0;
	fz = 0.0;
}


float mesh::calcularEnergia()
{
	return 0.0;
}

void mesh::calcularMomentoNeto(float fMomento[3])
{

}

void mesh::actualizarNodos(float **nodos)
{
	float x, y, z;
	for(int u = 0 ; u< nNodos ; u++)
	{
		x = nodos[u][0];
		y = nodos[u][1];
		z = nodos[u][2];

		VERTEX(u, 0) = x;
		VERTEX(u, 1) = y;
		VERTEX(u, 2) = z;
	}
}


void mesh::encontrarNodosProblema()
{
	nodosProblema = new int[12];
	int seguidor = 0;
	for(int i = 0 ; i<nNodos;i++)
	{
		if(CARASPORNODO(i, 0)<6)
		{
			nodosProblema[seguidor] = i;
			seguidor++;
		}
	}
}


bool mesh::esNodoProblema(int nodo)
{
	bool flag = false;
	for(int i = 0; i<12; i++)
	{
		if(nodo == nodosProblema[i])
		{
			flag = true;
		}
	}
	return flag;
}

float mesh::darAreaAlrededorPorNodo(int nodo)
{
	float area = 0.0;
	for(int i = 1; i <= CARASPORNODO(nodo, 0); i++)
	{
		area += darAreaElemento(CARASPORNODO(nodo, i));
	}
	return area;
}

float mesh::darK1PorNodo(int nodo)
{
	float kh = darKhPorNodo(nodo);
	float kg = darKgPorNodo(nodo);
	return (kh + sqrt(pow(kh,2) - kg));
}

float mesh::darK2PorNodo(int nodo)
{
	float kh = darKhPorNodo(nodo);
	float kg = darKgPorNodo(nodo);
	return (kh - sqrt(pow(kh,2) - kg));
}
