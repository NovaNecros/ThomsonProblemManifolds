#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

//individuos con fitness
typedef struct
{
	double * parametros;
	//parametros[0] = tamPob;
	//parametros[1] = probaCruza;
	//parametros[2] = probaMut;
	double fitness;
}vector_DE;

//Metodo de Fisher-Yates para permutar un arreglo
void shuffle(int * L, int n) 
{
	for(int i=n-1; i>0; --i) 
	{
		int j = rand() % (i+1);
		int aux = L[i];
		L[i] = L[j];
		L[j] = aux;
	}
}

//funcion objetivo
void funcionFitness(int numParam, vector_DE * vec)
{
		char * args[numParam+2];
	
    args[0] = "ThomsonGenetico.c";
		for(int i=1; i<=numParam; ++i)
		{
			sprintf(args[i], "%lf", vec->parametros[i]);
		}
		args[numParam+1] = NULL;

    execvp(args[0], args);
    //perror("execvp");

    FILE * energyFile = fopen("energyval.txt", "r");
    fscanf(energyFile, "%lf", &vec->fitness);
    fclose(energyFile);
    free(energyFile);
}

//Genera poblacion inicial
void poblacionInicial(int tamPob, int numParam, double * inf, double * sup, vector_DE  * pob)
{
	for(int i=0; i<tamPob; ++i)
	{
		for(int j=0; j<numParam; ++j)
		{
			pob[i].parametros[j] = (double)rand()/RAND_MAX*(sup[j]-inf[j])+inf[j];
		}
		funcionFitness(numParam, &pob[i]); 
	}
}

//cruza de 3 
void mutation(int tamPob, int numParam, double tasaMut, double * inf, double * sup, vector_DE * objetivo, vector_DE * ruidosos)
{
	int * rnd = (int *)malloc(sizeof(int)*tamPob);
	int op;
	for(int i=0; i<tamPob; ++i)
	{
		rnd[i] = i;
	}
	
	for(int i=0; i<tamPob; ++i)
	{
		shuffle(rnd, tamPob);
		
		for(int j=0; j<numParam; ++j)
		{
			ruidosos[i].parametros[j] = objetivo[rnd[0]].parametros[j] + tasaMut*(objetivo[rnd[1]].parametros[j]-objetivo[rnd[2]].parametros[j]);
			while(ruidosos[i].parametros[j] < inf[j])
			{
				op = rand() % 3;
				if(op == 0) //Dirichlet
				{
					ruidosos[i].parametros[j] = inf[j];
				}
				else if(op == 1) //Neumann 
				{
					ruidosos[i].parametros[j] = 2.0*inf[j] - ruidosos[i].parametros[j];
				}
				else //Periodicas
				{
					ruidosos[i].parametros[j] = 2.0*sup[j] + ruidosos[i].parametros[j];
				}	
			}
			while(ruidosos[i].parametros[j] > sup[j])
			{
				op = rand() % 3;
				if(op == 0) //Dirichlet
				{
					ruidosos[i].parametros[j] = sup[j];
				}
				else if(op == 1) //Neumann 
				{
					ruidosos[i].parametros[j] = 2.0*sup[j] - ruidosos[i].parametros[j];
				}
				else //Periodicas
				{
					ruidosos[i].parametros[j] = 2.0*inf[j] + ruidosos[i].parametros[j];
				}
			}
		}
	}
	
	free(rnd);
}

//seleccion aleatoria gen por gen
void recombination(int tamPob, int numParam, double tasaCr, vector_DE * objetivo, vector_DE * ruidosos, vector_DE * prueba)
{
	for(int i=0; i<tamPob; ++i)
	{
		for(int j=0; j<numParam; ++j)
		{
			if((double)rand()/RAND_MAX < tasaCr)
			{
				prueba[i].parametros[j] = ruidosos[i].parametros[j];
			}
			else
			{
				prueba[i].parametros[j] = objetivo[i].parametros[j];
			}
		}
		
		funcionFitness(numParam, &prueba[i]); 
	}
}

//determinista
void selection(int tamPob, int numParam, vector_DE * objetivo, vector_DE * prueba)
{
	for(int i=0; i<tamPob; ++i)
	{
		if(prueba[i].fitness < objetivo[i].fitness) //minimizacion
		{
			for(int j=0; j<numParam; ++j)
			{
				objetivo[i].parametros[j] = prueba[i].parametros[j];
			}
			objetivo[i].fitness = prueba[i].fitness;
		}
	}
}

//cuerpo principal evolucion diferencial
void differentialEvolution(int numGen, int tamPob, double tasaCr, double tasaMut, int numParam, double * inf, double * sup, vector_DE * mejorSolucion)
{
	srand(time(NULL));
	
	vector_DE * objetivo = (vector_DE *)malloc(sizeof(vector_DE)*tamPob);
	vector_DE * ruidosos = (vector_DE *)malloc(sizeof(vector_DE)*tamPob);
	vector_DE * prueba = (vector_DE *)malloc(sizeof(vector_DE)*tamPob);
	
	for(int i=0; i<tamPob; ++i)
	{
		objetivo[i].parametros = (double *)malloc(sizeof(double)*numParam);
		ruidosos[i].parametros = (double *)malloc(sizeof(double)*numParam);
		prueba[i].parametros = (double *)malloc(sizeof(double)*numParam);
	}
	
	
	poblacionInicial(tamPob, numParam, inf, sup, objetivo);
	
	
	for(int t=0; t<numGen; ++t)
	{
		//Operadores evolutivos
		mutation(tamPob, numParam, tasaMut, inf, sup, objetivo, ruidosos);
		recombination(tamPob, numParam, tasaCr, objetivo, ruidosos, prueba);
		selection(tamPob, numParam, objetivo, prueba);
	}
	

	for(int i=0; i<numParam; ++i)
	{
		if(mejorSolucion->fitness > objetivo[i].fitness) //minimizacion
		{
			for(int j=0; j<numParam; ++j)
			{
				mejorSolucion->parametros[j] = objetivo[i].parametros[j];	
			}
			mejorSolucion->fitness = objetivo[i].fitness;
		}
	}
	
	free(objetivo);
	free(ruidosos);
	free(prueba);
}

int main() 
{
	//Parametros algoritmo
	int numGeneraciones = 200;
	int tamPoblacion = 20;
	double tasaCruza = 0.5;
	double tasaMutacion = 0.5;

	//Parametros problema
	int numParametros = 2;
	double * limiteInferior = (double *)malloc(sizeof(double)*numParametros);
	double * limiteSuperior = (double *)malloc(sizeof(double)*numParametros);

	limiteInferior[0] = 1.0;
	limiteSuperior[1] = 16.0;
	for(int i=1; i<3; ++i)
	{
		limiteInferior[i] = 0.0;
		limiteSuperior[i] = 1.0;
	}
	
	//Salida
	vector_DE s;
	s.parametros = (double *)malloc(sizeof(double)*numParametros);
	s.fitness = INFINITY; //minimizacion
	
	differentialEvolution(numGeneraciones, tamPoblacion, tasaCruza, tasaMutacion, numParametros, limiteInferior, limiteSuperior, &s);
	
	printf("\nTam Poblacion:\t%d",2*(int)s.parametros[0]);
	printf("\nProba Cruza:\t%f", s.parametros[1]);
	printf("\nProba Mutacion:\t%f", s.parametros[2]);
	printf(")\nEnergia\t%f", s.fitness);

	/*FILE * file = fopen("parametros.txt", "w");
	if(file == NULL)
	{
		printf("Error opening file.\n");
		return 420;
	}
 
	fprintf(file, "%d %lf %lf", 2*(int)s.parametros[0], s.parametros[1], s.parametros[2]);
	fclose(file);*/
	
	free(limiteInferior);
	free(limiteSuperior);
	
  return 69;
}
