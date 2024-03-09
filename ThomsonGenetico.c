#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <unistd.h>

// tipo de dato electron (gen)
typedef struct 
{
  double * surface;
  double * cartesian;
} Electron;

// tipo de dato conifguracion de electrones (individuo)
typedef struct 
{
  Electron * electrons;
  double energy;
} Configuration;

// verifica que los parametros u,v esten en los rangos dados
void normalize(Electron * e, int dimV, double * UMIN, double * UMAX, double * URANGO) 
{
  for(int i=0; i<dimV; ++i)
  {
    while(e->surface[i] < UMIN[i])
	{
      e->surface[i] += URANGO[i];
    }
    while(e->surface[i] > UMAX[i])
	{
      e->surface[i] -= URANGO[i];
    }
  }
}

// covierte puntos en la superficie a puntos en R3
void cartesianCoordinates(Electron * e) 
{
  // Esfera
  /*const double R = 1.0;
  e->cartesian[0] = R*sin(e->surface[0])*cos(e->surface[1]);
  e->cartesian[1] = R*sin(e->surface[0])*sin(e->surface[1]);
  e->cartesian[2] = R*cos(e->surface[0]);*/

  // Cilindro
  /*const double R = 1.0, H = 4.0;
  e->cartesian[0] = R*cos(e->surface[0]);
  e->cartesian[1] = R*sin(e->surface[0]);
  e->cartesian[2] = e->surface[1];*/

  // Toro
  const double A = 4.0, B = 1.0;
  e->cartesian[0] = (A + B*cos(e->surface[0]))*cos(e->surface[1]);
  e->cartesian[1] = (A + B*cos(e->surface[0]))*sin(e->surface[1]);
  e->cartesian[2] = B*sin(e->surface[0]);
}

// funcion objetivo
void configurationEnergy(Configuration * conf, int dimE, int numElectrones) 
{
  conf->energy = 0.0;
  double di, d;

  for(int i=1; i<numElectrones; ++i) 
  {
    for(int j=0; j<i; ++j) 
	{
      d = 0.0;
      for(int k=0; k < dimE; ++k)
	  {
        di = conf->electrons[i].cartesian[k] - conf->electrons[j].cartesian[k];
        d += (di*di);
      }
      if(d == 0.0) 
	  {
        conf->energy = INFINITY; // division entre cero, dos electrones en el mismo punto
        break;
      }

      conf->energy += (1.0 / sqrt(d));
    }
  }
}

// Metodo de Marsaglia para generar aleatorios normales
void normalRand(double * L)
{
  double u1, u2, w;

  do 
  {
    u1 = 2.0*rand()/RAND_MAX - 1.0;
    u2 = 2.0*rand()/RAND_MAX - 1.0;
    w = u1*u1 + u2*u2;
  } while (w >= 1.0 || w <= 0.0);

  L[0] = u1*sqrt(-2.0*log(w)/w);
  L[1] = u2*sqrt(-2.0*log(w)/w);
}

// Metodo de Fisher-Yates para permutar un arreglo
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

//poblacion inicial
void poblacionInicial(Configuration * poblacion, int dimE, int dimV, int tamPob, int numElectrones, double * UMIN, double * URANGO)
{
  for(int i=0; i<tamPob; ++i)
  {
    for(int j=0; j<numElectrones; ++j) 
	{
      for(int k=0; k<dimV; ++k) 
	  {
        poblacion[i].electrons[j].surface[k] = UMIN[k] + ((double)rand()/RAND_MAX)*URANGO[k];
      }
      cartesianCoordinates(&poblacion[i].electrons[j]);
    }
    configurationEnergy(&poblacion[i], dimE, numElectrones);
  }
}

//cruza con media con pesos aleatorios //mejorar
void recombination(Configuration * poblacion, Configuration * descendencia, int dimV, int tamPob, double probaCruza, int numElectrones)
{
  int * ind = (int *)malloc(sizeof(int)*tamPob);
  double peso;

  for(int i=0; i<tamPob; ++i) 
  {
    ind[i] = i;
  }

  shuffle(ind, tamPob);

  for(int i=0; i<tamPob; i+=2)
  {
    for(int j=0; j<numElectrones; ++j)
	{
      if((double)rand() / RAND_MAX < probaCruza)
	  {
        for(int k=0; k<2; ++k)
		{
          for(int p=0; p<dimV; ++p)
		  {
            peso = (double)rand()/RAND_MAX;
            descendencia[i+k].electrons[j].surface[p] = (peso)*poblacion[ind[i]].electrons[j].surface[p] + (1.0-peso)*poblacion[ind[i+1]].electrons[j].surface[p];
          }
        }
      } 
	  else
	  {
        for(int k=0; k<2; ++k)
		{
          for(int p=0; p<dimV; ++p)
		  {
            descendencia[i+k].electrons[j].surface[p] = poblacion[ind[i+k]].electrons[j].surface[p];
          }
        }
      }
    }
  }
  free(ind);
}

//mutacion con perturbaciones gaussianas
void mutation(Configuration * poblacion, int dimV, int tamPob, double probaMutacion, int t, int numGen, int numElectrones, double * UMIN, double * UMAX, double * URANGO)
{
  double gaux = (double) 1.0 / t;
  double * rnd = (double *)malloc(sizeof(double)*2);
  double * stdDev = (double *)malloc(sizeof(double)*dimV);

  for(int i=0; i<dimV; ++i)
  {
    stdDev[i] = URANGO[i] * gaux;
  }

  for(int i=0; i<tamPob; ++i)
  {
    for(int j=0; j<numElectrones; ++j)
	{
      for(int k=0; k<dimV; ++k)
	  {
        if((double)rand()/RAND_MAX < probaMutacion)
		{
          normalRand(rnd);
          poblacion[i].electrons[j].surface[k] = poblacion[i].electrons[j].surface[k] + stdDev[k]*rnd[k];
        }

        normalize(&poblacion[i].electrons[j], dimV, UMIN, UMAX, URANGO);
      }
    }
  }
  free(rnd);
  free(stdDev);
}

// seleccion por ruleta con elitismo
void selection(Configuration * poblacion, Configuration * descendencia, int dimE, int dimV, int tamPob, int tamPobCombinada, int numElectrones) 
{
  double normalizacion = 0.0;
  double * probabilidades = (double *)malloc(sizeof(double)*tamPobCombinada);
  Configuration * pobCombinada = (Configuration *)malloc(sizeof(Configuration)*tamPobCombinada);

  for(int i=0; i<tamPob; ++i)
  {
    pobCombinada[i].electrons = (Electron *)malloc(sizeof(Electron) * numElectrones);
    pobCombinada[i + tamPob].electrons = (Electron *)malloc(sizeof(Electron) * numElectrones);

    for(int j=0; j<numElectrones; ++j)
	{
      cartesianCoordinates(&descendencia[i].electrons[j]);

      for(int k=0; k<tamPobCombinada; k+=tamPob)
	  {
        pobCombinada[i+k].electrons[j].surface = (double *)malloc(sizeof(double)*dimV);
        pobCombinada[i+k].electrons[j].cartesian = (double *)malloc(sizeof(double)*dimE);
      }

      for(int k=0; k<dimV; ++k)
	  {
        pobCombinada[i].electrons[j].surface[k] = poblacion[i].electrons[j].surface[k];
        pobCombinada[i+tamPob].electrons[j].surface[k] = descendencia[i].electrons[j].surface[k];
      }

      for(int k=0; k<dimE; ++k)
	  {
        pobCombinada[i].electrons[j].cartesian[k] = poblacion[i].electrons[j].cartesian[k];
        pobCombinada[i+tamPob].electrons[j].cartesian[k] = descendencia[i].electrons[j].cartesian[k];
      }
    }

    configurationEnergy(&descendencia[i], dimE, numElectrones);
    pobCombinada[i].energy = poblacion[i].energy;
    pobCombinada[i+tamPob].energy = descendencia[i].energy;
  }

  // Elitismo
  int indBest = 0;
  for(int i=tamPob; i<tamPobCombinada; ++i) 
  {
    if(pobCombinada[i].energy < pobCombinada[indBest].energy) // minimizacion
    {
      indBest = i;
    }
  }

  for(int i=0; i<numElectrones; ++i)
  {
    for(int j=0; j<dimV; ++j)
	{
      poblacion[0].electrons[i].surface[j] = pobCombinada[indBest].electrons[i].surface[j];
      pobCombinada[indBest].electrons[i].surface[j] = pobCombinada[0].electrons[i].surface[j];
      pobCombinada[0].electrons[i].surface[j] = poblacion[0].electrons[i].surface[j];
    }

    for (int j=0; j<dimE; ++j)
	{
      poblacion[0].electrons[i].cartesian[j] = pobCombinada[indBest].electrons[i].cartesian[j];
      pobCombinada[indBest].electrons[i].cartesian[j] = pobCombinada[0].electrons[i].cartesian[j];
      pobCombinada[0].electrons[i].cartesian[j] = poblacion[0].electrons[i].cartesian[j];
    }
  }

  poblacion[0].energy = pobCombinada[indBest].energy;
  pobCombinada[indBest].energy = pobCombinada[0].energy;
  pobCombinada[0].energy = poblacion[0].energy;

  // Ruleta
  for(int i=0; i<tamPobCombinada; ++i)
  {
    probabilidades[i] = 1.0/pobCombinada[i].energy;
    normalizacion += probabilidades[i];
  }

  probabilidades[0] /= normalizacion;
  for(int i=1; i<tamPobCombinada; ++i)
  {
    probabilidades[i] /= normalizacion;
    probabilidades[i] += probabilidades[i-1];
  }

  for (int i=1; i<tamPob; ++i)
  {
    int p = 0;
    double rn = (double)rand()/RAND_MAX;

    while (probabilidades[p]<rn)
	{
      ++p;
    }
    for(int j=0; j<numElectrones; ++j)
	{
      for (int k=0; k<dimV; ++k)
	  {
        poblacion[i].electrons[j].surface[k] = pobCombinada[p].electrons[j].surface[k];
      }

      for (int k=0; k<dimE; ++k)
	  {
        poblacion[i].electrons[j].cartesian[k] = pobCombinada[p].electrons[j].cartesian[k];
      }
    }
    poblacion[i].energy = pobCombinada[p].energy;
  }

  for(int i=0; i<tamPobCombinada; ++i)
  {
    for(int j=0; j<numElectrones; ++j)
    {
      free(pobCombinada[i].electrons[j].surface);
      free(pobCombinada[i].electrons[j].cartesian);
    }
    free(pobCombinada[i].electrons);
  }
  free(pobCombinada);
  free(probabilidades);
}

// cuerpo principal algoritmo genetico
void geneticAlgorithm(Configuration * bestSol, int * tBest, int dimE, int dimV, double * UMIN, double * UMAX, int numGen, int tamPob, int numElectrones, double probaCruza, double probaMutacion)
{
  // Parametros auxiliares
  int tamPob2 = 2*tamPob;
  double * URANGO = (double *)malloc(sizeof(double)*dimV);
  for (int i=0; i<dimV; ++i)
  {
    URANGO[i] = UMAX[i] - UMIN[i];
  }

  // salidas
  double energyBest = INFINITY; // problema de minimizacion
  Configuration * poblacion = (Configuration *)malloc(sizeof(Configuration)*tamPob);
  Configuration * descendencia = (Configuration *)malloc(sizeof(Configuration)*tamPob);

  for(int i=0; i<tamPob; ++i)
  {
    poblacion[i].electrons = (Electron *)malloc(sizeof(Electron)*numElectrones);
    descendencia[i].electrons = (Electron *)malloc(sizeof(Electron)*numElectrones);

    for(int j=0; j<numElectrones; ++j)
	{
      poblacion[i].electrons[j].surface = (double *)malloc(sizeof(double)*dimV);
      descendencia[i].electrons[j].surface = (double *)malloc(sizeof(double)*dimV);
      poblacion[i].electrons[j].cartesian = (double *)malloc(sizeof(double)*dimE);
      descendencia[i].electrons[j].cartesian = (double *)malloc(sizeof(double)*dimE);
    }
  }

  poblacionInicial(poblacion, dimE, dimV, tamPob, numElectrones, UMIN, URANGO);

  for(int i=0; i<tamPob; ++i)
  {
    if(poblacion[i].energy < energyBest)
	{
      energyBest = poblacion[i].energy;
    }
  }

  printf("\tGeneracion\t\tMenor Energia\n");
  printf("\t     %d     \t\t%.9f\n", 0, energyBest);

  for(int t=1; t<numGen+1; ++t)
  {
    // Operadores geneticos
    recombination(poblacion, descendencia, dimV, tamPob, probaCruza, numElectrones);
    mutation(descendencia, dimV, tamPob, probaMutacion, t, numGen, numElectrones, UMIN, UMAX, URANGO);
    selection(poblacion, descendencia, dimE, dimV, tamPob, tamPob2, numElectrones);

    if(poblacion[0].energy<energyBest)
	{
      energyBest = poblacion[0].energy;
      *tBest = t;
    }

    if(t%(numGen/10) == 0)
	{
      printf("\t     %d     \t\t%.9f\n", t, energyBest);
    }
  }

  bestSol->energy = poblacion[0].energy;
  for(int i=0; i<numElectrones; ++i)
  {
    for(int j=0; j<dimV; ++j)
	{
      bestSol->electrons[i].surface[j] = poblacion[0].electrons[i].surface[j];
    }

    for(int j=0; j<dimE; ++j)
	{
      bestSol->electrons[i].cartesian[j] = poblacion[0].electrons[i].cartesian[j];
    }
  }

  // Liberar memoria
  for(int i=0; i<tamPob; ++i)
  {
    for(int j=0; j<numElectrones; ++j)
    {
      free(poblacion[i].electrons[j].surface);
      free(descendencia[i].electrons[j].surface);
      free(poblacion[i].electrons[j].cartesian);
      free(descendencia[i].electrons[j].cartesian);
    }
    free(poblacion[i].electrons);
    free(descendencia[i].electrons);
  }
  free(URANGO);
  free(poblacion);
  free(descendencia);
}

// llamada a python para graficar
void generateGraph()
{
  char * pythonScript;

  // Esfera
  /*"U_MIN = 0.0\n"
  "U_MAX = np.pi\n"
  "V_MIN = 0.0\n"
  "V_MAX = 2.0*np.pi\n\n"

  "R = 1.0\n"
  "u, v = np.mgrid[U_MIN:U_MAX:500j, V_MIN:V_MAX:500j]\n\n"

  "x = R*np.sin(v)*np.cos(u)\n"
  "y = R*np.sin(v)*np.sin(u)\n"
  "z = R*np.cos(v)\n\n"*/

  // Cilindro
  /*"U_MIN = 0\n"
  "U_MAX = 2*np.pi\n"
  "V_MIN = 0\n"
  "V_MAX = 4\n\n"

  "R = 1.0\n"
  "H = V_MAX - V_MIN\n"
  "u, v = np.mgrid[U_MIN:U_MAX:100j, V_MIN:V_MAX:100j]\n\n"

    "x = R*np.cos(u)\n"
    "y = R*np.sin(u)\n"
    "z = v\n\n"*/

  // Toro
  /*"U_MIN = 0\n"
  "U_MAX = 2*np.pi\n"
  "V_MIN = 0\n"
  "V_MAX = 2*np.pi\n\n"

  "A = 4.0\n"
  "B = 1.0\n"
  "u, v = np.mgrid[U_MIN:U_MAX:100j, V_MIN:V_MAX:100j]\n\n"

    "x = (A+B*np.cos(u))*np.cos(v)\n"
    "y = (A+B*np.cos(u))*np.sin(v)\n"
    "z = B*np.sin(u)\n\n"*/

	pythonScript =  "import numpy as np\n"
					"import matplotlib.pyplot as plt\n\n"
					
					"U_MIN = 0\n"
					"U_MAX = 2*np.pi\n"
					"V_MIN = 0\n"
					"V_MAX = 2*np.pi\n\n"
					
					"A = 4.0\n"
					"B = 1.0\n"
					"u, v = np.mgrid[U_MIN:U_MAX:250j, V_MIN:V_MAX:250j]\n\n"
					
					"x = (A+B*np.cos(u))*np.cos(v)\n"
					"y = (A+B*np.cos(u))*np.sin(v)\n"
					"z = B*np.sin(u)\n\n"
					
					"data = np.loadtxt('coordinates.txt')\n"
					"ex = data[:, 0]\n"
					"ey = data[:, 1]\n"
					"ez = data[:, 2]\n\n"
					
					"field = np.arctan(np.exp(sum(1.0/np.sqrt((x-ex[i])**2 + (y-ey[i])**2 + (z-ez[i])**2) for i in range(len(ex)))))\n"
					"field = (field - field.min()) / (field.max() - field.min())\n"
					"col = plt.cm.jet(field)\n\n"
					
					"fig = plt.figure()\n"
					"ax = fig.add_subplot(111, projection='3d')\n"
					"ax.plot_surface(x, y, z, facecolors = col, rstride = 1, cstride = 1, alpha = 1, antialiased = True, shade = False)\n\n"
					
					"ax.set_box_aspect([4, 4, 1])\n"
					"ax.grid(False)\n"
					"ax.axis('off')\n\n"
					
					"plt.rcParams['figure.figsize'] = [7.00, 7.00]\n"
					"plt.rcParams['figure.autolayout'] = True\n\n"
					
					"plt.show()";
					
					FILE * file = fopen("graficaThomson.py", "w");
					fprintf(file, "%s", pythonScript);
					fclose(file);
					
					system("python graficaThomson.py");
					free(pythonScript);
}

int main(int argc, char * argv[])
{
  srand(time(NULL));

  // Numero de electrones, no es parametro
  int numElectrones = 32;

  // Parametros superficie
  int dimEspacio = 3;
  int dimVariedad = 2;
  double * UMIN = (double *)malloc(sizeof(double)*dimVariedad);
  double * UMAX = (double *)malloc(sizeof(double)*dimVariedad);
  // Esfera
  // UMAX[0] = M_PI; UMAX[1] = 2.0*M_PI;
  // UMIN[0] = 0.0; UMIN[1] = 0.0;
  // Cilindro
  // UMAX[0] = 2.0*M_PI; UMAX[1] = 4.0;
  // UMIN[0] = 0.0; UMIN[1] = 0.0;
  // Toro
  UMAX[0] = 2.0*M_PI;
  UMAX[1] = 2.0*M_PI;
  UMIN[0] = 0.0;
  UMIN[1] = 0.0;

  // Parametros algoritmo
  int numGeneraciones = 10000;
  int tamanioPoblacion;
  double probaCruza;
  double probaMutacion;

  // Salidas
  int mejorGen;
  double tiempo;
  Configuration bestSol;
  bestSol.electrons = (Electron *)malloc(sizeof(Electron)*numElectrones);
  for(int i=0; i<numElectrones; ++i)
  {
    bestSol.electrons[i].surface = (double *)malloc(sizeof(double)*dimVariedad);
    bestSol.electrons[i].cartesian = (double *)malloc(sizeof(double)*dimEspacio);
  }

  if(argc>1) // Calibracion
  {
    tamanioPoblacion = 2*(int)atof(argv[1]);
    probaCruza = atof(argv[2]);
    probaMutacion = atof(argv[3]);

    geneticAlgorithm(&bestSol, &mejorGen, dimEspacio, dimVariedad, UMIN, UMAX, numGeneraciones, tamanioPoblacion, numElectrones, probaCruza, probaMutacion);

    FILE * energyFile = fopen("energyval.txt", "w");
    fprintf(energyFile, "%lf", bestSol.energy);
    fclose(energyFile);
  } 
  else // Corrida normal
  {
    tamanioPoblacion = 8;
    probaCruza = 0.5;
    probaMutacion = 0.5;

    clock_t start = clock();
    geneticAlgorithm(&bestSol, &mejorGen, dimEspacio, dimVariedad, UMIN, UMAX, numGeneraciones, tamanioPoblacion, numElectrones, probaCruza, probaMutacion);
    tiempo = (double)(clock() - start) / CLOCKS_PER_SEC;

    printf("\n\nMejor solucion encontrada:\n");
    printf("\nEnergia:   \t\t %.12lf", bestSol.energy);
    printf("\nGeneracion:\t\t %d", mejorGen);
    printf("\nTiempo de Ejecucion:\t %f segundos\n\n", tiempo);

    if(dimEspacio==3 && dimVariedad==2)
	{
      FILE * file = fopen("coordinates.txt", "w");

      for(int i=0; i<numElectrones; ++i)
	  {
        fprintf(file, "%lf %lf %lf\n", bestSol.electrons[i].cartesian[0], bestSol.electrons[i].cartesian[1], bestSol.electrons[i].cartesian[2]);
      }
      fclose(file);
      generateGraph();
    }
  }

  for(int i=0; i<numElectrones; ++i)
  {
    free(bestSol.electrons[i].surface);
    free(bestSol.electrons[i].cartesian);
  }
  free(UMIN);
  free(UMAX);
  free(bestSol.electrons);

  return 69;
}

