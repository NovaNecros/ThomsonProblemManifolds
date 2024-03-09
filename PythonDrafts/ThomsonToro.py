import numpy as np
import numpy.random as rn
import numpy.linalg as la
import matplotlib.pyplot as plt
import time as tm

%matplotlib qt

A = 4
B = 1

#clase para representar un electrón
class Electron:
    def __init__(self, u, v):
        self.u = u
        self.v = v
        self.normalize()

	  # normalizamos los valores para asegurarse de que están dentro de los rangos adecuados
    def normalize(self):
        while self.u > 2*np.pi:
            self.u -= 2*np.pi
        while self.u < 0:
            self.u += 2*np.pi
        while self.v > 2*np.pi:
            self.v -= 2*np.pi
        while self.v < 0:
            self.v += 2*np.pi
 

		#convertimos las coordenadas cilíndricas a coordenadas cartesianas
    def cartesian(self):
        x = (A+B*np.cos(self.u))*np.cos(self.v)
        y = (A+B*np.cos(self.u))*np.sin(self.v)
        z = B*np.sin(self.u)
        return x,y,z

    def __str__(self):
        x,y,z = self.cartesian()
        return "[{:.6f}, {:.6f}, {:.6f}]".format(x, y, z)

# Función para calcular la distancia entre dos electrones
def distancia(e1, e2):
    return la.norm(np.array(e1.cartesian())-np.array(e2.cartesian()))

# graficamos los electrones en un espacio tridimensional
def grafica(electrons):
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, 2*np.pi, 100)
    x = np.outer((A+B*np.cos(u)),np.cos(v))
    y = np.outer((A+B*np.cos(u)),np.sin(v))
    z = np.outer(B*np.sin(u),np.ones(np.size(v)))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(x, y, z, color="red", alpha=0.2, edgecolor = (0,0,0,0.1), linewidth = 0.4)

    for e in electrons:
        ax.scatter(e.cartesian()[0], e.cartesian()[1], e.cartesian()[2], color="lime", s=100)

    ax.set_box_aspect([4, 4, 1])

    ax.grid(False)
    ax.axis("off")

    plt.rcParams["figure.figsize"] = [7.00, 7.00]
    plt.rcParams["figure.autolayout"] = True
    
    plt.show()

    return

# generamos una población inicial
def poblacionInicial(tam_pob, num_e):
    electrones = [[] for _ in range(tam_pob)]

    for i in range(tam_pob):
        for j in range(num_e):
            electrones[i].append(Electron(rn.uniform(0, 2*np.pi), rn.uniform(0, 2*np.pi)))

    return electrones


# Función para calcular la energía de un conjunto de electrones
def energia(electrones):
    s = 0
    for i in range(len(electrones)):
        for j in range(i):
            d = distancia(electrones[i], electrones[j])
            if d == 0:
                return np.inf
            s += 1 / d
    return s

def recombinacion(A, PC):
    m, n = len(A), len(A[0])
    H = [[] for _ in range(m)]
    ind = rn.permutation(m)

    for i in range(0, m, 2):
        for j in range(n):
            if rn.uniform(0, 1) < PC:
                H[i].append(Electron(A[ind[i+1]][j].u, A[ind[i]][j].v))
                H[i+1].append(Electron(A[ind[i]][j].u, A[ind[i+1]][j].v))
            else:
                H[i].append(Electron(A[ind[i]][j].u, A[ind[i]][j].v))
                H[i+1].append(Electron(A[ind[i+1]][j].u, A[ind[i+1]][j].v))
    return H

def mutacion(A, PM, u_sd, v_sd):
    m, n = len(A), len(A[0])
    M = [[] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            if rn.uniform(0, 1) < PM:
                M[i].append(Electron(A[i][j].u + rn.normal(0, u_sd),
                                     A[i][j].v + rn.normal(0, v_sd)))
            else:
                M[i].append(Electron(A[i][j].u, A[i][j].v))

    return M

def seleccion(A, fA):
    m = len(A) // 2
    fS = [min(fA)]
    S = [A[fA.index(fS[0])]]

    gA = [1 / fa for fa in fA]
    suma = sum(gA)
    proba = [ga / suma for ga in gA]

    for i in range(1, 2 * m):
        proba[i] += proba[i - 1]

    for i in range(1, m):
        ind = np.searchsorted(proba, rn.uniform(0, 1))
        fS.append(fA[ind])
        S.append(A[ind])

    return S, fS

def algoritmoGenetico(num_gen, tam_pob, num_electrones, proba_cruza, proba_mutacion):
    A = poblacionInicial(tam_pob, num_electrones)
    fA = [energia(x) for x in A]
    t_best = 0
    fs = min(fA)

    print("\n\tProblema de Thomson")
    print("\tNumero de generaciones =", num_gen)
    print("\n\tAlgoritmo Genetico")
    print("\tTamano de la poblacion =", tam_pob)
    print("\tNumero de electrones =", num_electrones)
    print("\tProbabilidad de cruza =", proba_cruza)
    print("\tProbabilidad de mutacion =", proba_mutacion)

    print("\n Generacion", 0)
    print("\t{:^18}".format("Energia"))
    for i in range(tam_pob):
        print("\t{:^18.12f}".format(fA[i]))

    for t in range(1,num_gen+1):
        B = recombinacion(A, proba_cruza)
        B = mutacion(B, proba_mutacion, 2*np.pi/t, 2*np.pi/t)
        fB = [energia(e) for e in B]
        A, fA = seleccion(A + B, fA + fB)

        if t % (num_gen // 4) == 0:
            print("\n\n Generacion", t)
            print("\t{:^18}".format("Energia"))
            for i in range(tam_pob):
                print("\t{:^18.12f}".format(fA[i]))

        if fA[0] < fs:
            fs = fA[0]
            t_best = t

    return t_best, A[0], fA[0]

rn.seed()
PC, PM = 0.5, 0.5
N, m, n = 1000, 10, 12

start = tm.time()
ts, s, fs = algoritmoGenetico(N, m, n, PC, PM)
stop = tm.time()



print("\n\tConfiguracion:")
for e in s:
    print("\t\t", e)


print("\n   Energia Minima Encontrada:\t", round(fs, 12))
print("\n   Encontrada en la Generacion:\t", ts)
print("\n\n\nTiempo de ejecucion:\t", round(stop - start, 4), "segundos")

grafica(s)