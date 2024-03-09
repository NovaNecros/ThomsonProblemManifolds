import numpy as np
import numpy.random as rn
import numpy.linalg as la
import matplotlib.pyplot as plt
import time as tm

%matplotlib qt
R = 1

# genes, puntos en coordenadas de la superficie y cartesianas
class Electron:
    def __init__(self, polar, elevation):
        self.polar = polar
        self.elevation = elevation
        self.normalization()

    def normalization(self):
        while self.polar > 2*np.pi:
            self.polar -= 2*np.pi
        while self.polar < 0:
            self.polar += 2*np.pi
        while self.elevation > 2*np.pi:
            self.elevation -= 2*np.pi
        while self.elevation < 0:
            self.elevation += 2*np.pi
        if self.elevation > np.pi:
            self.elevation = 2*np.pi - self.elevation
            self.polar = -self.polar
        

    def cartesian(self):
        return [R * np.cos(self.elevation) * np.sin(self.polar), R * np.sin(self.elevation) * np.sin(self.polar), R * np.cos(self.polar)]
			
    def __str__(self):
        x,y,z = self.cartesian()
        return "[{:.6f}, {:.6f}, {:.6f}]".format(x, y, z)

# función objetivo
def energy(electrons):
    s = 0
    for i in range(len(electrons)):
        for j in range(i):
            d = distance(electrons[i], electrons[j])
            if d == 0:
                return np.inf
            s += 1 / d
    return s

# distancia euclidiana
def distance(e1, e2):
    return la.norm(np.array(e1.cartesian())-np.array(e2.cartesian()))

# graficar resultados
def graficar(electrons):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x, y, z, color="red", alpha=0.2, edgecolor = (0,0,0,0.1), linewidth = 0.4)

    for e in electrons:
        ax.scatter(e.cartesian()[0], e.cartesian()[1], e.cartesian()[2], color="lime", s=100)

    ax.set_box_aspect([1, 1, 1])

    ax.grid(False)
    ax.axis("off")

    plt.rcParams["figure.figsize"] = [7.00,7.00]
    plt.rcParams["figure.autolayout"] = True
    
    plt.show()

    return

# población inicial
def population_initialization(pop_size, num_electrons):
    electrons = [[Electron(0, 0)] for _ in range(pop_size)]
    
    for i in range(pop_size):
        for j in range(1, num_electrons):
            electrons[i].append(Electron(rn.uniform(0, 2*np.pi), rn.uniform(0, np.pi)))
    
    return electrons


# cruza gen por gen con media aritmetica y armonica
def recombination(A, PC):
    m, n = len(A), len(A[0])
    ind = rn.permutation(m)
    H = [[Electron(0, 0)] for i in range(m)]
    
    for i in range(0, m, 2):
        for j in range(1, n):
            if rn.uniform(0, 1) < PC:
                H[i].append(Electron((A[ind[i]][j].polar + A[ind[i+1]][j].polar) / 2, (A[ind[i]][j].elevation + A[ind[i+1]][j].elevation) / 2))
                H[i+1].append(Electron(2 / (1/A[ind[i]][j].polar + 1/A[ind[i+1]][j].polar), 2 / (1/A[ind[i]][j].elevation + 1/A[ind[i+1]][j].elevation)))	
            else:
                H[i].append(Electron(A[ind[i]][j].polar, A[ind[i]][j].elevation))
                H[i+1].append(Electron(A[ind[i+1]][j].polar, A[ind[i+1]][j].elevation))
    
    return H

# mutaciones aleatorias con perturbacion gaussiana
def mutation(A, PM, sd_polar, sd_elevation):
    m, n = len(A), len(A[0])
    M = [[Electron(0, 0)] for _ in range(m)]
    	
    for i in range(m):
        for j in range(1, n):
            if rn.uniform(0, 1) < PM:
                M[i].append(Electron(A[i][j].polar + rn.normal(0, sd_polar), A[i][j].elevation + rn.normal(0, sd_elevation)))
            else:
                M[i].append(Electron(A[i][j].polar, A[i][j].elevation))
    
    return M



# seleccion aleatoria con elitismo
def selection(A, fA):
    m = len(A) // 2
    fS = [min(fA)]
    S = [A[fA.index(fS[0])]]
    
    gA = [1 / fa for fa in fA]
    suma = sum(gA)
    proba = [ga / suma for ga in gA]
    
    for i in range(1, 2 * m):
        proba[i] += proba[i-1]
    
    for i in range(1, m):
        ind = np.searchsorted(proba, rn.uniform(0, 1))
        fS.append(fA[ind])
        S.append(A[ind])
    
    return S, fS

# cuerpo principal del algoritmo genetico 
def genetic_algorithm(num_gen, pop_size, num_electrons, prob_recombination, prob_mutation):
    A = population_initialization(pop_size, num_electrons)
    fA = [energy(x) for x in A]
    t_best = 0
    fs = min(fA)
    
    print("\n\tProblema de Thomson")
    print("\tNumero de generaciones =", num_gen)
    print("\n\tAlgoritmo Genetico")
    print("\tTamano de la poblacion =", pop_size)
    print("\tNumero de electrones =", num_electrons)
    print("\tProbabilidad de cruza =", prob_recombination)
    print("\tProbabilidad de mutacion =", prob_mutation)
    
    print("\n Generacion", 0)
    print("\t{:^18}".format("Energia"))
    for i in range(pop_size):
        print("\t{:^18.12f}".format(fA[i]))
    
    for t in range(1,num_gen+1):
        sd_polar = 2*np.pi / t
        sd_elevation = np.pi / t
        B = recombination(A, prob_recombination)
        B = mutation(B, prob_mutation, sd_polar, sd_elevation)
        fB = [energy(e) for e in B]
        A, fA = selection(A + B, fA + fB)
        
        if t % (num_gen // 4) == 0:
            print("\n\n Generacion", t)
            print("\t{:^18}".format("Energia"))
            for i in range(pop_size):
                print("\t{:^18.12f}".format(fA[i]))
        
        if fA[0] < fs:
            fs = fA[0]
            t_best = t
    
    return t_best, A[0], fA[0]

rn.seed()
PC, PM = 0.5, 0.5
N, m, n = 10000, 10, 12

start = tm.time()
ts, s, fs = genetic_algorithm(N, m, n, PC, PM)
stop = tm.time()

print("\n\n   Energia Minima Encontrada:")
print("\n\tConfiguracion:")
for e in s:
    print("\t\t", e)

print("\n   Energia Minima:\t", round(fs, 12))
print("\n   Encontrada en la Generacion:\t", ts)
print("\n   Radio:\t ", R)
print("\n\n\nTiempo de ejecucion:\t", round(stop - start, 4), "segundos")

graficar(s)
