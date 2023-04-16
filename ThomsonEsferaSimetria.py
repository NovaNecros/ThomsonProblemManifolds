#ESFERA UNITARIA CON SIMETRIA

from numpy.random import uniform, normal, shuffle, seed
from numpy import pi, sqrt, sin, cos, inf, mgrid
from matplotlib import pyplot as plt
from time import time

_2pi = 2*pi
polar_sd = pi / 12
azimuth_sd = pi / 24

class Electron:
	def __init__(self, polar, azimuth):
		self.polar = polar
		self.azimuth = azimuth
		
		while self.polar > _2pi:
			self.polar -= _2pi
		while self.polar < 0:
			self.polar += _2pi
		while self.azimuth > _2pi:
			self.azimuth -= _2pi
		while self.azimuth < 0:
			self.azimuth += _2pi
		if self.azimuth > pi:
			self.azimuth = _2pi - self.azimuth

	def cartessian(self):
		return [cos(self.azimuth)*sin(self.polar), sin(self.azimuth)*sin(self.polar), cos(self.polar)]

	def ElectronToString(self):
		return "(" + str(round(self.polar,6)) + ", " + str(round(self.azimuth,6)) + ")"


def energia(electrones):
	s = 0
	for i in range(len(electrones)):
		for j in range(i):
			d = distancia(electrones[i], electrones[j])
			if d == 0:
				return inf
			s += 1/d
	return s

def ListToString(L):
	n = len(L)
	s = "("
	for j in range(n):
		s += str(round(L[j],6))
		if j < n-1:
			s += ", "
		else:
			s += ")"
	return s

def busquedaBinaria(L, x):
	return BBrec(L,x,0,len(L))

def BBrec(L, x, inicio, fin):
	m = (inicio + fin) // 2
	if x == L[m] or fin-inicio < 2:
		return m
	elif x < L[m]:
		return BBrec(L, x, inicio, m)
	elif x > L[m]:
		return BBrec(L, x, m, fin)
		
def distancia(e1, e2):
	return sqrt(2-2*(sin(e1.polar)*sin(e2.polar)*cos(e1.azimuth-e2.azimuth)+cos(e1.polar)*cos(e2.polar)))

def grafica(electrones):
	plt.rcParams["figure.figsize"] = [7.00, 3.50]
	plt.rcParams["figure.autolayout"] = True
	fig = plt.figure()
	ax = fig.add_subplot(projection = "3d")
	u, v = mgrid[0:2*pi:30j, 0:pi:20j]
	x = cos(u)*sin(v)
	y = sin(u)*sin(v)
	z = cos(v)
	ax.plot_surface(x, y, z, alpha = 0.2, color = "red")

	for e in electrones:
		ax.scatter(e[0], e[1], e[2], color = "lime", s = 100)
	
	ax.grid(False)
	ax.axis("off")
	plt.show()

def poblacionInicial(tam_pob, num_e):
	electrones = [[Electron(0,0)] for i in range(tam_pob)]
	
	for i in range(tam_pob):
		for j in range(1,num_e):
			electrones[i].append(Electron(uniform(0,_2pi), uniform(0,pi)))
			
	return electrones

def recombinacion(A, PC):
	m,n = len(A), len(A[0])
	ind = [i for i in range(m)]
	shuffle(ind)
	H = [[Electron(0,0)] for i in range(m)]

	for i in range(0,m,2):
		for j in range(1,n):
			if uniform(0,1) < PC:
				H[i].append(Electron((A[ind[i]][j].polar+A[ind[i+1]][j].polar)/2, (A[ind[i]][j].azimuth+A[ind[i+1]][j].azimuth)/2))
				H[i+1].append(Electron(A[ind[i+(j%2)]][j].polar, A[ind[i+(j%2)]][j].azimuth))
			else: 
				H[i].append(Electron(A[ind[i]][j].polar, A[ind[i]][j].azimuth))
				H[i+1].append(Electron(A[ind[i+1]][j].polar, A[ind[i+1]][j].azimuth))
						
	return H

def mutacion(A, PM):
	m,n = len(A),len(A[0])
	M = [[Electron(0,0)] for i in range(m)]
	for i in range(m):
		for j in range(1,n):
			if uniform(0,1) < PM:
				M[i].append(Electron(A[i][j].polar+normal(0,polar_sd), A[i][j].azimuth+normal(0,azimuth_sd)))
			else:
				M[i].append(Electron(A[i][j].polar, A[i][j].azimuth))
		
	return M

def seleccion(A, fA):
	m = len(A)//2
	fS = [min(fA)]
	S = [A[fA.index(fS[0])]]
	
	gA = [1/fa for fa in fA]
	suma = sum(gA)
	proba = [ga / suma for ga in gA]

	for i in range(1, 2*m):
		proba[i] += proba[i-1]
	
	for i in range(1,m):
		ind = busquedaBinaria(proba, uniform(0,1))
		fS.append(fA[ind])
		S.append(A[ind])
	
	return S,fS

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
	
	for t in range(num_gen):
		B = recombinacion(A, proba_cruza)
		B = mutacion(B, proba_mutacion)
		fB = [energia(e) for e in B]
		A,fA = seleccion(A+B,fA+fB)

		#if (t+1) % (num_gen//4) == 0:
		print("\n\n Generacion", t+1)
		print("\t{:^18}".format("Energia"))
		for i in range(tam_pob):
			print("\t{:^18.12f}".format(fA[i]))

		if fA[0] < fs:
			fs = fA[0]
			t_best = t+1
	
	return t_best, A[0], fA[0]


seed()

PC, PM = 0.5, 0.5
N, m, n = 1000, 10, 6
start = time()
ts, s, fs = algoritmoGenetico(N,m,n,PC,PM)
stop = time()
p = []

print("\n\n   Mejor Solucion:")
print("\n\tCoordenadas Esfericas:")
for e in s:
	print("\t\t", e.ElectronToString())
	
print("\n\tCoordenadas Cartesianas:")
for e in s:
	p.append(e.cartessian())
	print("\t\t", ListToString(p[-1]))
	
print("\n   Mejor valor F.O.:\t", round(fs,12))
print("\n   Encontrado en la Generacion:\t", ts)
print("\n\n\nTiempo de ejecucion:\t", round(stop-start,4), "segundos")

grafica(p)
