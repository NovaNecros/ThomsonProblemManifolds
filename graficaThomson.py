import numpy as np
import matplotlib.pyplot as plt

U_MIN = 0
U_MAX = 2*np.pi
V_MIN = 0
V_MAX = 2*np.pi

A = 4.0
B = 1.0
u, v = np.mgrid[U_MIN:U_MAX:250j, V_MIN:V_MAX:250j]

x = (A+B*np.cos(u))*np.cos(v)
y = (A+B*np.cos(u))*np.sin(v)
z = B*np.sin(u)

data = np.loadtxt('coordinates.txt')
ex = data[:, 0]
ey = data[:, 1]
ez = data[:, 2]

field = np.arctan(np.exp(sum(1.0/np.sqrt((x-ex[i])**2 + (y-ey[i])**2 + (z-ez[i])**2) for i in range(len(ex)))))
field = (field - field.min()) / (field.max() - field.min())
col = plt.cm.jet(field)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, facecolors = col, rstride = 1, cstride = 1, alpha = 1, antialiased = True, shade = False)

ax.set_box_aspect([4, 4, 1])
ax.grid(False)
ax.axis('off')

plt.rcParams['figure.figsize'] = [7.00, 7.00]
plt.rcParams['figure.autolayout'] = True

plt.show()