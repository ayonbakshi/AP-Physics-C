from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import numpy as np
import matplotlib.pyplot as plt

# USER DEFINABLE PARAMETERS BEGIN HERE
wire_resolution = 100 # how fine to render the wire
field_resolution = 7 # how large the field is
ds = 0.1 # integration step size

# write x, y, z as parametric functions of t
def x(t):
    return 3 * np.sin(t * 200)
def y(t):
    return 3 * np.cos(t * 200)
def z(t):
    return 0.5 - t
# current in the wire
I = 1
# induction constant ( = mu_0/4pi)
C = 5 #10**-7

# rendering boundaries
x_lo, x_hi = -3, 3
y_lo, y_hi = -3, 3
z_lo, z_hi = -1, 1

# END USER DEFINABLE PARAMETERS
# time complexity: O(WF^3) where W, F are the wire and field resolutions respectively

# solenoid
t = np.linspace(0, 1, wire_resolution)
tt = np.linspace(0, 1, wire_resolution*100)
xx = x(tt)
yy = y(tt)
zz = z(tt)

# B field = (BU, BV, BW)
BU = np.zeros((field_resolution, field_resolution, field_resolution))
BV = np.zeros((field_resolution, field_resolution, field_resolution))
BW = np.zeros((field_resolution, field_resolution, field_resolution))

X, Y, Z = np.meshgrid(np.linspace(x_lo, x_hi, field_resolution),
                      np.linspace(y_lo, y_hi, field_resolution),
                      np.linspace(z_lo, z_hi, field_resolution))

print("computing... 0%")
prg, cprg = 0, 0.1
for i in range(field_resolution):
    for j in range(field_resolution):
        for k in range(field_resolution):
            prg += 1
            if prg/(field_resolution**3)>cprg:
                print("computing... %d%%"%(cprg*100))
                cprg+=0.1
            point = [X[i][j][k], Y[i][j][k], Z[i][j][k]]
            for prv, nxt in np.dstack((t[:-1], t[1:]))[0]:
                mag = np.linalg.norm(np.array([x(nxt) - x(prv),
                                               y(nxt) - y(prv),
                                               z(nxt) - z(prv)]))
                npi = int(np.ceil(mag / ds))
                if npi <= 2:
                    print("error: integration step is very large!")
                lx = np.linspace(x(prv), x(nxt), npi)
                ly = np.linspace(y(prv), y(nxt), npi)
                lz = np.linspace(z(prv), z(nxt), npi)
                rx = lx[:-1] - point[0]
                ry = ly[:-1] - point[1]
                rz = lz[:-1] - point[2]
                dlx = np.diff(lx)
                dly = np.diff(ly)
                dlz = np.diff(lz)
                r = np.linalg.norm(np.array([rx, ry, rz]))
                dU = C * (dly * rz - dlz * ry)/r/r/r
                dV = C * (dlz * rx - dlx * rz)/r/r/r
                dW = C * (dlx * ry - dly * rx)/r/r/r
                BU[i][j][k] += np.sum(dU)
                BV[i][j][k] += np.sum(dV)
                BV[i][j][k] += np.sum(dW)

plt.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

# plot solenoid
ax.plot(xx, yy, zz, label='solenoid 1')
ax.legend()

# plot field
ax.quiver(X, Y, Z, BU, BV, BW)#, normalize=True)
#SX, SY, SZ = np.meshgrid(np.linspace(-5, 5, 20),
#                         np.linspace(-5, 5, 20),
#                         np.linspace(-1, -1, 1))
#ax.streamplot(X, Y, Z, BU, BV, BW, SX, SY, SZ)
# TODO: figure out how to plot streams in 3d in matplotlib
# 3d stream lines were lost in translation between MATLAB and python :(
plt.show()
