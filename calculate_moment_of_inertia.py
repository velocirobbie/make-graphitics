import sys
import math

xyzfile = sys.argv[1]

axis = [0,0,0]

with open(xyzfile) as f:
    N = int(f.readline())
    f.readline()
    MoI = 0
    mass = 0
    for i in range(N):
        element,x,y,z = f.readline().split()
        rx = float(x) - axis[0]
        ry = float(y) - axis[1]
        rz = float(z) - axis[2]
        r2 = rx**2 + ry**2 + rz**2
        if element == 'C': m = 12
        elif element == 'H': m = 1
        else: print 'Unkown element: ', atom
        moment = r2 * m
        MoI += moment

        mass += m
print 'Moment of inertia = ',MoI
print 'Total mass = ',mass




