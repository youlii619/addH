import os
import sys

with open(sys.argv[1], 'r') as f:
     B = f.readlines()

hydrC = sys.argv[2]
num_atoms = int(len(B)-2)
atom_list = []

for i in range(num_atoms):
    k = B[i+2].split()
    temp = []
    temp.append(k)
    temp2 = []
    t0 = temp[0][0].strip('\n')
    temp2.append(t0)
    t1 = temp[0][1].strip('\n')
    temp2.append(t1)
    t2 = temp[0][2].strip('\n')
    temp2.append(t2)
    t3 = temp[0][3].strip('\n')
    temp2.append(t3)
    atom_list.append(temp2)
#print(atom_list)
f.close

def dist(i,j):
    return ((float(i[1])-float(j[1]))**2 + (float(i[2])-float(j[2]))**2 + (float(i[3])-float(j[3]))**2)**0.5

##### Making neighbour map (ngmap) including H #####

ngmap = {}

for i in atom_list:
    ngmap[int(atom_list.index(i))] = []
    for j in atom_list:
        if atom_list.index(j) != atom_list.index(i) and 1.2 < dist(i,j) < 1.8:
            ngmap[int(atom_list.index(i))].append(int(atom_list.index(j)))

cent = []
tri = []

##### Put H on the atom 'hydrC' #####

v = int(hydrC) - 1
x0 = float(atom_list[v][1])
y0 = float(atom_list[v][2])
z0 = float(atom_list[v][3])
C = [x0, y0, z0]
xk1 = float(atom_list[ngmap[v][0]][1])
yk1 = float(atom_list[ngmap[v][0]][2])
zk1 = float(atom_list[ngmap[v][0]][3])
xk2 = float(atom_list[ngmap[v][1]][1])
yk2 = float(atom_list[ngmap[v][1]][2])
zk2 = float(atom_list[ngmap[v][1]][3])
xk3 = float(atom_list[ngmap[v][2]][1])
yk3 = float(atom_list[ngmap[v][2]][2])
zk3 = float(atom_list[ngmap[v][2]][3])

Vk1 = [xk1, yk1, zk1]
Vk2 = [xk2, yk2, zk2]
Vk3 = [xk3, yk3, zk3]
Vk = [Vk1, Vk2, Vk3]

V01 = [float(xk1)-float(x0), float(yk1)-float(y0), float(zk1)-float(z0)]
V02 = [float(xk2)-float(x0), float(yk2)-float(y0), float(zk2)-float(z0)]
V03 = [float(xk3)-float(x0), float(yk3)-float(y0), float(zk3)-float(z0)]

cent.append(C) # position of the centre atom
tri.append(Vk) # position of the neighbour 3 atoms

a = 1.10
for i in range(len(cent)):
    cx = (float(tri[i][0][0]) + float(tri[i][1][0]) + float(tri[i][2][0])) / 3   # centre of the triangle
    cy = (float(tri[i][0][1]) + float(tri[i][1][1]) + float(tri[i][2][1])) / 3
    cz = (float(tri[i][0][2]) + float(tri[i][1][2]) + float(tri[i][2][2])) / 3

    cx0 = float(cent[i][0]) - cx   # vector from the triangle centre to central atom -> OK
    cy0 = float(cent[i][1]) - cy
    cz0 = float(cent[i][2]) - cz
    
    lenc0 = (cx0**2 + cy0**2 + cz0**2)**0.5

    Hx = (cx0 / lenc0 * a) + float(cent[i][0])
    Hy = (cy0 / lenc0 * a) + float(cent[i][1])
    Hz = (cz0 / lenc0 * a) + float(cent[i][2])

    hydr = []
    hydr.append('H')
    hydr.append('{:18.011f}'.format(Hx))
    hydr.append('{:18.011f}'.format(Hy))
    hydr.append('{:18.011f}'.format(Hz))
    atom_list.append(hydr)

##### Make a new xyz file #####

print(str(len(atom_list)) + '\n')
for i in atom_list:
    p = [i[0], float(i[1]), float(i[2]), float(i[3])]
    print('\t'.join(str(n) for n in p))

