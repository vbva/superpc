from math import pi, sin, exp, sqrt
import numpy as np
import meshio

size = 4 # number of out_*.txt files

with open("out_0.txt") as f:
    cell = f.readline().split()


Nx = int(cell[0])
Ny = int(cell[1])
Nz = int(cell[2])


hx = 1 / (Nx - 1)
hy = 1 / (Ny - 1)
hz = 1 / (Nz - 1)

u = []
x = []
y = []
z = []
u_an = []
    
for j in range(0, size):
    filename = "out_%d.txt" % j
    u_j = list(np.loadtxt(filename, skiprows=1)[:, 0]) # U
    x_j = list(np.loadtxt(filename, skiprows=1)[:, 1]) # X
    y_j = list(np.loadtxt(filename, skiprows=1)[:, 2]) # Y
    z_j = list(np.loadtxt(filename, skiprows=1)[:, 3]) # Z
    
    u += u_j
    x += x_j
    y += y_j
    z += z_j

    u_an_j = []
    diff = []
    diff_sum = 0

    # Give a list of analytical solutions
    for i in range(0, len(u_j)):
        n = sin(pi * x_j[i] * hx) * sin(pi * y_j[i] * hy) * sin(pi * z_j[i] * hz) * (1 - exp(-0.5 * pi*pi))
        u_an_j.append(n)
    u_an += u_an_j
    # Give a list of difference between analytical solutions and list of U
    for i in range(0, len(u_j)):
        diff.append(abs(u_an_j[i] - u_j[i]))

    #Sum of all differences
    for i in range(0, len(diff)):
        diff_sum += diff[i]

    #standard deviation
    dev = sqrt(diff_sum * diff_sum / len(diff))


    print("Max diff in %s =" % filename, max(diff))
    print("Sum of diff in %s =" % filename, diff_sum)
    print("Standart deviation in %s =" % filename, round(dev,2))
    print("End if %s\n" % filename)

zipped_list = list(zip(z, y, x, u, u_an))

for i in range(0, Nx):
    for j in range(0, Ny):
        for k in range(0, Nz):
            if (i == 0 or j == 0 or k == 0 or i == (Nx - 1) or j ==(Ny - 1) or k == (Nz - 1)):
                zipped_list.append((k, j, i, 0, 0))

zipped_list.sort()

z, y, x, u, u_an = list(zip(*zipped_list))

""""""""""""""""""""""""""""""""""""""""""
diff = []
diff_sum = 0

# Give a list of difference between analytical solutions and list of U
for i in range(0, len(u)):
    diff.append(abs(u_an[i] - u[i]))

#Sum of all differences
for i in range(0, len(diff)):
    diff_sum += diff[i]

#standard deviation
dev = sqrt(diff_sum * diff_sum / len(diff))


print("Max diff =", max(diff))
print("Sum of diff =", diff_sum)
print("Standart deviation =", round(dev,2))

""""""""""""""""""""""""""""""""""""""""""

points_all = []
for i in range(0, len(x)):
    points_all.append(x[i] * hx)
    points_all.append(y[i] * hy)
    points_all.append(z[i] * hz)
points = [points_all[i:i+3] for i in range(0, len(points_all), 3)]

cells = []

for i in range(0, Nx - 1):
    for j in range(0, Ny - 1):
        for k in range(0, Nz - 1):
            cells_list = []
            cells_list.append(i + 1 + j * Nx + k * Nx * Ny)
            cells_list.append(i + 1 + (j + 1) * Nx + k * Nx * Ny)
            cells_list.append(i + 1 + (j + 1) * Nx + (k + 1) * Nx * Ny)
            cells_list.append(i + 1 + j * Nx + (k + 1) * Nx * Ny)
            cells_list.append(i + j * Nx + k * Nx * Ny)
            cells_list.append(i + (j + 1) * Nx + k * Nx * Ny)
            cells_list.append(i + (j + 1) * Nx + (k + 1) * Nx * Ny)
            cells_list.append(i + j * Nx + (k + 1) * Nx * Ny)
            cells += [( "hexahedron", [cells_list])]
            
#cells = [("quad8", [[22, 1, 0, 21, 463, 442, 441, 462]])]
mesh = meshio.Mesh(
    points,
    cells,
    point_data={"U": u},
)

mesh.write("numeric.vtk", file_format="vtk42")

mesh = meshio.Mesh(
    points,
    cells,
    point_data={"U": u_an},
)

mesh.write("analytic.vtk", file_format="vtk42")
