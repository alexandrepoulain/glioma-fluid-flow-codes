import sys 
from dolfin import *
import time
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('--mesh',default="DTI_16", type=str)      
parser.add_argument('--time_steps',default=90, type=int)      
parser.add_argument('--final_time', default=9, type=float) # default is 9 hours :w
parser.add_argument('--lumped',default="lumped", type=str)      
parser.add_argument('--label', default="std", type=str)
args = parser.parse_args()

numpy.set_printoptions(threshold=sys.maxsize)

# read mesh and subdomain data 
# We assume the following subdomains: 1:gray, 2:white matters and 4:tumor
# marked as in all.sh 
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(),args.mesh, "r")
hdf.read(mesh, "/mesh", False)  
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")

boundaries  = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")
dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# print vertices, cells and volume 
print ("mesh num vertices ", mesh.num_vertices())
print ("mesh num cells", mesh.num_cells())
brain_volume = assemble(Constant(1)*dx)  
print( "Brain volume (mm)^3 ", brain_volume) 

# create the diffusion tensor (constants on each element) 
T = TensorFunctionSpace(mesh, "DG", 0)
D = Function(T) 

# determine the diffusion coefficients (by "hand" since we do not access to DTI for Patient 1)
coeff_diffusion_wm = 1.4e-4
coeff_diffusion_gm = 1.6e-4
coeff_diffusion_tumor = 0.33 * coeff_diffusion_wm

# various parameters   
time_steps = args.time_steps 
time_final = args.final_time    
dt         = 3600*(time_final / time_steps) # units are mm and s
u_0 = Constant(0.0) 
f   = Constant(0.0)

# output file 
vtkfile = File('diffusion-mritracer-{}/solution.pvd'.format(args.label) )

# function space for solution 
V    = FunctionSpace(mesh, 'Lagrange', 1)

# bc
ud_str = 't < 3600*10 ? 1.0 : 0.0'   
u_d    = Expression(ud_str,degree=1,t=dt)
bc   = DirichletBC(V, u_d, 'on_boundary')

# previous time step 
u_n = interpolate(u_0, V)
bc.apply(u_n.vector())

# Rename u_n to Concentration and have unit M (Molar)
u_n.rename("Concentration","M")

# dump initial condition 
vtkfile << (u_n,0.00) 

# trial and test function 
u = TrialFunction(V)
v = TestFunction(V)


# in case of mass lumping a separate lumped mass matrix 
M = u_n*dx
mass_form = v*u*dx
Mass_matrix = assemble(mass_form) # consistent mass matrix

# the variational forms (without the mass matrix) 
F = dt*coeff_diffusion_gm*dot(grad(u),grad(v))*dx(1)+    dt*coeff_diffusion_wm*dot(grad(u), grad(v) )*dx(2)  + \
     dt*coeff_diffusion_tumor*dot(grad(u), grad(v)  )*dx(4) - \
 (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# (Re-)define u as the solution at the current time 
u = Function(V)    # AU

# Rename u_n to Concentration and have unit M (Molar)
u_n.rename("Concentration","M")

# Create system matrix 
A = assemble(a)
A = A+Mass_matrix 

# compute the volumes of the various subdomains 
vol1 = assemble(1.0*dx(1)) 
vol2 = assemble(1.0*dx(2))
vol4 = assemble(1.0*dx(4))
print( "Volume of domains gray matter, white matter, tumor: ", vol1 ,vol2,vol4) 

# initialize lists for the concentrations of the different subdomains 
unit1   = []
unit2 = []
unit4 = []

# time loop 
time_points = []
t = 0
for n in range(time_steps):
   # update time
   t += dt
   u_d.t=t
   # solve the system 
   b = assemble(L)
   bc.apply(A,b)
   solve(A,u.vector(),b, "gmres", "amg")

   # store results every hour 
   if t/(60*60) - int(t/(60*60)) < dt/(60*60):  vtkfile << (u, t)

   # update previous solution 
   u_n.assign(u)

   # compute the amount in the various subdomains 
   unit1 += [assemble(u*dx(1))/vol1]
   unit2 += [assemble(u*dx(2))/vol2]
   unit4 += [assemble(u*dx(4))/vol4]

   time_points += [t] 

   # print to screen 
   print( "time [hours]", t/(60*60), " concentrations at 1 2 4 are: ",  unit1[-1], unit2[-1], unit4[-1])  


# write the various results to csv files 
import csv
with  open('time_{}.csv'.format(args.label), 'w') as outfile:
   writer = csv.writer(outfile)
   writer.writerows(map(lambda x: [x], time_points))

with  open('tracer1_{}.csv'.format(args.label), 'w') as outfile:
   writer = csv.writer(outfile)
   writer.writerows(map(lambda x: [x], unit1))

with  open('tracer2_{}.csv'.format(args.label), 'w') as outfile:
   writer = csv.writer(outfile)
   writer.writerows(map(lambda x: [x], unit2))

with  open('tracer4_{}.csv'.format(args.label) , 'w') as outfile:
   writer = csv.writer(outfile)
   writer.writerows(map(lambda x: [x], unit4))

















