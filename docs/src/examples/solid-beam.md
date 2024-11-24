# Example of a solid beam using Von-Mises model

A finite element analysis of a cantilever beam fixed at the left end and a prescribed displacement at the right end.

```@example
using Amaru

# Data
th = 0.05  # thickness (m)
L  = 1.0   # beam length (m)
h  = 0.1   # beam height (m)
E  = 210e6 # Young modulus (kPa)
nu = 0.3   # Poisson ratio
fy = 240e3 # Yied strength (kPa)
H  = 0.0   # Hardening modulus (kPa)

# Geometry model
geo = GeoModel()
sz = 0.05 # mesh size

# List of points
p1 = addpoint!(geo, 0, 0, 0, size=sz)
p2 = addpoint!(geo, L, 0, 0, size=sz)
p3 = addpoint!(geo, L, 0, h/2, size=sz)
p4 = addpoint!(geo, L, 0, h, size=sz)
p5 = addpoint!(geo, 0, 0, h, size=sz)

# Define a closed loop to create a surface
addpath!(geo, :M, p1, :L, p2, :L, p3, :L, p4, :L, p5, :L, p1)

# Pull the surface to create a solid
pull!(geo, geo.surfaces, axis=[0, 1, 0], length=th)

# Finite element mesh
mesh= Mesh(geo)

# List of element types and constitutive model
mat = [ :all => MechSolid => VonMises => (E=E, nu=nu, fy=fy, H=H) ]

# A mechanical analysis context
ctx = MechContext()

# A finite element model
model = FEModel(mesh, mat, ctx)

# A finite element analysis object
ana = MechAnalysis(model)

# List of data loggers
loggers = [
    (x==L, z==h/2) => NodeSumLogger("file.dat")
]
addloggers!(ana, loggers)

# List of monitors
monitors = [
    (x==L, y==0, z==h/2) => NodeMonitor(:fz)
]
addmonitors!(ana, monitors)

# List of boundary conditions
bcs = [
    x==0 => NodeBC(ux=0, uy=0, uz=0),
    (x==L, z==h/2) => NodeBC(uz=-0.08),
]

# Adds a load stage
addstage!(ana, bcs, nincs=20, nouts=1)

# Run the analysis
solve!(ana, autoinc=true)
```
