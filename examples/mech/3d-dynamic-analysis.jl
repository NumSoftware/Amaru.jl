using Amaru

# Mesh generation
blocks = [
    Block( [0 0 0; 0.2 2.0 0.2], nx=2, ny=12, nz=2, shape=HEX8, tag="solids"),
]

msh = Mesh(blocks)

# Finite element modeling
materials = [
             "solids" => MechSolid => LinearElastic => (E=36e6, nu=0.2, rho=24.0),
            ]

# Finite element model
ana = MechAnalysis()
model = Model(msh, materials, ana)

loggers = [
           :(x==0.2 && y==2 && z==0.2) => NodeLogger("end-node.dat"),
          ]

setloggers!(model, loggers)

# List of boundary conditions
bcs = [
       :(y==0 && z==0)   => NodeBC(ux=0, uy=0, uz=0),
       :(y==2 && z==0)   => NodeBC(uz=0),
       :(y==1 && z==0.2) => NodeBC(fz=-10),
]

addstage!(model, bcs, nincs=10)

# Perform the finite element analysis
solve!(model, tspan=1.0, nincs=4)
