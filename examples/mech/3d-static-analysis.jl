using Amaru

# Mesh generation
blocks = [
    Block( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, cellshape=HEX8, tag="solids"),
]

msh = Mesh(blocks);


# Finite element modeling
materials = [
             "solids" => MechSolid => LinearElastic => (E=100.0, nu=0.2),
            ]

# Finite element model
ana = MechAnalysis()
model = Model(msh, materials, ana)
addlogger!(model, :(z==1) => FaceLogger("top-face.dat"))

# List of boundary conditions
bcs = [
       :(z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(z==1) => SurfaceBC(tz=:(-10*x)),   # triangular load
]
addstage!(model, bcs, nincs=4, nouts=1)

# Perform the finite element analysis
solve!(model)