using Amaru

# Mesh generation
bl  = Block( [0 0 0; 3.0 3.0 3.0], nx=5, ny=5, nz=5, cellshape=HEX8, tag="solids")
bl1 = BlockInset( [1.5 1.5 3.0; 1.5 1.5 1.5], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl2 = BlockInset( [1.2 1.5 1.6; 0.0 1.5 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl3 = BlockInset( [1.8 1.5 1.6; 3.0 1.5 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl4 = BlockInset( [1.5 1.2 1.6; 1.5 0.0 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl5 = BlockInset( [1.5 1.8 1.6; 1.5 3.0 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bls = [bl, bl1, bl2, bl3, bl4, bl5]

mesh = Mesh(bls)

gw = 10.0     # water specific weight
A  = 0.01     # drain area
k  = 1e-10    # specific permeability
kb = 100
kj = 100
Q  = 4        # volume em metro cubico

# FEM analysis

mats = [
    "solids" << LinSeep(k=k, S=0.0),
    "joints" << SeepJoint1D(A=A) << Joint1DLinSeep(k=kj),
    "drains" << DrainPipe(A=A) << LinDrainPipe(k=kb),
]

ana = HydroAnalysis(gammaw=10)
model = FEModel(mesh, mats, ana)

# Stage 1: pore-pressure stabilization
bcs = [
       :(z==3.0) << NodeBC(uw=0),
      ]
addstage!(model, bcs, tspan=100)
solve!(model, tol=1e-2)


model.env.t = 0.0

tag!(model.elems.lines.nodes[:(x==1.5 && y==1.5 && z==3.0)], "input")

# Stage 2: volume application
bcs = [
       :(x==0.0 && y==1.5 && z==0.0) << NodeBC(uw=0),
       :(x==3.0 && y==1.5 && z==0.0) << NodeBC(uw=0),
       :(x==1.5 && y==0.0 && z==0.0) << NodeBC(uw=0),
       :(x==1.5 && y==3.0 && z==0.0) << NodeBC(uw=0),
       "input" << NodeBC(fw=Q),
      ]
addstage!(model, bcs, tspan=100)
solve!(model, tol=1e-2)
