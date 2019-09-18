using Amaru

# Mesh generation
bl  = Block( [0 0 0; 3.0 3.0 3.0], nx=5, ny=5, nz=5, cellshape=HEX8, tag="solids")
bl1 = BlockInset( [1.5 1.5 3.0; 1.5 1.5 1.5], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl2 = BlockInset( [1.2 1.5 1.6; 0.0 1.5 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl3 = BlockInset( [1.8 1.5 1.6; 3.0 1.5 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl4 = BlockInset( [1.5 1.2 1.6; 1.5 0.0 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bl5 = BlockInset( [1.5 1.8 1.6; 1.5 3.0 0.0], curvetype="polyline", cellshape=LIN3, tag="drains", jointtag="joints")
bls = [bl, bl1, bl2, bl3, bl4, bl5]

mesh = Mesh(bls, verbose=true)

gw = 10.0     # water specific weight
A  = 0.01     # drain area
k  = 1e-10    # specific permeability
kb = 100 
kj = 100 
Q  = 4        # volume em metro cubico

# FEM analysis

mats = [
    "solids" => LinSeep(k=k, gammaw=gw, S=0.0),
    "joints" => Joint1DLinSeep(k=kj, gammaw=gw, A=A),
    "drains" => LinDrainPipe(k=kb, gammaw=gw, A=A),
]

dom = Domain(mesh, mats,gammaw=10)

# Stage 1: pore-pressure stabilization
bcs = [
       :(z==3.0) => NodeBC(uw=0),
      ]

hm_solve!(dom, bcs, end_time=100.0, nincs=1, tol=1e-2, nouts=0, verbose=true)


dom.env.t = 0.0

tag!(dom.elems[:lines][:nodes][:(x==1.5 && y==1.5 && z==3.0)], "input")

# Stage 2: volume application
bcs = [
       :(x==0.0 && y==1.5 && z==0.0) => NodeBC(uw=0),
       :(x==3.0 && y==1.5 && z==0.0) => NodeBC(uw=0),
       :(x==1.5 && y==0.0 && z==0.0) => NodeBC(uw=0),
       :(x==1.5 && y==3.0 && z==0.0) => NodeBC(uw=0),
       "input" => NodeBC(fw=Q),
      ]


hm_solve!(dom, bcs, end_time=100.0, nincs=1, tol=1e-2, nouts=0, verbose=true)
