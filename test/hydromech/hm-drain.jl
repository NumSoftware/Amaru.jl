using Amaru

# Mesh generation
bl  = Block( [0 0; 2.0 2.0], nx=2, ny=2, cellshape=QUAD8, tag="solids")
bl1 = BlockInset( [1.5 2.0; 1.5 1.5], curvetype="polyline", tag="drains", cellshape=LIN2, jointtag="joints") 
bl2 = BlockInset( [1.5 1.4; 0.001 0.001], curvetype="polyline", tag="drains", cellshape=LIN2, jointtag="joints")
bls = [bl, bl1, bl2]

mesh = Mesh(bls)

E  = 5000;
nu = 0.25;
gw = 10.0     # water specific weight
A  = 0.01
p = 2*pi*(A/pi)^0.5
k  = 1e-4
kb = 1 
kj = 100 
Q  = 1       # volume em metro cubico

# FEM analysis

mats = [
    "solids" => HMSolid => LinearElasticSeep => (E=E, nu=nu, k=k, alpha=1.0, S=0.0),
    "joints" => SeepJoint1D => Joint1DConstPermeability => (k=kj, p=p),
    "drains" => DrainPipe => LinDrainPipe => (k=kb, A=A),
]

ana = HydromechAnalysis(gammaw=10)
model = FEModel(mesh, mats, ana)

changequadrature!(model.elems.lines, 3)

# Stage 1: pore-pressure stabilization
bcs = [
       :(y==0.0) => NodeBC(ux=0,uy=0),
       :(y==2.0) => NodeBC(uw=0),
      ]
addstage!(model, bcs, tspan=100, nincs=2, nouts=2)
solve!(model, tol=1e-2)

model.env.t = 0.0

# Stage 2: volume application
bcs = [
       :(y==0.0) => NodeBC(ux=0,uy=0),
       :(x==0.0 && y==0.0) => NodeBC(uw=0),
       :(x==1.5 && y==2.0) => NodeBC(fw=Q),
      ]

addstage!(model, bcs, tspan=600, nincs=2, nouts=2)
solve!(model, tol=1e-2)
