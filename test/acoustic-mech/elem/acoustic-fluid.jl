using Amaru

# Mesh generation
L = 4.0
c = 1500

bl   = Block( [0 0; L 1], cellshape=QUAD4, nx=4, ny=2, tag="fluido" )
mesh = Mesh(bl)

# Finite element analysis
mats = [
    "fluido" => AcousticFluid => LinearAcousticFluid => (c=c,)
]

ana = AcousticMechAnalysis(thickness=1.0)
model = FEModel(mesh, mats, ana)

addlogger!(model, :(x==4 && y==0.5) => NodeLogger("node.dat") )
addlogger!(model, :(y==0.5) => NodeGroupLogger("nodes.dat") )

t0 = 0.001
q0 = 0.0001  # kN/m/m2, kg/s2/m2 (por superficie)

tmax = 0.005
# q = :( t<=$tmax/5 ? $q0 : 0.0 )
q = :( t<=$tmax/5 ? $q0*t : $q0*$tmax/5 )

bcs = [
    :(x==0) => SurfaceBC(tq=q)   # tq.  [kN/m/m2] [kg/s2/m2]
]

addstage!(model, bcs, nincs=20, nouts=20, tspan=tmax)
solve!(model, quiet=false, tol=1e-1)