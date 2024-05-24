using Amaru

# Mesh generation
L = 4.0
c = 1500

bl   = Block( [0 0; L 1], cellshape=QUAD4, nx=4, ny=2, tag="fluido" )
mesh = Mesh(bl)

# Finite element analysis
mats = [
    "fluido" => AcousticFluid => LinearAcousticFluid => (rho=1000.0, c=c)
]

ana = AcousticModalAnalysis(thickness=1.0)
model = FEModel(mesh, mats, ana)

bcs = [
    :(x==0) => SurfaceBC(up=0)   # tq.  [kN/m/m2] [kg/s2/m2]
]

addstage!(model, bcs)

solve!(model, nmodes=10, quiet=false)