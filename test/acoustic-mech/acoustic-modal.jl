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

ctx = AcousticContext()
model = FEModel(mesh, mats, ctx, thickness=1.0)
ana = AcousticModalAnalysis(model)

bcs = [
    :(x==0) => SurfaceBC(up=0)   # tq.  [kN/m/m2] [kg/s2/m2]
]

addstage!(ana, bcs)

solve!(ana, nmodes=10, quiet=false)