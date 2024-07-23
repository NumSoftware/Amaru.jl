using Amaru
using Test

# Mesh
bl = Block([0 0; 2 1], tag="fluid", nx=20, ny=10, shape=QUAD4)
mesh = Mesh(bl)

# Finite element model
mats = [ "fluid" => MechFluid => LinearElasticFluid => (K=2.2e6, rho=1.0, gamma=10.0) ]  # K in kPa

tmax = 10
wy = :(-10*t/$tmax)


bcs =
    [
     y==0 => NodeBC(uy=0),
     x==0 => NodeBC(ux=0),
     x==2 => NodeBC(ux=0),
     "fluid" => BodyC(wy=wy)
    ]

ana = DynAnalysis()
model = FEModel(mesh, mats, ana)

addstage!(model, bcs, nincs=10, nouts=10, tspan=tmax)

solve!(model, quiet=false)

