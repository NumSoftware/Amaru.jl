using Amaru
using Test

# mesh
bls = [
       Block( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2, tag="solids"),
      ]
msh= Mesh(bls)

# fem domain
mats = [
       "solids" << DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1)
      ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)
 
tag!(model.elems.ips[1], "ip")
log1 = IpLogger()
loggers = [
           "ip" << log1
          ]
setloggers!(model, loggers)

# boundary conditions
bcs = [
    :(z==0.0)         << NodeBC(ux=0, uy=0, uz=0),
    :(z==0.5)         << NodeBC(uz=-0.033),
    :(x==0 || x==1.0) << NodeBC(ux=0, uy=0),
    :(y==0 || y==1.0) << NodeBC(ux=0, uy=0),
]
addstage!(model, bcs, nincs=10)
bcs[2] = :(z==0.5) << NodeBC(uz=+0.008)
addstage!(model, bcs, nincs=10)

@test solve!(model, tol=1e-2, autoinc=true).success

# boundary conditions
# bcs[2] = :(z==0.5) << NodeBC(uz=+0.008)
# @test solve!(model, autoinc=true, nincs=10, tol=1e-2).success


