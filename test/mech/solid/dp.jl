using Amaru
using Test

# mesh 
bls = [
       Block( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2, tag="solids"),
      ]
msh= Mesh(bls, verbose=true)

# fem domain
mats = [ 
       "solids" => DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1)
      ]

dom = Domain(msh, mats)

tag!(dom.elems.ips[1], "ip")
log1 = IpLogger()
loggers = [
           "ip" => log1
          ]
setloggers!(dom, loggers)

# boundary conditions
bcs = [
    :(z==0.0)         => NodeBC(ux=0, uy=0, uz=0),
    :(z==0.5)         => NodeBC(uz=-0.033),
    :(x==0 || x==1.0) => NodeBC(ux=0, uy=0),
    :(y==0 || y==1.0) => NodeBC(ux=0, uy=0),
]

@test solve!(dom, bcs, autoinc=true, nincs=10, tol=1e-2)

# boundary conditions
bcs[2] = :(z==0.5) => NodeBC(uz=+0.008)
@test solve!(dom, bcs, autoinc=true, nincs=10, tol=1e-2)


