using Amaru
using Test

# mesh
bls = [
       Block( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2, tag="solids"),
      ]
msh= Mesh(bls, silent=true)

# fem domain
mats = [
        "solids" => VonMises(E=2.0e8, nu=0.28, σy=5.0e5)
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
    :(z==0)           => NodeBC(uz=0),
    :(z==0.5)         => NodeBC(uz=+0.033),
    :(x==0 || x==1.0) => NodeBC(ux=0),
    :(y==0 || y==1.0) => NodeBC(uy=0),
]

@test solve!(dom, bcs, autoinc=true, nincs=40, tol=1e-2).success

if Amaru.config.makeplots
    using PyPlot
    tab = log1.table
    plot( tab[:ezz], tab[:szz], "-o")
    show()
    plot( tab[:j1], tab[:srj2d], "-o")
    show()
end

