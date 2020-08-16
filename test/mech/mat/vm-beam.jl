using Amaru
using Test

# mesh
bls = [
       Block( [0 0 -0.05; 0.05 1.0 0.05], nx=1, ny=50, nz=2, cellshape=HEX20),
      ]
msh= Mesh(bls, silent=true
iptag!(msh.elems[end], "ip")

# fem domain
mat = MaterialBind(:all, VonMises(E=210e6, nu=0.3, σy=0.24e6) )
#mat = MaterialBind(:all, DruckerPrager(E=210e6, nu=0.3, alpha=0.05e6, kappa=0.1 ) )

mon = NodeLogger("ip")


# boundary conditions
bcs = [
    NodeBC(:(y==0), :(uy=0)),
    NodeBC(:(y==0 && z==0), :(uz=0)),
    FaceBC(:(x==0 && y==0 && z==0), :(ux=0)),
    FaceBC(:(y==1 && z==0), :(uz=-0.1)),
]

dom = Domain(msh, mat)
mon = NodeLogger(dom.edges[:(y==1 && z==0)])
setlogger!(dom, mon)

@test solve!(dom, bcs, autoinc=true, nincs=6, nouts=1, tol=1e-2, verbose=false).success

if Amaru.config.makeplots
    using PyPlot
    tab = mon.table
    plot( tab[:uz], tab[:fz], "-o")
    show()
end

