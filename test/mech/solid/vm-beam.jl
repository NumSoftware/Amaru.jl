using Amaru
using Base.Test

# mesh 
bls = [
       Block3D( [0 0 -0.05; 0.05 1.0 0.05], nx=1, ny=50, nz=2, shape=HEX20),
      ]
msh= Mesh(bls, verbose=true)
iptag!(msh.cells[end], "ip")

# fem domain
mat = MaterialBind(:all, VonMises(E=210e6, nu=0.3, σy=0.24e6) )
#mat = MaterialBind(:all, DruckerPrager(E=210e6, nu=0.3, alpha=0.05e6, kappa=0.1 ) )

mon = Logger(:ip, "ip")


# boundary conditions
bcs = [
    BC(:node, :(y==0), :(uy=0)),
    BC(:node, :(y==0 && z==0), :(uz=0)),
    BC(:node, :(x==0 && y==0 && z==0), :(ux=0)),
    BC(:edge, :(y==1 && z==0), :(uz=-0.1)),
]

dom = Domain(msh, mat)
mon = Logger(dom.edges[:(y==1 && z==0)])
setlogger!(dom, mon)

@test solve!(dom, bcs, autoinc=true, nincs=6, nouts=1, tol=1e-2)

# boundary conditions
#top    = BC(:node, :(z==0.5) , uz=+0.008)
#bcs[2] = top
#@test solve!(dom, bcs, autoinc=true, nincs=10, tol=1e-2)


if Amaru.Debug.makeplots
    using PyPlot
    tab = mon.table
    plot( tab[:uz], tab[:fz], "-o")
    show()
end

