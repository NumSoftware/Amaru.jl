using Amaru
using Test

# mesh
bls = [
       Block( [0 0 -0.05; 0.05 1.0 0.05], nx=1, ny=50, nz=2, cellshape=HEX20),
      ]
msh= Mesh(bls)
iptag!(msh.elems[end], "ip")

# fem domain
mat = MaterialBind(:all, VonMises << (E=210e6, nu=0.3, fy=0.24e6) )
#mat = MaterialBind(:all, DruckerPrager << (E=210e6, nu=0.3, alpha=0.05e6, kappa=0.1 ) )

mon = NodeLogger("ip")


# boundary conditions
bcs = [
    NodeBC(:(y==0), :(uy=0)),
    NodeBC(:(y==0 && z==0), :(uz=0)),
    SurfaceBC(:(x==0 && y==0 && z==0), :(ux=0)),
    SurfaceBC(:(y==1 && z==0), :(uz=-0.1)),
]

model = FEModel(msh, mat)
mon = NodeLogger(model.edges[:(y==1 && z==0)])
setlogger!(model, mon)

@test solve!(model, autoinc=true, nincs=6, nouts=1, tol=1e-2).success

if @isdefined(makeplots) && makeplots
    using PyPlot
    tab = mon.table
    plot( tab[:uz], tab[:fz], "-o")
    show()
end

