using Amaru
using Test

# mesh generation

bl  = Block( [0 0 0; 0.1 0.6 0.1], nx=1, ny=33, nz=6)
msh = Mesh(bl)
iptag!(msh.elems[2], 100)
tag!(msh.elems[22], 100)
tag!(msh.elems, 100)
tag!(msh.edges[:(y==0.2 && z==0.1)], 10)

# finite element analysis

E = 27.e6

mats = [
    MaterialBind(:solids, ElasticSolid(E=E, nu=0.2)),
    MaterialBind(100, SmearedCrack(E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk" ) ),
]

# Loggers
log_edge = EdgeLogger(10)

#model = Model(msh, mats, modeltype="plane-stress", thickness=1.0)
model = Model(msh, mats, log_edge)

# Boundary conditions
bcs = [
       NodeBC(:(y==0 && z==0), :(ux=0, uy=0, uz=0 )),
       NodeBC(:(y==0.6 && z==0), :(ux=0, uz=0)),
       EdgeBC(:(y==0.2 && z==0.1), :(uz=-1e-3)),
       EdgeBC(:(y==0.4 && z==0.1), :(uz=-1e-3)),
      ]

try
    @test solve!(model, autoinc=true, nincs=1000, maxits=4, tol=0.01, report=false, scheme=:FE, nouts=50, maxincs=0).success
catch
end

if @isdefined(makeplots) && makeplots
    using PyPlot
    tab = log_edge.table
    plot(-tab[:uz], -tab[:fz], marker="o", color="blue")
    show()
    #plot(tab[:w], tab[:s1], marker="o", color="blue")
    #show()
end

