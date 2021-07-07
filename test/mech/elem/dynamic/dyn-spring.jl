using Amaru
using Test

#coords = [ 1.0; 0.0 ]
#conns  = [ [ 1, 2 ], [1] ]
coords = [ 0.0; 1.0 ]
conns  = [ [ 1, 2 ], [2] ]

msh = Mesh(coords, conns)


mat = [
       1 => ElasticSpring(kx=100.0) ,
       2 => LumpedMass(m=100.0) ,
      ]

dom = Domain(msh, mat)

log1 = NodeLogger()
logs = [
        :(x==1) => log1
       ]
setloggers!(dom, logs)


bcs = [
       :(x==0) => NodeBC(ux=0),
       :(x==1) => NodeBC(fx=1),
       #:(x==1) => NodeBC(fx=:( t<0.1 ? -1 : 0 )),
      ]

#solve!(dom, bcs)
@test dynsolve!(dom, bcs, time_span=7, nincs=14, verbosity=0)

log1.table;

if Amaru.config.makeplots
    using PyPlot
    t = log1.table
    plot(t[:t], t[:ux], "-o")
    show()
end



