using Amaru
using Test

coords = [ 0.0; 1.0 ]
conns  = [ [ 1, 2 ], [2] ]

msh = Mesh(coords, conns)


mat = [
       1 => ElasticSpring(kx=100.0) ,
       2 => LumpedMass(m=100.0) ,
      ]

model = Model(msh, mat)

log1 = NodeLogger()
logs = [
        :(x==1) => log1
       ]
setloggers!(model, logs)


bcs = [
       :(x==0) => NodeBC(ux=0),
       :(x==1) => NodeBC(fx=1),
       #:(x==1) => NodeBC(fx=:( t<0.1 ? -1 : 0 )),
      ]

addstage!(model, bcs, tspan=7, nincs=14)

dyn_solve!(model; quiet=false)


if @isdefined(makeplots) && makeplots
    using PyPlot
    t = log1.table
    plot(t[:t], t[:ux], "-o")
    show()
end



