using Amaru
using Test

# mesh
bls = [
       Block( [0 0 0; 0.1 0.1 0.2], nx=2, ny=2, nz=2, tag="solids"),
      ]
msh= Mesh(bls)

# fem domain
mats = [
       "solids" => MechSolid => WillamWarnke => (E=30e6, nu=0.25, fc=-30e3, ft=3e3, fb=-33e3, H=0)
      ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)
 
loggers = [
           [0, 0, 0.1] => IpLogger()
          ]
setloggers!(model, loggers)

# boundary conditions
bcs = [
    # and(x==0, y==0, z==0) => NodeBC(ux=0, uy=0),
    x==0 => NodeBC(ux=0),
    y==0 => NodeBC(uy=0),
    z==0 => NodeBC(uz=0),
    z==0.2 => NodeBC(uz=-0.001),
    x==0.1 => NodeBC(ux=-0.0005),
    # y==0.1 => NodeBC(uy=-0.001),
    # z==0.2 => NodeBC(uz=+0.0001),
]
addstage!(model, bcs, nincs=10, nouts=10)

# bcs[2] = :(z==0.5) => NodeBC(uz=+0.001)
# addstage!(model, bcs, nincs=10)
solve!(model, tol=1e-1, autoinc=true, quiet=false).success

# boundary conditions
# bcs[2] = :(z==0.5) => NodeBC(uz=+0.008)
# @test solve!(model, autoinc=true, nincs=10, tol=1e-2).success


# Plotting

chart = Chart(
    xlabel = L"\varepsilon_{zz}\times 1000",
    ylabel = L"\sigma_{zz}",
    xmult = 1e3,
)
table = model.loggers[1].table

series = [ 
    LineSeries(table.ezz, table.szz),
]

addseries!(chart, series)
save(chart, "chart.pdf")