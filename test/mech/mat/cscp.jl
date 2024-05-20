using Amaru
using Test

# mesh
bls = [
       Block( [0 0 0; 0.1 0.1 0.2], nx=1, ny=1, nz=2, tag="solids"),
      ]
msh= Mesh(bls)

fc = -30e3
ft = 3e3
fb = -33e3

# fc_fun = PathFunction(:M, 0.0, 0.4*fc, :L, 0.01, 1.0*fc)
fc_fun = PathFunction(:M, 0.0, 0.4*fc, :Q, 0, 0.7*fc, 0.02, fc)
ft_fun = PathFunction(:M, 0.0, 1.0*ft, :L, 0.02, ft)

fic_fun = PathFunction(:M, 0.0, fc, :L, 0.002, 0.8*fc)

# fem domain
mats = [
       "solids" => MechSolid => CSCP => (E=30e6, nu=0.25, fc=fc, ft=ft, fb=fb,  fc_fun=fc_fun, ft_fun=ft_fun, fic_fun=fic_fun)
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
    z==0.2 => NodeBC(uz=-0.0005),
    # x==0.1 => NodeBC(ux=-0.0005),
    # y==0.1 => NodeBC(uy=-0.001),
    # z==0.2 => NodeBC(uz=+0.0001),
]
addstage!(model, bcs, nincs=10, nouts=20)

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
    LineSeries(-table.ezz, -table.szz),
]

addseries!(chart, series)
save(chart, "chart.pdf")