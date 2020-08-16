using Amaru
using Test

bl  = Block( [0 0 0; 1. 1. 1.], nx=2, ny=2, nz=2, cellshape=HEX8, tag="solids")

# mesh generation
msh = Mesh(bl, silent=true)

mats = [
        "solids" => Mazars(E=30000, nu=0.2, eps0=1.e-4, At=0.9, Bt=5000., Ac=1.0, Bc=1500.0)
       ]

dom = Domain(msh, mats, silent=true)


tag!(dom.elems.ips[1], "ip")
log1 = IpLogger()
log2 = FaceLogger()
loggers = [
           "ip" => log1
           :(z==1) => log2
          ]
setloggers!(dom, loggers)

bcs = [
       :(x==0) => NodeBC(ux=0),
       :(y==0) => NodeBC(uy=0),
       :(z==0) => NodeBC(uz=0),

       :(z==0) => NodeBC(ux=0, uy=0),

       #:(z==1) => NodeBC(uz=-1.2e-2),
       :(z==1) => FaceBC(uz=+2e-3),
      ]

@test solve!(dom, bcs, autoinc=true, nincs=100, maxits=4, tol=0.1, nouts=10, verbose=false).success

if Amaru.config.makeplots
    using PyPlot
    tab = log2.table
    #plot(tab[:ezz], tab[:szz], marker="o")
    plot(tab[:uz], tab[:fz], marker="o")
    show()
end
