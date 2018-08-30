using Amaru
using Test

bl  = Block3D( [0 0 0; 1. 1. 1.], nx=5, ny=5, nz=5, shape=HEX8)

# mesh generation
msh = Mesh(bl, verbose=false)

mat = MaterialBind(:solids, Mazars(E=30000, nu=0.2, eps0=1.e-4, At=0.9, Bt=5000., Ac=1.0, Bc=1500.0))

dom = Domain(msh, mat) 

tag!(dom.elems[:ips][1], "ip")
#logger = IpLogger("ip")
logger = FaceLogger(:(z==1))
setlogger!(dom, logger)

bcs = [
       NodeBC(:(x==0), ux=0)
       NodeBC(:(y==0), uy=0)
       NodeBC(:(z==0), uz=0)

       NodeBC(:(z==0), ux=0, uy=0)

       #NodeBC(:(z==1), uz=-1.2e-2)
       FaceBC(:(z==1), uz=+2e-3)
      ]

@test solve!(dom, bcs, autoinc=true, nincs=100, maxits=4, tol=0.002, nouts=10, verbose=true)

if Amaru.Debug.makeplots
    using PyPlot
    tab = logger.table
    #plot(tab[:ezz], tab[:szz], marker="o")
    plot(tab[:uz], tab[:fz], marker="o")
    show()
end
