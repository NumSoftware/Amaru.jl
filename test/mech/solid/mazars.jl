using Amaru
using Base.Test

bl  = Block3D( [0 0 0; 1. 1. 1.], nx=1, ny=1, nz=1, shape=HEX8)

# mesh generation
msh = Mesh(bl, verbose=false)

mat = MaterialBind(:solids, Mazars(E=30000, nu=0.2, eps0=1.e-4, At=0.9, Bt=5000., Ac=1.0, Bc=1500.0))

dom = Domain(msh, mat) 

tag!(dom.elems[:ips][1], "ip")
logger = Logger(:ip, "ip")
setlogger!(dom, logger)

bcs = [
       BC(:node, :(z==0), ux=0, uy=0, uz=0 )
       BC(:face, :(z==1), uz=-1.5e-2)
       #BC(:face, :(z==1), uz=+2e-3)
      ]

@test solve!(dom, bcs, autoinc=true, nincs=20, maxits=4, tol=0.01, verbose=true)

if Amaru.Debug.makeplots
    using PyPlot
    tab = logger.table
    plot(tab[:ezz], tab[:szz], marker="o", color="blue")
    show()
end
