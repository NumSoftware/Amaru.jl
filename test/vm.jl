using Amaru
using Base.Test

# mesh 
bls = [
       Block3D( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2),
      ]
msh= Mesh(bls, verbose=true)
iptag!(msh.cells[end], "ip")

# fem domain
mat = MaterialBind(:all, VonMises(E=2.0e8, nu=0.28, σy=5.0e5) )

mon = Logger(:ip, "ip")

# boundary conditions
bcs = [
    #BC(:node, :(x==0 && y==0 && z==0) , ux=0, uy=0, uz=0),
    BC(:node, :(z==0), uz=0),
    BC(:node, :(z==0.5) , uz=+0.033),
    BC(:node, :(x==0 || x==1.0), ux=0),
    BC(:node, :(y==0 || y==1.0), uy=0),
]

dom = Domain(msh, mat, mon)
@test solve!(dom, bcs, auto_inc=true, nincs=40, tol=1e-2)

# boundary conditions
#top    = BC(:node, :(z==0.5) , uz=+0.008)
#bcs[2] = top
#@test solve!(dom, bcs, auto_inc=true, nincs=10, tol=1e-2)


tab = mon.table

using PyPlot
plot( tab[:ezz], tab[:szz], "-o")
show()
plot( tab[:j1], tab[:srj2d], "-o")
show()

