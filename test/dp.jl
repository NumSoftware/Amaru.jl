using Amaru
using Base.Test

# mesh 
bls = [
       Block3D( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2),
      ]
msh= Mesh(bls, verbose=true)

# fem domain
mat = MaterialBind(:all, DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1), iptag="ip" )

mon = Logger(:ip, "ip")

# boundary conditions
bcs = [
    BC(:node, :(z==0.0) , ux=0, uy=0, uz=0),
    BC(:node, :(z==0.5) , uz=-0.033),
    BC(:node, :(x==0 || x==1.0), ux=0, uy=0),
    BC(:node, :(y==0 || y==1.0), ux=0, uy=0),
]

dom = Domain(msh, mat, mon)
@test solve!(dom, bcs, autoinc=true, nincs=10, tol=1e-2)

# boundary conditions
top    = BC(:node, :(z==0.5) , uz=+0.008)
bcs[2] = top
@test solve!(dom, bcs, autoinc=true, nincs=10, tol=1e-2)


