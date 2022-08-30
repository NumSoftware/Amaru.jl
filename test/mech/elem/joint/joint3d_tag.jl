using Amaru
using Test

# mesh generation

bl1  = Block( [0 0 0; 0.1 0.1 0.1] , nx=1, ny=1, nz=1, cellshape=HEX8, tag="a")
bl2  = Block( [0 0.1 0;0.1 0.2 0.1], nx=1, ny=1, nz=1, cellshape=HEX8, tag="b")
bl3  = Block( [0 0.2 0;0.1 0.3 0.1], nx=1, ny=1, nz=1, cellshape=HEX8, tag="b")
msh = Mesh(bl1,bl2,bl3)
generate_joints_by_tag!(msh, tag="joints")

# finite element analysis

E = 27.e6

mats = [
    "a"      => ElasticSolid(E=E, nu=0.2),
    "b"      => ElasticSolid(E=E, nu=0.2),
    "joints" => ElasticJoint(E=E, nu=0.2, zeta=5),
]


model = Model(msh, mats)


# Boundary conditions
bcs = [
       :(y==0)   => SurfaceBC(ux=0, uy=0, uz=0),
       :(y==0.3) => SurfaceBC(ty=100)
      ]

@test solve!(model, autoinc=true, nincs=20, maxits=3, tol=0.01, quiet=true, scheme=:FE, nouts=10).success


