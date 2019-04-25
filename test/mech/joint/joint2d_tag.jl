using Amaru
using Test

# mesh generation

bl1  = Block2D( [0 0; 0.1 0.1], nx=1, ny=1, cellshape=QUAD4, tag="a")
bl2  = Block2D( [0.1 0; 0.2 0.1], nx=1, ny=1, cellshape=QUAD4, tag="b")
bl3  = Block2D( [0.2 0; 0.3 0.1], nx=1, ny=1, cellshape=QUAD4, tag="b")
msh = Mesh(bl1,bl2,bl3, verbose=true)
generate_joints_by_tag!(msh, tag="joints")

# finite element analysis

E = 27.e6

mats = [
    "a"      => ElasticSolid(E=E, nu=0.2),
    "b"      => ElasticSolid(E=E, nu=0.2),
    "joints" => ElasticJoint(E=E, nu=0.2, alpha=5),
]


dom = Domain(msh, mats, model_type=:plane_stress, thickness=1.0)

# Boundary conditions
bcs = [
       :(x==0)   => FaceBC(ux=0, uy=0 )
       :(x==0.3) => FaceBC(ux=2.0*1.7e-4)
      ]

@test solve!(dom, bcs, autoinc=true, nincs=20, maxits=3, tol=0.01, verbose=true, scheme=:ME, nouts=10)


