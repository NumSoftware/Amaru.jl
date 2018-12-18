using Amaru
using Test

# mesh generation

bl1  = Block3D( [0 0 0; 0.1 0.1 0.1], nx=1, ny=1, cellshape=HEX8, tag="a")
bl2  = Block3D( [0 0.1 0;0.1 0.2 0.1], nx=1, ny=1, cellshape=HEX8, tag="b")
bl3  = Block3D( [0 0.2 0;0.1 0.3 0.1], nx=1, ny=1, cellshape=HEX8, tag="b")
msh = Mesh(bl1,bl2,bl3, verbose=true)
generate_joints_tag!(msh)

# finite element analysis

E = 27.e6

mats = [
    MaterialBind(:solids, ElasticSolid(E=E, nu=0.2)),
    MaterialBind(:joints, ElasticJoint(E=E, nu=0.2, alpha=5)),
]


dom = Domain(msh, mats)


# Boundary conditions
bcs = [
       FaceBC(:(y==0), :(ux=0, uy=0, uz=0)),
       FaceBC(:(y==0.3), :(ty=100)),
      ]

try
    @test solve!(dom, bcs, autoinc=true, nincs=20, maxits=3, tol=0.01, verbose=true, scheme=:FE, nouts=10)
catch
end


