@show ARGS

using Amaru
using Base.Test

# mesh generation

bl  = Block2D( [0 0; 0.2 0.1], nx=2, ny=4, shape=QUAD4)
msh = Mesh(bl, verbose=true)
generate_joints!(msh)
iptag!(msh.cells[:joints], 100)

# finite element analysis

E = 27.e6

mats = [
    MaterialBind(:solids, ElasticSolid(E=E, nu=0.2)),
    MaterialBind(:joints, MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5, model="hordijk" ) ),
    #MaterialBind(:joints, ElasticJoint(E=E, nu=0.2, alpha=5), iptag="jnt_ip" ),
]

# Loggers
logger = Logger(:ip, 100)

dom = Domain(msh, mats, logger, model_type=:plane_stress, thickness=1.0)

# Boundary conditions
bcs = [
       BC(:face, :(x==0), :(ux=0, uy=0 )),
       BC(:face, :(x==0.2), :(ux=2.0*1.7e-4)),
      ]

@test solve!(dom, bcs, autoinc=true, nincs=20, maxits=3, tol=0.01, verbose=true, scheme=:ME, nouts=10)

save(dom, "dom1.vtk")
