using Amaru
using Base.Test

bl  = Block2D( [0 0; 0.2 0.1], nx=2, ny=1, shape=QUAD4)

# mesh generation
msh = Mesh(bl, verbose=true)
generate_joints!(msh)


E = 27.e6

mats = [
    MaterialBind(:solids, ElasticSolid(E=E, nu=0.2)),
    MaterialBind(:joints, MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5 ), iptag="jnt_ip" ),
    #MaterialBind(:joints, ElasticJoint(E=E, nu=0.2, alpha=5), iptag="jnt_ip" ),
]

# Loggers
logger = Logger(:ip, "jnt_ip")

dom = Domain(msh, mats, logger, model_type=:plane_stress, thickness=1.0)

# Boundary conditions
bcs = [
       BC(:face, :(x==0), :(ux=0, uy=0 )),
       BC(:face, :(x==0.2), :(ux=2.0*1.7e-4)),
      ]

@test solve!(dom, bcs, autoinc=true, nincs=20, maxits=3, tol=0.01, verbose=true, scheme=:ME)

save(dom, "dom1.vtk")
