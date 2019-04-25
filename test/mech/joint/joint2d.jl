using Amaru
using Test

# mesh generation

bl  = Block2D( [0 0; 0.2 0.1], nx=2, ny=6, cellshape=QUAD4, tag="solids")
msh = Mesh(bl, verbose=true)
generate_joints!(msh, tag="joints")
#iptag!(msh.cells[:joints], "jips")

# finite element analysis

E = 27.e6

mats = [
    "solids" => ElasticSolid(E=E, nu=0.2),
    "joints" => MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk" ),
    #"joints" => ElasticJoint(E=E, nu=0.2, alpha=5)
]

dom = Domain(msh, mats, modeltype=:plane_stress, thickness=1.0)

# Loggers
tag!(dom.elems["joints"].ips, "jips")
log1 = IpLogger()
loggers = [
           "jips" => log1
          ]
setloggers!(dom, loggers)

# Boundary conditions
bcs = [
       :(x==0)   => FaceBC(ux=0, uy=0 ),
       :(x==0.2) => FaceBC(ux=2.0*1.7e-4),
      ]

@test solve!(dom, bcs, autoinc=true, nincs=20, maxits=3, tol=0.01, verbose=true, scheme=:ME, nouts=10)

if Amaru.debug.makeplots
    using PyPlot
    table = log1.table
    plot(table[:upa], table[:s1])
end
