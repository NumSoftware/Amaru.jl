using Amaru
using Test

# Mesh generation
# ===============

bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=3, tag="solid")
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=false, tag="bars", jointtag="joint")
bl2 = copy(bl1)
move!(bl2, dx=0.6)
bls = [ bl, bl1, bl2 ]

msh = Mesh(bls, verbosity=1)

# FEM analysis
# ============

mats = [
        "solid" => ElasticSolid(E=1.e4, nu=0.25),
        "bars"  => PPRod(E=1.e8, A=0.005, sig_y=500e3),
        "joint" => ElasticJoint1D(ks=1e8, kn=1e8, A=0.01),
       ]


dom = Domain(msh, mats)

bcs = [
       :(y==0 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(y==6 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(z==1) => FaceBC(tz=-1000),
      ]

solve!(dom, bcs, nincs=20, verbosity=1)
save(dom, "domain.vtk")

mplot(dom, "beam.pdf", field="sa", fieldmult=1e-3, axis=false, opacity=0.1, dist=6, colorbarlabel="axial stress in bars")
