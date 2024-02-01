using Amaru
using Test

# Mesh:

bls = [
       Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1, tag="solids"),
       BlockInset( [0.5 3 0.2; 0.5 6.0 0.2], curvetype="polyline", tag="bars", jointtag="joints"),
       BlockInset( [0.5 3 0.8; 0.5 6.0 0.8], curvetype="polyline", tag="bars", jointtag="joints", tipjoint=:both, tipjointtag="tips"),
      ]

msh = Mesh(bls)
save(msh, "mesh.vtk")


# Finite elements:
mats = [
        "solids" => MechSolid => LinearElastic => (E=24e3, nu=0.2),
        "bars"   => MechBar => LinearElastic => (E=200e6, A=0.00011),
        "joints" => MechRSJoint => CebRSJoint => (taumax=12, taures=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta=0.5, ks=(12/0.001)*5, kn=5000, p=0.25),
        "tips"   => MechTipJoint => TipJoint => (k=1e8, fixed=false)
       ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)

bcs = [
        :(y==0 && z==0) => NodeBC(ux=0, uy=0, uz=0),
        :(y==6 && z==0) => NodeBC(uz=0),
        :(y==3 && z==0) => EdgeBC(qz=-1.0),
      ]

tol = 0.01
scheme = :Ralston
nincs=20
maxits=3

addstage!(model, bcs, nincs=nincs)

@test solve!(model, autoinc=true, scheme=scheme, tol=tol, maxits=maxits).success
