using Amaru
using Test

# Mesh generation
# ===============

bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=3, tag="solids")
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=true, tag="bars")
bl2 = copy(bl1)
move!(bl2, dx=0.6)
bls = [ bl, bl1, bl2 ]

msh = Mesh(bls, verbose=true)

# FEM analysis
# ============

mats = [
        "solids" => ElasticSolid(E=1.e4, nu=0.25),
        "bars"  => PPRod(E=1.e8, A=0.005, sig_y=500e3),
       ]

dom = Domain(msh, mats)

#setstate!(dom.elems["solids"], ElasticSolidState())
#states = [
          #"solids" => ElasticSolidState()
         #]

#setstates!(dom, states)

bcs = [
       :(y==0 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(y==6 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(z==1) => FaceBC(tz=-1000),
      ]

solve!(dom, bcs, nincs=20, verbose=true)
save(dom, "domain.vtk")

println("Available data fields: \n", datafields(dom))
mplot(dom, "beam.pdf", field="sa", fieldscale=1e-3, axis=false, alpha=0.1, dist=6, colorbarlabel="axial stress in bars")
