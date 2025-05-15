using Amaru
using Test

# Mesh generation
b = 0.5
h = 0.6
ℓ = 4.0
bl  = Block( [0 0 0; b ℓ h], nx=1, ny=25, nz=3, tag="solid")
# bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=true, tag="embedded")
# bl2 = copy(bl1)
# move!(bl2, dx=0.6)
# bls = [ bl, bl1, bl2 ]

# mesh = Mesh(bls)
geo = GeoModel()
addblock!(geo, bl)

p1 = addpoint!(geo, 0.2, 0.1, 0.1)
p2 = addpoint!(geo, 0.2, 3.6, 0.1)
px = addpoint!(geo, 0.2, 3.6, 0.4)
p3 = addpoint!(geo, 0.2, 3.9, 0.4)
# addsubpath!(geo, :M, p1, :L, p2, embedded=true, tag="embedded")
addsubpath!(geo, :M, p1, :L, p2, :A, px, p3, embedded=true, tag="embedded")

mesh = Mesh(geo)


# FEM analysis
mats = [
        "solid" => MechSolid => LinearElastic => (E=1.e4, nu=0.25),
        "embedded" => MechBar => VonMises => (E=1.e8, fy=500e3, A=0.005),
        # "embedded" => MechBar => LinearElastic => (E=1.e8, A=0.005),
       ]

ctx = MechContext()
model = FEModel(mesh, mats, ctx)
ana = MechAnalysis(model)

bcs = [
       (y==0, z==0) => NodeBC(ux=0, uy=0, uz=0),
       (y==ℓ, z==0) => NodeBC(ux=0, uy=0, uz=0),
       z==h => SurfaceBC(tz=-1000),
      ]
addstage!(ana, bcs, nincs=1)

solve!(ana, autoinc=true)

save(model, "embedded.vtu")

mplot = MeshPlot(model,
    view_mode  = :wireframe,
    elevation = 3,
    azimut    = 5,
)
save(mplot, "embedded.pdf")