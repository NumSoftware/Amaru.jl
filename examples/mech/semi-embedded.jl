using Amaru

# Mesh generation
bl  = Block( [0 0 0; 0.5 6.0 0.5], nx=1, ny=10, nz=3, tag="solids")
bl1 = BlockInset( [0.05 0.05 0.05; 0.05 5.95 0.05], curvetype="polyline", tag="bars", jointtag="interface")
bl2 = copy(bl1, dx=0.4)

msh = Mesh(bl, bl1, bl2)


h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0.0
nu = 0.3



# FEM analysis
mats = [
        "solids" => MechSolid => LinearElastic => (E=24e3, nu=0.25),
        "bars"  => MechTruss => VonMises => (E=200e6, A=0.0001, fy=500e3),
        "interface" => MechLSInterface => LinearLSInterface => (ks=1e9, kn=1e9, p=0.02),
]

ctx = MechContext()
model = FEModel(msh, mats, ctx)
ana = MechAnalysis(model)

bcs = [
       (y==0, z==0) => NodeBC(ux=0, uy=0, uz=0),
       (y==6, z==0) => NodeBC(ux=0, uz=0),
       z==0.5 => SurfaceBC(tz=-0.001),
]
addstage!(ana, bcs, nincs=20)

solve!(ana, autoinc=true)

# mplot(model, "beam.pdf", 
#     field="sa", 
#     fieldmult=1e-3,
#     axis=false,
#     opacity=0.1,
#     dist=6,
#     colorbarlabel="axial stress in bars"
# )
