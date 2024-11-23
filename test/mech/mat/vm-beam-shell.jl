using Amaru
using Test

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0.0
nu = 0.3

# mesh
bl = Block( [0 0;L h], nx=50, ny=2, cellshape=QUAD8)
msh= Mesh(bl, ndim=3)

save(msh, "mesh.vtk")

# fem domain
mat = [ :bulks => MechShell => VonMises => (E=E, nu=nu, fy=fy, H=H, thickness=th) ]

ctx = MechContext()
model = FEModel(msh, mat, ctx)
ana = MechAnalysis(model)

log = NodeLogger()
addlogger!(ana, :(y==$h/2 && x==1) => log)
addmonitor!(ana, :(y==$h/2 && x==1) => NodeMonitor(:fy))

# boundary conditions
bcs = [
    x==0 => NodeBC(ux=0, rx=0, ry=0, rz=0),
    (x==0, y==h/2) => NodeBC(uy=0),
    (x==1.0, y==h/2) => NodeBC(uy = -0.08),
]

addstage!(ana, bcs, nincs=30, nouts=10)

solve!(ana, autoinc=true)

println(@test log.table.fy[end]≈-30 atol=0.4)

if makeplots
    using PyPlot
    tab = log.table
    plot( -tab[:uy], -tab[:fy], "-o")
end

