using Amaru
using Test

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

# mesh
bl = Block( [0 0;L h], nx=50, ny=2, cellshape=QUAD8)
msh= Mesh(bl)

# fem domain
mat = [ :bulks => MechSolid => VonMises => (E=E, nu=nu, fy=fy, H=H) ]

ana = MechAnalysis(stressmodel=:planestress, thickness=th)
model = FEModel(msh, mat, ana)

log = NodeLogger()
addlogger!(model, :(y==$h/2 && x==1) => log)
addmonitor!(model, :(y==$h/2 && x==1) => NodeMonitor(:fy))

# boundary conditions
bcs = [
    :(x==0) => NodeBC(ux=0),
    :(x==0 && y==$h/2) => NodeBC(uy=0),
    :(x==1.0 && y==$h/2) => NodeBC(uy = -0.08),
]

addstage!(model, bcs, nincs=30, nouts=1)

solve!(model, autoinc=true)

println(@test log.table.fy[end]≈-30 atol=0.4)

if makeplots
    using PyPlot
    tab = log.table
    plot( -tab[:uy], -tab[:fy], "-o")
end

