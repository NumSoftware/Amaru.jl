using Amaru
using Test


th = 0.05
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0
nu = 0.3

# mesh
bl = Block( [0 0 -0.05; 0.05 1.0 0.05], nx=1, ny=50, nz=2, cellshape=HEX20)
msh= Mesh(bl)

# fem domain
mat = [ :bulks => MechSolid => VonMises => (E=E, nu=nu, fy=fy, H=H) ]

ana = MechAnalysis(stressmodel="3d")
model = FEModel(msh, mat, ana)

log = NodeLogger()
addlogger!(model, :(x==$th/2 && y==1 && z==0) => log)
addmonitor!(model, :((x==$th/2 && y==1 && z==0)) => NodeMonitor(:fz))

# boundary conditions
bcs = [
    :(y==0) => NodeBC(uy=0),
    :(y==0 && z==0) => NodeBC(uz=0),
    :(x==$th/2 && y==0 && z==0) => NodeBC(ux=0),
    :(x==$th/2  && y==1 && z==0) => NodeBC(uz=-0.08),
]

addstage!(model, bcs, nincs=20, nouts=1)
solve!(model, autoinc=true)

println(@test log.table.fz[end]≈-30 atol=0.7)

if makeplots
    using PyPlot
    tab = log.table
    plot( -tab[:uz], -tab[:fz], "-o")
end

