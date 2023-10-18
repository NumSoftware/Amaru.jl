using Amaru
using Test

h  = 0.1
th = 0.05
L  = 1.0
E  = 210e6 # kPa
fy = 240e3 # kPa
H  = 0.0
# H  = 1e-3
nu = 0.3

# mesh
bl = Block( [0 0 0; L 0 0], nx=50, cellshape=LIN3)
msh = Mesh(bl, ndim=3)

# fem domain
mat = [ :lines => MechBeam => VonMises => (E=E, nu=nu, fy=fy, H=H, thy=th, thz=h) ]

ana = MechAnalysis()
model = FEModel(msh, mat, ana)

log = NodeLogger()
addlogger!(model, :(x==$L) => log)
addmonitor!(model, :(x==$L) => NodeMonitor(:fz))

# boundary conditions
bcs = [
    :(x==0) => NodeBC(ux=0, uy=0, uz=0, rx=0, ry=0, rz=0),
    :(x==$L) => NodeBC(uz = -0.08),
]

addstage!(model, bcs, nincs=30, nouts=1)

solve!(model, autoinc=true)

println(@test log.table.fz[end]≈-30 atol=5.0)

if makeplots
    using PyPlot
    tab = log.table
    plot( -tab[:uz], -tab[:fz], "-o")
    show()
end

