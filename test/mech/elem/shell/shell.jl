using Amaru

L  = 6
R  = 3
nx = 20
n  = 20

bl1 = Block( [-R 0 0; -R L 0], nx=nx, cellshape=LIN3, tag="shell")
msh = Mesh(bl1, quiet=true, reorder=false)
msh = revolve(msh, angle=180, n=n, base=[0,0,0], axis=[0,1,0])

# Finite element model
mats = [ "shell" => MechShell => LinearElastic => (E=3e4, nu=0.3, thickness=0.03) ]

ana = MechAnalysis(stressmodel="3d")

model = FEModel(msh, mats, ana)

addlogger!(model, :(y==$L/2 && z==$R) => NodeGroupLogger())
addmonitor!(model, :(y==$L/2 && z==$R) => NodeMonitor(:(uz)))

bcs =
    [
    :(y==0) => NodeBC(ux=0, uz=0, rz=0),
    :(y==$L) => NodeBC(ux=0, uz=0, rz=0),
    :(z==0.0) => NodeBC(uz=0, ry=0),
    :(y==$L/2 && z==$R) => NodeBC(fz=-0.01),
    ]
addstage!(model, bcs, nincs=1)

solve!(model, tol=0.1)