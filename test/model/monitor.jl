using Amaru
using Test

bl = Block( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, cellshape=HEX8, tag="solids")
mesh = Mesh(bl)

mats = [ "solids" => MechSolid => LinearElastic => (E=100.0, nu=0.2) ]

ana = MechAnalysis()
model = FEModel(mesh, mats, ana)

# Monitors
monitors = [
    and(x==1, y==1, z==1) => NodeMonitor(:uy, "node.table")
    and(x==1, y==0) => NodeMonitor(:uy)
    y==1 => NodeSumMonitor(:uy, "nodessum.book")

    and(x>0.5, y>0.5, z>0.5) => IpMonitor(:syy, "ip.table")
    and(x>0.5, y>0.5) => IpGroupMonitor(:syy, "ips.book")
    [0.5,0.5,0.0] => IpMonitor(:syy)
]

setloggers!(model, loggers)

bcs = [
       :(z==0) => NodeBC(ux=0, uy=0, uz=0 ),
       :(z==1) => SurfaceBC(tz=-10.0),
      ]
addstage!(model, bcs, nincs=4, nouts=4)

@test solve!(model).success