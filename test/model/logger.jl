using Amaru
using Test

bl = Block( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, cellshape=HEX8, tag="solids")
mesh = Mesh(bl)

mats = [ "solids" => MechSolid => LinearElastic => (E=100.0, nu=0.2) ]

ana = MechAnalysis()
model = FEModel(mesh, mats, ana)

# Loggers
loggers = [
    and(x==1, y==1, z==1) => NodeLogger("node.table")
    and(x==1, y==0) => NodeLogger()
    and(x==1, y==1) => NodeGroupLogger("nodes.book")

    z==1 => FaceLogger("face.table")
    z==1 => EdgeLogger("edge.table")

    and(x>0.5, y>0.5, z>0.5) => IpLogger("ip.table")
    and(x>0.5, y>0.5) => IpGroupLogger("ips.book")
    [0.5,0.5,0.0] => IpLogger()
    [0.5,0.5,0.0] => PointLogger("point.table")
    [0.5 0.5 0.0; 1 1 0] => SegmentLogger("segment.book")
]

setloggers!(model, loggers)

bcs = [
       :(z==0) => NodeBC(ux=0, uy=0, uz=0 ),
       :(z==1) => SurfaceBC(tz=-10.0),
      ]
addstage!(model, bcs, nincs=4, nouts=4)

@test solve!(model).success