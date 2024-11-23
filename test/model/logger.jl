using Amaru
using Test

bl = Block( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, cellshape=HEX8, tag="solids")
mesh = Mesh(bl)

mats = [ "solids" => MechSolid => LinearElastic => (E=100.0, nu=0.2) ]

ctx = MechContext()
model = FEModel(mesh, mats, ctx)
ana = MechAnalysis(model)

# Loggers
loggers = [
    (x==1, y==1, z==1) => NodeLogger("node.table")
    (x==1, y==0) => NodeLogger()
    (x==1, y==1) => NodeGroupLogger("nodes.book")

    z==1 => FaceLogger("face.table")
    z==1 => EdgeLogger("edge.table")

    (x>0.5, y>0.5, z>0.5) => IpLogger("ip.table")
    (x>0.5, y>0.5) => IpGroupLogger("ips.book")
    [0.5,0.5,0.0] => IpLogger()
    [0.5,0.5,0.0] => PointLogger("point.table")
    [0.5 0.5 0.0; 1 1 0] => SegmentLogger("segment.book")
]

setloggers!(ana, loggers)

bcs = [
       z==0 => NodeBC(ux=0, uy=0, uz=0 ),
       z==1 => SurfaceBC(tz=-10.0),
      ]
addstage!(ana, bcs, nincs=4, nouts=4)

@test solve!(ana).success