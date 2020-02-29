using Amaru
using Test

bl = Block( [0 0 0; 1 1 1], nx=4, ny=4, nz=4, cellshape=HEX8, tag="solids")
mesh = Mesh(bl, verbose=false)

mats = [ "solids" => ElasticSolid(E=100.0, nu=0.2) ]

dom = Domain(mesh, mats)

# Loggers
node_log1 = NodeLogger("node_log1.dat")
node_log2 = NodeLogger()
nodes_log1 = NodeGroupLogger("nodes_log1.dat")

faces_log1 = FaceLogger("faces_log1.dat")
faces_log2 = FaceLogger()
edges_log1 = EdgeLogger()
edges_log2 = EdgeLogger()

ip_log1 = IpLogger()
ips_log1 = IpGroupLogger()

loggers = [
           :(x==1 && y==1 && z==1) => node_log1,
           :(x==1 && y==0 ) => node_log2,
           :(x==1 && y==1 ) => nodes_log1,

           :(z==1) => faces_log1,
           :(z==1) => faces_log2,
           :(z==1) => edges_log1,
           :(z==1) => edges_log2,
           #"barra1" => GroupElemLogger(), # TODO

           :(x>0.5 && y>0.5 && z>0.5) => ip_log1,
           :(x>0.5 && y>0.5) => ips_log1,
          ]

setloggers!(dom, loggers)

bcs = [
       :(z==0)         => NodeBC(ux=0, uy=0, uz=0 ),
       :(x==0 || x==1) => NodeBC(ux=0),
       :(z==1)         => FaceBC(tz=-10.0),
      ]

@test solve!(dom, bcs, nincs=4, verbose=false)
println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].vals)
