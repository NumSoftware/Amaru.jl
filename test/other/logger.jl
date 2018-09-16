using Amaru
using Test

bl = Block3D( [0 0 0; 1 1 1], nx=4, ny=4, nz=4, cellshape=HEX8)
mesh = Mesh(bl, verbose=false)

mat = MaterialBind(:solids, ElasticSolid(E=100.0, nu=0.2) )

# Loggers
node_mon1 = NodeLogger(:(x==1 && y==1 && z==1) )
node_mon2 = NodeLogger(:(x==1 && y==0 ) )
node_mon3 = NodeLogger(1)

nodes_mon1 = GroupLogger(:node, :(x==1 && y==1 ) )

faces_mon1 = FaceLogger(:(z==1) )
faces_mon2 = FaceLogger(:(z==1) )
edges_mon1 = EdgeLogger(:(z==1) )
edges_mon2 = EdgeLogger(:(z==1) )
#barra1_mon = GroupLogger(:elem, "barra1") TODO

ip_mon1   = IpLogger(1)
ip_mon2   = IpLogger(:(x>0.5 && y>0.5 && z>0.5) )

ips_mon1   = GroupLogger(:ip, :(x>0.5 && y>0.5) )

logs = [ node_mon1, node_mon2, node_mon3, nodes_mon1, faces_mon1, faces_mon2,
        edges_mon1, edges_mon2, ip_mon1, ip_mon2, ips_mon1 ]

dom = Domain(mesh, mat, logs)

bcs = [
       NodeBC(:(z==0), :(ux=0, uy=0, uz=0 )),
       NodeBC(:(x==0 || x==1), :(ux=0)),
       FaceBC(:(z==1), :(tz=-10.0)),
      ]

@test solve!(dom, bcs, nincs=1, verbose=true)
println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].vals)
