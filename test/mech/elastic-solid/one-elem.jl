# From the book: 
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru

# Mesh generation

block = Block3D( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, shape=HEX8) 

mesh = Mesh(block, verbose=true)

#tag!(mesh.points[1], 1)       # node 1
#tag!(mesh.faces[:(z==0)], 10) # top face
#tag!(mesh.faces[:(z==1)], 20) # bottom face


# Domain definition

materials = [
    MaterialBind(0, ElasticSolid(E=1.0, nu=0.3) ),
]

#loggers = [
    #Logger(:node, 1)  # by tag
    #Logger(:ip, :(x>0))  # by expression
    #GroupLogger(:nodes, "nodes")  # by tag
    #GroupLogger(:ips, :(x>0))  # by expression
#]


# Load cases

bcs1 = [
    BC(:node, :(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    BC(:node, :(x==1 && y==0 && z==0), :(uy=0) ),
    BC(:node, :(x==0 && y==1 && z==0), :(ux=0) ),
    BC(:node, :(z==0), :(uz=0) ),
    BC(:node, :(z==1), :(fz=1) ),
]

bcs2 = [
    BC(:node, :(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    BC(:node, :(x==1 && y==0 && z==0), :(uy=0) ),
    BC(:node, :(x==0 && y==1 && z==0), :(ux=0) ),
    BC(:node, :(z==0), :(uz=0) ),
    BC(:edge, :(y==1 && z==1), :(ty=2) ),
]

bcs3 = [
    BC(:node, :(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    BC(:node, :(x==1 && y==0 && z==0), :(uy=0) ),
    BC(:node, :(x==0 && y==1 && z==0), :(ux=0) ),
    BC(:node, :(z==0), :(uz=0) ),
    BC(:face, :(x==1), :(tx=3*z) ),
]

bcs4 = [
    BC(:node, :(x==0 && y==0 && z==0), :(ux=0, uy=0) ),
    BC(:node, :(x==1 && y==0 && z==0), :(uy=0) ),
    BC(:node, :(x==0 && y==1 && z==0), :(ux=0) ),
    BC(:node, :(z==0), :(uz=0) ),
    BC(:element, :all, :(tz=-1) ),
]

for bcs in [bcs1, bcs2, bcs3, bcs4]

    dom = Domain(mesh, materials)
    solve!(dom, bcs, nincs=1, nouts=1, verbose=true)

    println("Displacements:")
    D = nodes_dof_vals(dom.nodes)[[:ux, :uy, :uz]]
    #D = [ node.dofdict[key].vals[key] for node in dom.nodes, key in (:ux, :uy, :uz) ]
    println(D)

    println("Stress:")
    S = elem_ip_vals(dom.elems[1])[[:sxx, :syy, :szz, :syz, :sxz, :sxy]]
    println(S)

    println("Support reactions:")
    F = nodes_dof_vals(dom.nodes[:(z==0)])[[:fx, :fy, :fz]]
    println(F)
end


