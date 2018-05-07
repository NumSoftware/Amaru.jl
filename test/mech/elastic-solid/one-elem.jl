# From the book: 
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru

# Mesh generation

block = Block3D( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, shape=HEX8) 
mesh = Mesh(block, verbose=true)

# Domain definition

materials = [
    MaterialBind(0, ElasticSolid(E=1.0, nu=0.3) ),
]

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
    println(D)

    println("Stress:")
    S = elems_ip_vals(dom.elems[1])[[:sxx, :syy, :szz, :syz, :sxz, :sxy]]
    println(S)

    println("Support reactions:")
    F = nodes_dof_vals(dom.nodes[:(z==0)])[[:fx, :fy, :fz]]
    println(F)
end


