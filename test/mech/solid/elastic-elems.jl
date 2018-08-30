using Amaru
using Test

dis = [ -0.012, -0.095 ]

for shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9) 
    printstyled(shape.name, color=:cyan); println()
    bl = Block2D( [0 0; 1 1], nx=5, ny=5, shape=shape)
    mesh = Mesh(bl, verbose=false)
    tag!(mesh.faces[:(y==0)], "bottom") # bottom face
    tag!(mesh.faces[:(y==1)], "top") # top face

    materials = [
        MaterialBind(0, ElasticSolid(E=100.0, nu=0.2) ),
    ]

    dom = Domain(mesh, materials, verbose=false)


    bcs = [
        FaceBC("bottom", :(ux=0, uy=0) ),
        FaceBC("top", :(ty=-10.) ),
    ]

    solve!(dom, bcs, nincs=1)

    top_node = dom.nodes[:(y==1)][1]
    ux = top_node.dofdict[:ux].vals[:ux]
    uy = top_node.dofdict[:uy].vals[:uy]
    @test [ux, uy] ≈ dis atol=4e-2

    println( nodes_dof_vals(dom.nodes[:(y==1)][1]) )

end

for shape in (TET4, TET10, HEX8, HEX20)
    printstyled(shape.name, color=:cyan); println()
    bl = Block3D( [0 0 0; 1 1 1], nx=5, ny=5, nz=5, shape=shape)
    mesh = Mesh(bl, verbose=false)
    tag!(mesh.faces[:(z==0)], "bottom") # bottom face
    tag!(mesh.faces[:(z==1)], "top") # top face
    tag!(mesh.faces[:(x==0 || x==1)], "sides") # lateral face

    materials = [
        MaterialBind(0, ElasticSolid(E=100.0, nu=0.2) ),
    ]

    dom = Domain(mesh, materials, verbose=false)

    bcs = [
        FaceBC("bottom", :(ux=0, uy=0, uz=0) ),
        FaceBC("sides", :(ux=0) ),
        FaceBC("top", :(tz=-10.) ),
    ]

    solve!(dom, bcs, nincs=1)

    top_node = dom.nodes[:(z==1)][1]
    uy = top_node.dofdict[:uy].vals[:uy]
    uz = top_node.dofdict[:uz].vals[:uz]
    @test [uy, uz] ≈ dis atol=1e-2
end
