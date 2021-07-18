using Amaru
using Test

dis = [ -0.012, -0.095 ]

for shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9)
    printstyled(shape.name, color=:cyan); println()
    bl = Block( [0 0; 1 1], nx=2, ny=2, cellshape=shape, tag="solids")
    mesh = Mesh(bl, printlog=false)
    tag!(mesh.faces[:(y==0)], "bottom") # bottom face
    tag!(mesh.faces[:(y==1)], "top") # top face

    materials = [
        "solids" => ElasticSolid(E=100.0, nu=0.2)
    ]

    dom = Domain(mesh, materials, printlog=false)


    bcs = [
        "bottom" => FaceBC(ux=0, uy=0),
        "top"    => FaceBC(ty=-10.)
    ]

    solve!(dom, bcs, nincs=1, printlog=false).success

    top_node = dom.nodes[:(y==1)][1]
    ux = top_node.dofs[:ux].vals[:ux]
    uy = top_node.dofs[:uy].vals[:uy]
    @test [ux, uy] ≈ dis atol=4e-2

    println( get_data(dom.nodes[:(y==1)][1]) )

end

#for shape in (TET4, TET10, HEX8, HEX20, HEX27)
for shape in (TET4, TET10, HEX8, HEX20, HEX27)
    printstyled(shape.name, color=:cyan); println()
    bl = Block( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, cellshape=shape, tag="solids")
    mesh = Mesh(bl, printlog=false)
    tag!(mesh.faces[:(z==0)], "bottom") # bottom face
    tag!(mesh.faces[:(z==1)], "top") # top face
    tag!(mesh.faces[:(x==0 || x==1)], "sides") # lateral face

    materials = [
        "solids" => ElasticSolid(E=100.0, nu=0.2)
    ]

    dom = Domain(mesh, materials, printlog=false)

    bcs = [
        "bottom" => FaceBC(ux=0, uy=0, uz=0),
        "sides"  => FaceBC(ux=0),
        "top"    => FaceBC(tz=-10.)
    ]

    solve!(dom, bcs, nincs=1, printlog=false).success

    top_node = dom.nodes[:(z==1)][1]
    uy = top_node.dofs[:uy].vals[:uy]
    uz = top_node.dofs[:uz].vals[:uz]

    println( get_data(dom.nodes[:(z==1)][1]) )

    @test [uy, uz] ≈ dis atol=1e-2
end
