# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
`extrude(block, [axis=[0,0,1],] [len=1.0,] [n=1])`

Extrudes a 2D `block` generating a 3D block based on a direction `axis`, a lenght `len` and a number `n` of divisions.
"""
function extrude(block::Block2D; axis=[0,0,1], len::Number=1.0, n::Int=1)::Block3D

    if axis==:x
        axis = Vec3(1,0,0)
    elseif axis==:y
        axis = Vec3(0,1,0)
    elseif axis==:z
        axis = Vec3(0,0,1)
    end

    V = axis/norm(axis)
    δ = len/n
    coords = get_coords(block.nodes)

    # check numbering order
    v1 = vec(coords[1,:])
    v2 = vec(coords[2,:])
    v3 = vec(coords[3,:])
    normal_order = dot(cross(v2-v1, v3-v2), axis) > 0.0

    # generate extra points
    npoints = size(coords,1)
    topcoords = coords .+ V'*len

    local newcoords::Array{Float64,2}

    if npoints==4
        if normal_order
            newcoords = vcat(coords, topcoords)
        else
            newcoords = vcat(topcoords, coords)
        end
    end

    if npoints==8
        midcoords = coords[1:4,:] .+ V'*len*0.5
        if normal_order
            newcoords = vcat(coords[1:4,:], topcoords[1:4,:], coords[5:8,:], topcoords[5:8,:], midcoords)
        else
            newcoords = vcat(coords[5:8,:], topcoords[5:8,:], coords[1:4,:], topcoords[1:4,:], midcoords)
        end
    end

    #shape = HEX8
    #if block.shape==TRI3
        #shape = WED6
    #elseif block.shape==TRI6
        #shape = WED15
    #elseif block.shape==QUAD8
        #shape = HEX20
    #end

    shape = block.cellshape==QUAD8 ? HEX20 : HEX8

    return Block3D( newcoords, nx=block.nx, ny=block.ny, nz=n, cellshape=shape)

end

function extrude(blocks::Array; axis=[0,0,1], len=1.0::Number, n=1::Int)::Array{Block3D,1}
    if axis==:x
        axis = Vec3(1,0,0)
    elseif axis==:y
        axis = Vec3(0,1,0)
    elseif axis==:z
        axis = Vec3(0,0,1)
    end

    blocks3D = []

    for bl in blocks
        bl3D = extrude(bl, axis=axis, len=len, n=n)
        push!(blocks3D, bl3D)
    end

    return blocks3D
end


# Generates a new mesh obtained by extrusion of a 2D mesh
"""
    extrude(mesh, [axis=[0,0,1],] [len=1.0,] [n=1,] [verbose=true,] [genedges=false])

Generates a 3D mesh by extruding a planar `mesh` using 
a direction `axis`, a `lenght` and a number of divisions `n`.
"""
function extrude(mesh::Mesh; axis=[0,0,1], length::Number=1.0, n::Int=1, verbose::Bool=true)
    if axis==:x
        axis = Vec3(1,0,0)
    elseif axis==:y
        axis = Vec3(0,1,0)
    elseif axis==:z
        axis = Vec3(0,0,1)
    end

    axis = Vec3(normalize(axis))

    verbose && printstyled("Mesh extrude:\n", bold=true, color=:cyan)
    Δl = length/n

    # Generate new cells
    cells = Cell[]

    for l in range(0, step=Δl, length=n)
        for cell in mesh.elems

            if cell.shape==LIN2
                newshape = QUAD4
                ls = [l, l+Δl]
                nidx = [1,2,2,1]
                lidx = [1,1,2,2]
            elseif cell.shape==LIN3
                newshape = QUAD8
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1,2,2,1,3,2,3,1]
                lidx = [1,1,3,3,1,2,3,2]
            elseif cell.shape==TRI3
                newshape = WED6
                ls = [l, l+Δl]
                nidx = [1,2,3,1,2,3]
                lidx = [1,1,1,2,2,2]
            elseif cell.shape==QUAD4
                newshape = HEX8
                ls = [l, l+Δl]
                nidx = [1:4;1:4]
                lidx = [1,1,1,1,2,2,2,2]
            elseif cell.shape==TRI6
                newshape = WED15
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:4;1:4]
                lidx = [1,1,1,1,2,2,2,2]
            elseif cell.shape==QUAD8
                newshape = HEX20
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:4;1:4;5:8;5:8;1:4]
                lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2]
            elseif cell.shape==QUAD9
                newshape = HEX27
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:4;1:4;5:8;5:8;1:4; 8;6;5;7;9;9;9]
                lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2, 2,2,2,2,1,3,2]
            else
                error("extrude: Cell shape $(cell.shape.name) is not supported")
            end

            nodes = Node[]
            for (node,li) in zip(cell.nodes[nidx], ls[lidx])
                coord = node.coord + li*axis
                push!(nodes, Node(coord))
            end

            newcell = Cell(newshape, nodes, tag=cell.tag)
            isinverted(newcell) && flip!(newcell)
            push!(cells, newcell)
        end
    end

    # Merge points
    node_dict = Dict{UInt64,Node}( hash(n) => n for c in cells for n in c.nodes )

    for cell in cells
        cell.nodes = Node[ node_dict[hash(n)] for n in cell.nodes ]
    end
    nodes = collect(values(node_dict))

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    fixup!(newmesh, reorder=true)

    if verbose
        @printf "  %5d points obtained\n" Base.length(newmesh.nodes)
        @printf "  %5d cells obtained\n" Base.length(newmesh.elems)
        @printf "  %5d faces obtained\n" Base.length(newmesh.faces)
        @printf "  %5d surface edges obtained\n" Base.length(newmesh.edges)
    end

    return newmesh

end
