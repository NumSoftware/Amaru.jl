# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    $(TYPEDSIGNATURES)

Gets a 3D `Block` by extruding a 2D `block` in the direction given by `axis` and a distance `length`.
It also sets the number of divisions `n` in the extruded direction.

# Examples

```julia
julia> block2d = Block([ 0 0; 1 1 ], nx=3, ny=4, cellshape=QUAD8);
julia> block3d = extrude(block2d, axis=[0,0,1], length=1, n=5)
Block
  nodes: 8-element Vector{Node}:
    1: Node  id=1
    2: Node  id=2
    3: Node  id=3
    4: Node  id=4
    5: Node  id=5
    6: Node  id=6
    7: Node  id=7
    8: Node  id=8
  shape: CellShape  name="HEX8"
  cellshape: CellShape  name="HEX20"
  nx: 3
  ny: 4
  nz: 5
  rx: 1.0
  ry: 1.0
  rz: 1.0
  tag: ""
```
"""
function extrude(block::Block; axis=[0,0,1], length::Number=1.0, n::Int=1)::Block

    if axis==:x
        axis = Vec3(1,0,0)
    elseif axis==:y
        axis = Vec3(0,1,0)
    elseif axis==:z
        axis = Vec3(0,0,1)
    end

    axis = Vec3(normalize(axis))
    Δl = length
    l = 0.0

    if block.shape==QUAD4
        newshape = HEX8
        ls = [l, l+Δl]
        nidx = [1:4;1:4]
        lidx = [1,1,1,1,2,2,2,2]
    elseif block.shape==QUAD8
        newshape = HEX20
        ls = [l, l+Δl/2, l+Δl]
        nidx = [1:4;1:4;5:8;5:8;1:4]
        lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2]
    else
        error("extrude: Block shape $(block.shape.name) is not supported")
    end

    nodes = Node[]
    for (node,li) in zip(block.nodes[nidx], ls[lidx])
        coord = node.coord + li*axis
        push!(nodes, Node(coord))
    end

    return Block(getcoords(nodes), nx=block.nx, ny=block.ny, nz=n, cellshape=newshape)

end

function extrude(blocks::Array; axis=[0,0,1], length=1.0::Number, n=1::Int)
    return extrude.(blocks, axis=axis, length=length, n=n)
end


# Generates a new mesh obtained by extrusion of a 2D mesh
"""
    extrude(mesh, [axis=[0,0,1],] [length=1.0,] [n=1,] [verbosity=1,] [genedges=false])

Generates a 3D mesh by extruding a planar `mesh` using 
a direction `axis`, a `length` and a number of divisions `n`.
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
