# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


#Base.copy(p::Node) = Node(p.coord.x, p.coord.y, p.coord.z, p.tag)
#Base.copy(c::Cell)  = Cell(c.shape, c.nodes, tag=c.tag, oelem=c.oelem)
Base.copy(c::Cell)  = Cell(c.shape, c.nodes, tag=c.tag, oelem=c.oelem)

function Base.copy(bls::Array{<:AbstractBlock,1})
    return [ copy(obj) for obj in bls ]
end


"""
`move(block, [dx=0.0,] [dy=0.0,] [dz=0.0])`

Changes the coordinates of a `block`. Also returns a reference.
"""
function move!(bl::AbstractBlock; dx=0.0, dy=0.0, dz=0.0)
    for p in bl.nodes
        p.coord.x += dx
        p.coord.y += dy
        p.coord.z += dz
    end
    return bl
end


"""
`move(blocks, [dx=0.0,] [dy=0.0,] [dz=0.0])`

Changes the coordinates of an array of blocks. Also returns a reference.
"""
function move!(blocks::Array; dx=0.0, dy=0.0, dz=0.0)
    for bl in blocks
        for p in bl.nodes
            p.coord.x += dx
            p.coord.y += dy
            p.coord.z += dz
        end
    end
    return blocks
end


function scale!(bl::AbstractBlock; factor=1.0, base=[0.,0,0], axis=nothing)
    coords =get_coords(bl)
    coords = ( base .+ factor*(coords' .- base) )'

    for (i,p) in enumerate(bl.nodes)
        p.coord.x, p.coord.y, p.coord.z = coords[i,:]
    end

    return bl
end


function scale!(blocks::Array{AbstractBlock,1}; factor=1.0, base=[0.0,0,0])
    for bl in blocks
        scale!(bl, factor=factor, base=base)
    end
    return blocks
end

function mirror(block::AbstractBlock; face=[0.0 0 0; 0 1 0; 0 0 1])
    nr, nc = size(face)
    if nc==2
        face = [ face zeros(nr) ]
    end
    if nr==2
        face = [ face; [face[1,1] face[1,2] 1.0] ]
    end

    p1 = face[1,:]
    p2 = face[2,:]
    p3 = face[3,:]
    normal = cross(p2-p1, p3-p1)
    normal = normal/norm(normal)

    bl = copy(block)
    coords =get_coords(bl.nodes)

    distances    = (coords .- p1')*normal       # d = n^.(xi - xp)
    coords = coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

    # fix coordinates in bl to keep anti-clockwise numbering
    npts = size(coords)[1]
    ndim = typeof(bl)==Block2D ? 2 : 3

    if npts==8 && ndim==2
        idxs = [ 4, 3, 2, 1, 7, 6, 5, 8 ]
    elseif npts==8 && ndim==3
        idxs = [ 5:8; 1:4 ]
    elseif npts==20 && ndim==3
        idxs = [ 5:8; 1:4; 13:16; 9:12; 17:20 ]
    else
        idxs = [ npts:-1:1; ]  # reverse
    end
    coords = coords[idxs,:]
    bl.nodes = [ Node(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]

    return bl
end

function mirror(blocks::Array{<:AbstractBlock,1}; face=[0.0 0 0; 0 1 0; 0 0 1])
    return [ mirror(bl, face=face) for bl in blocks ]
end


function mirror(mesh::Mesh; face=[0.0 0 0; 0 1 0; 0 0 1])
    nr, nc = size(face)
    if nc==2
        face = [ face zeros(nr) ]
    end
    if nr==2
        face = [ face; [face[1,1] face[1,2] 1.0] ]
    end
    p1 = face[1,:]
    p2 = face[2,:]
    p3 = face[3,:]
    normal = cross(p2-p1, p3-p1)
    normal = normal/norm(normal)

    # copy mesh
    newmesh = copy(mesh)

    # mirror
    coords =get_coords(newmesh.nodes)
    distances = (coords .- p1')*normal       # d = n^.(xi - xp)
    coords    = coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

    # updating points
    for (i,p) in enumerate(newmesh.nodes)
        p.coord.x = coords[i,1]
        p.coord.y = coords[i,2]
        p.coord.z = coords[i,3]
    end

    # fix connectivities
    for c in newmesh.elems
        if c.shape==HEX8
            idxs = [ 5:8; 1:4 ]
            c.nodes = c.nodes[idxs]
        end
    end

    fixup!(newmesh)
    return newmesh
end


"""
    array(block; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)

Creates a list with copies of `block` along x, y and z axes.
"""
function array(bl::AbstractBlock; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)
    blocks = [ bl ]
    for k=0:nz-1
        for j=0:ny-1
            for i=0:nx-1
                i==j==k==0 && continue
                cp = copy(bl)
                move!(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                push!(blocks, cp)
            end
        end
    end
    return blocks
end


"""
    rotate!(block, base=[0,0,0], axis=[0,0,1], angle=90.0)`

Rotate `block` according to the provided `base` point, `axis` vector and `angle`.
"""
function rotate!(bl::AbstractBlock; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )
    # see also: https://lucidar.me/en/quaternions/quaternions-rotations/

    length(axis)==2 && ( axis=vcat(axis, 0.0) )
    length(base)==2 && ( base=vcat(base, 0.0) )

    # unit vector
    axis = axis/norm(axis)
    a, b, c = axis
    d = sqrt(b^2+c^2)

    # unit vector for rotation
    l = cos(angle*pi/180)
    m = sin(angle*pi/180)

    # Rotation matrices
    if d != 0.0
        Rx  = [  1.0   0.0   0.0
                 0.0   c/d  -b/d
                 0.0   b/d   c/d ]

        Rxi = [  1.0   0.0   0.0
                 0.0   c/d   b/d
                 0.0  -b/d   c/d ]
    end

    Ry  = [   d    0.0  -a
             0.0    1.0  0.0
              a    0.0   d ]

    Ryi = [   d    0.0   a
             0.0    1.0  0.0
             -a    0.0   d ]

    Rz  = [   l   -m   0.0
              m    l   0.0
             0.0   0.0   1.0 ]

    # all rotations matrix
    if d != 0.0
        R = Rxi*Ryi*Rz*Ry*Rx
    else
        R = Ryi*Rz*Ry
    end

    coords =get_coords(bl.nodes)

    # equation: p2 = base + R*(p-base)
    coords = ( base .+ R*(coords' .- base) )'

    setcoords!(bl.nodes, coords)

    return bl
end

function rotate!(blocks::Array{T,1}; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 ) where T <: AbstractBlock
    for bl in blocks
        rotate!(bl, base=base, axis=axis, angle=angle)
    end
    return blocks
end


"""
`polar(block, [base=[0,0,0],] [axis=[0,0,1],] [angle=360.0,] [n=2])`

Creates `n-1` copies of a `block` and places them using polar distribution based on
a `base` point, an `axis` vector, a total `angle`.
"""
function polar(bl::T; base=[0.0,0,0], axis=[0.0,0,1], angle=360, n=2 ) where T <: AbstractBlock
    blocks::Array{T,1} = [ bl ]
    angle = angle/n
    for i=1:n-1
        bli = copy(bl)
        rotate!(bli, base=base, axis=axis, angle=angle*i)
        push!(blocks, bli)
    end
    return blocks
end

function polar(blocks::Array{T,1}; base=[0.0,0,0], axis=[0.0,0,1], angle=360, n=2 ) where T <: AbstractBlock
    rblocks::Array{T,1} = []

    for bl in blocks
        bls = polar(bl, base=base, axis=axis, angle=angle, n=n)
        append!(rblocks, bls)
    end

    return rblocks
end


"""
`move(mesh, [dx=0.0,] [dy=0.0,] [dz=0.0])`

Moves a Mesh object `mesh`. Also returns a reference.
"""
function move!(mesh::Mesh; dx=0.0, dy=0.0, dz=0.0)
    for p in mesh.nodes
        p.coord.x += dx
        p.coord.y += dy
        p.coord.z += dz
    end
    return mesh
end



function scale!(msh::Mesh; factor=1.0, base=[0.0,0,0])
    for p in msh.nodes
        p.coord.x, p.coord.y, p.coord.z = base + ([p.coord.x, p.coord.y, p.coord.z] - base)*factor
    end
    return msh
end


"""
`rotate(mesh; base=[0,0,0], axis=[0,0,1], angle=90.0)`

Rotates a Mesh object `mesh` according to a `base` point, an `axis` vector and an `angle`.
"""
function rotate!(mesh::Mesh; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )

    length(axis)==2 && ( axis=vcat(axis, 0.0) )
    length(base)==2 && ( base=vcat(base, 0.0) )
    norm(axis) < 1e-10 && error("rotate: Invalid axis $axis")

    # unit vector
    axis = axis/norm(axis)
    a, b, c = axis
    d = sqrt(b^2+c^2)
    d==0.0 && ( d=1.0 )

    # unit vector for rotation
    l = cos(angle*pi/180)
    m = sin(angle*pi/180)

    # Rotation matrices
    Rx  = [  1.0    0.0    0.0
             0.0   c/d  -b/d
             0.0   b/d   c/d ]

    Rxi = [  1.0    0.0    0.0
             0.0   c/d   b/d
             0.0  -b/d   c/d ]

    Ry  = [   d    0.0  -a
             0.0    1.0  0.0
              a    0.0   d ]

    Ryi = [   d    0.0   a
             0.0    1.0  0.0
             -a    0.0   d ]

    Rz  = [   l   -m   0.0
              m    l   0.0
             0.0   0.0   1.0 ]

    # all rotations matrix
    R = Rxi*Ryi*Rz*Ry*Rx

    if axis==[1.0, 0.0, 0.0]
        R  = [ 1.0   0.0   0.0
               0.0   l     -m
               0.0   m     l ]
    end

    for p in mesh.nodes
        p.coord.x, p.coord.y, p.coord.z = base + R*([p.coord.x, p.coord.y, p.coord.z] - base)
    end

    return mesh
end


# Roll Axes
# =========

function rollaxes!(bl::AbstractBlock)
    for p in bl.nodes
        p.coord.x, p.coord.y, p.coord.z = p.coord.z, p.coord.x, p.coord.y
    end
    return nothing
end

rollaxes!(bls::Array{<:AbstractBlock,1}) = (rollaxes!.(bls); nothing)

function rollaxes!(mesh::Mesh)
    if mesh.ndim==2
        for p in mesh.nodes
            p.coord.x, p.coord.y = p.coord.y, p.coord.x
        end
    else
        for p in mesh.nodes
            p.coord.x, p.coord.y, p.coord.z = p.coord.z, p.coord.x, p.coord.y
        end
    end

    if length(mesh.node_data)>0 || length(mesh.elem_data)>0
        notify("rollaxes!: mesh associated data was not reordered according to new axes.")
    end
end


function changeaxes!(bl::AbstractBlock, order::String)
    @assert length(order)==3
    idxs = [ char-'w' for char in order ]
    for p in bl.nodes
        p.coord[1:3] = p.coord[idxs]
    end
end

changeaxes!(bls::Array{<:AbstractBlock,1}, order::String) = (changeaxes!.(bls, order); nothing)


"""
    changeaxes!(mesh, order)

Changes the coordinates axes of a `mesh` according to a new `order` given as a string.

# Example

```
julia> mesh = Mesh(Block([0 0; 1 1], nx=2, ny=2));
julia> changeaxes!(mesh, "zxy")
```
"""
function changeaxes!(mesh::Mesh, order::String)
    @assert length(order)==3

    idxs = [ char-'w' for char in order ]
    for p in mesh.nodes
        p.coord[1:3] = p.coord[idxs]
    end

    for elem in mesh.elems
        isinverted(elem) && flip!(elem)
    end

    if length(mesh.node_data)>0 || length(mesh.elem_data)>0
        notify("changeaxes!: mesh associated data was not reordered according to new axes.")
    end
end


function isinverted(elem::AbstractCell)
    elem.shape.family==LINE_SHAPE && return false

    if elem.shape.ndim==2
        coords = get_coords(elem)
        #dx, dy, dz = std(coords, dims=1)
        dz = mean(abs, coords[:,3])
        tol = 1e-8
        dz>tol && return false

        if elem.shape.basic_shape == TRI3
            nidx = [1,2,3]
        elseif elem.shape.basic_shape == QUAD4
            nidx = [1,2,4]
        else
            error("isinverted: Cell shape $(cell.shape.name) is not supported")
        end

        X1, X2, X3 = getfield.(elem.nodes[nidx], :coord)
        V1 = X2-X1
        V2 = X3-X1
        V3 = Vec3(0,0,1)
    elseif elem.shape.ndim==3
        if elem.shape.basic_shape in (TET4,WED6)
            nidx = [1,2,3,4]
        elseif elem.shape.basic_shape in (PYR5,HEX8)
            nidx = [1,2,4,5]
        else
            error("isinverted: Cell shape $(cell.shape.name) is not supported")
        end

        X1, X2, X3, X4 = getfield.(elem.nodes[nidx], :coord)
        V1 = X2-X1
        V2 = X3-X1
        V3 = X4-X1
    else
        error("isinverted: Cell shape $(cell.shape.name) is not supported")
    end

    return dot(cross(V1,V2),V3)<0
end


function flip!(elem::AbstractCell)
    elem.shape.family==LINE_SHAPE && return elem

    if elem.shape==TRI3
        nidx = [1,3,2]
    elseif elem.shape==QUAD4
        nidx = [1,4,3,2]
    elseif elem.shape==TRI6
        nidx = [1,3,2,6,5,4]
    elseif elem.shape==QUAD8
        nidx = [1,4,3,2,8,7,6,5]
    elseif elem.shape==QUAD9
        nidx = [1,4,3,2,8,7,6,5,9]
    elseif elem.shape==TET4
        nidx = [1,3,2,4]
    elseif elem.shape==PYR5
        nidx = [1,4,3,2,5]
    elseif elem.shape==WED6
        nidx = [1,3,2,4,6,5]
    elseif elem.shape==HEX8
        nidx = [1,4,3,2,5,8,7,6]
    elseif elem.shape==TET10
        nidx = [1,3,2,4, 7,6,5, 8,10,9]
    elseif elem.shape==WED15
        nidx = [1,3,2,4,6,5,9,8,7,12,11,10,13,15,14]
    elseif elem.shape==HEX20
        nidx = [1,4,3,2,5,8,7,6, 12,11,10,9, 16,15,14,13, 17,20,19,18]
    else
        error("flip!: Cell shape $(cell.shape.name) is not supported")
    end

    elem.nodes = elem.nodes[nidx]
    return elem
end
