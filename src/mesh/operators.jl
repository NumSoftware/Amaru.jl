# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


"""
    $(SIGNATURES)

Moves a `block` position by updating its coordinates according
to `dx`, `dy` and `dz`.

# Examples

```julia
julia> block = Block([ 0 0 0; 1 1 1], nx=1, ny=2, nz=3)

julia> move!(block, dx=0.5, dy=1.0)

julia> block
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
    cellshape: CellShape  name="HEX8"
    nx: 1
    ny: 2
    nz: 3
    rx: 1.0
    ry: 1.0
    rz: 1.0
    tag: ""
```

"""
function move!(block::AbstractBlock; dx=0.0, dy=0.0, dz=0.0)
    for p in block.nodes
        p.coord.x += dx
        p.coord.y += dy
        p.coord.z += dz

        round!(p.coord, digits=8)
    end
    return block
end


"""
    $(SIGNATURES)

Moves each `Block` in the array `blocks` according to `dx`, `dy` and `dz`.
"""
function move!(blocks::Array{<:AbstractBlock,1}; dx=0.0, dy=0.0, dz=0.0)
    move!.(blocks, dx=dx, dy=dy, dz=dz)
    return blocks
end


"""
    $(SIGNATURES)

Scales a `block` from the point `base` using the given `factor`.
If `axis` is provided, the scaling is performent only in the 
`axis` direction.

# Examples

```julia
julia> block = Block([ 0 0; 1 1 ], nx=2, ny=2);

julia> getcoords(block)
4×3 Matrix{Float64}:
 0.0  0.0  0.0
 1.0  0.0  0.0
 1.0  1.0  0.0
 0.0  1.0  0.0

julia> scale!(block, factor=0.5, base=[ 0, 0 ], axis=[ 1, 0 ]);

julia> getcoords(block)
4×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.5  0.0  0.0
 0.5  0.0  0.0
 0.0  0.0  0.0
```
"""
function scale!(block::AbstractBlock; factor=1.0, base=[0,0,0], axis=nothing)
    coords = getcoords(block)
    base = Vec3(base)

    if axis===nothing
        coords = ( base .+ factor*(coords' .- base) )'
    else
        axis = Vec3(normalize(axis))
        coords = ( base .+ factor*axis*axis'*(coords' .- base))'
    end

    coords[abs.(coords) .< eps()] .= 0.0

    for (i,p) in enumerate(block.nodes)
        p.coord.x, p.coord.y, p.coord.z = coords[i,:]
    end

    return block
end

"""
    $(SIGNATURES)

Scales each `Block` in the array `blocks` according to `factor`, `base`.
If `axis` is provided, the scaling is performent only in the `axis` direction.
```
"""
function scale!(blocks::Array{<:AbstractBlock,1}; factor=1.0, base=[0.0,0,0], axis=nothing)
    for bl in blocks
        scale!(bl, factor=factor, base=base, axis=axis)
    end
    return blocks
end

# function mirror(block::AbstractBlock; face=[0.0 0 0; 0 1 0; 0 0 1])
#     nr, nc = size(face)
#     if nc==2
#         face = [ face zeros(nr) ]
#     end
#     if nr==2
#         face = [ face; [face[1,1] face[1,2] 1.0] ]
#     end

#     p1 = face[1,:]
#     p2 = face[2,:]
#     p3 = face[3,:]
#     normal = cross(p2-p1, p3-p1)
#     normal = normal/norm(normal)

#     bl = copy(block)
#     coords =getcoords(bl.nodes)

#     distances    = (coords .- p1')*normal       # d = n^.(xi - xp)
#     coords = coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

#     # fix coordinates in bl to keep anti-clockwise numbering
#     npts = size(coords)[1]
#     ndim = typeof(bl)==Block2D ? 2 : 3

#     if npts==8 && ndim==2
#         idxs = [ 4, 3, 2, 1, 7, 6, 5, 8 ]
#     elseif npts==8 && ndim==3
#         idxs = [ 5:8; 1:4 ]
#     elseif npts==20 && ndim==3
#         idxs = [ 5:8; 1:4; 13:16; 9:12; 17:20 ]
#     else
#         idxs = [ npts:-1:1; ]  # reverse
#     end
#     coords = coords[idxs,:]
#     bl.nodes = [ Node(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]

#     return bl
# end

"""
    $(TYPEDSIGNATURES)

Creates a new Block by mirroring `block` according to a plane defined by 
a normal `axis` and `base` point.

# Examples

```
julia> block1 = Block([ 0 0; 2 1; 1 1; 0 1], nx=2, ny=2);
julia> block2 = mirror(block1, axis=[1,1], base=[2,1])
Block
  nodes: 4-element Vector{Node}:
    1: Node  id=1
    2: Node  id=2
    3: Node  id=3
    4: Node  id=4
  shape: CellShape  name="QUAD4"
  cellshape: CellShape  name="QUAD4"
  nx: 2
  ny: 2
  nz: 1
  rx: 1.0
  ry: 1.0
  rz: 1.0
  tag: ""
```
"""
function mirror(block::AbstractBlock;  axis=[0.0, 0, 1], base=[0.0, 0, 0] )
    newblock = copy(block)
    axis = normalize(Vec3(axis))
    base = Vec3(base)
    digs = 8

    L = Vec3()
    X = Vec3()

    for node in newblock.nodes
        L .= node.coord .- base
        dist = dot(L, axis) # dist = n^.(xi - xp)

        X .= node.coord .- (2*dist).*axis .+ 0.0 # xi = xi - 2*d*n^
        node.coord .= round.(X, digits=digs) 
    end

    #! Inversion is required to generate valid cells (det>0)
    isinverted(newblock) && flip!(newblock)

   return newblock
end


function mirror(blocks::Array{<:AbstractBlock,1}; axis=[0.0, 0, 1], base=[0.0, 0, 0] )
    return [ mirror(bl, axis=axis, base=base) for bl in blocks ]
end


function mirror(mesh::Mesh; axis=[0.0, 0, 1], base=[0.0, 0, 0])
    axis = normalize(Vec3(axis))
    base = Vec3(base)

    # copy mesh
    newmesh = copy(mesh)

    # mirror
    L = Vec3()
    X = Vec3()
    for node in newmesh.nodes
        L = node.coord .- base
        dist = dot(L, axis) # dist = n^.(xi - xp)

        X .= node.coord .- (2*dist).*axis .+ 0.0 # xi = xi - 2*d*n^
        node.coord .= round.(X, digits=8) 
    end

    for elem in mesh.elems
        isinverted(elem) && flip!(elem)
    end

    fixup!(newmesh, reorder=false)
    return newmesh
end

# function mirror2(mesh::Mesh; face=[0.0 0 0; 0 1 0; 0 0 1])
#     nr, nc = size(face)
#     if nc==2
#         face = [ face zeros(nr) ]
#     end
#     if nr==2
#         face = [ face; [face[1,1] face[1,2] 1.0] ]
#     end
#     p1 = face[1,:]
#     p2 = face[2,:]
#     p3 = face[3,:]
#     normal = cross(p2-p1, p3-p1)
#     normal = normal/norm(normal)

#     # copy mesh
#     newmesh = copy(mesh)

#     # mirror
#     coords =getcoords(newmesh.nodes)
#     distances = (coords .- p1')*normal       # d = n^.(xi - xp)
#     coords    = coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

#     # updating points
#     for (i,p) in enumerate(newmesh.nodes)
#         p.coord.x = coords[i,1]
#         p.coord.y = coords[i,2]
#         p.coord.z = coords[i,3]
#     end

#     # fix connectivities
#     for c in newmesh.elems
#         if c.shape==HEX8
#             idxs = [ 5:8; 1:4 ]
#             c.nodes = c.nodes[idxs]
#         end
#     end

#     fixup!(newmesh)
#     return newmesh
# end


"""
    $(SIGNATURES)

Creates a rectangular array of blocks of size `nx`×`ny`×`nz`
using copies of `block` spaced at `dx`, `dy` and `dz`.
The original `bĺock` is considered as part of the result.

# Examples

```julia
julia> block = Block([ 0 0; 1 1 ], nx=3, ny=3);

julia> blocks = array(block, nx=2, ny=2, dx=1, dy=1);

julia> length(blocks)
4
    ```
"""
function array(block::AbstractBlock; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)
    blocks = [ block ]
    for k=0:nz-1
        for j=0:ny-1
            for i=0:nx-1
                i==j==k==0 && continue
                cp = copy(block)
                move!(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                push!(blocks, cp)
            end
        end
    end
    return blocks
end


"""
    $(TYPEDSIGNATURES)

Rotates `block` an `angle` (default 90 degrees) around an `axis` that passes by a `base` point.

# Examples

```
julia> block = Block([ 0 0; 1 1 ], nx=2, ny=2);

julia> getcoords(block)
4×3 Matrix{Float64}:
 0.0  0.0  0.0
 1.0  0.0  0.0
 1.0  1.0  0.0
 0.0  1.0  0.0

julia> rotate!(block, base=[ 0, 0, 0 ], axis=[ 0, 0, 1], angle=45);

julia> getcoords(block)
4×3 Matrix{Float64}:
 0.0       0.0       0.0
 0.707107  0.707107  0.0
 0.0       1.41421   0.0
 0.0       0.707107  0.0
```
"""
function LinearAlgebra.rotate!(bl::AbstractBlock; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )
    # see also: https://lucidar.me/en/quaternions/quaternions-rotations/

    axis = normalize(Vec3(axis))
    base = Vec3(base)
    θ    = angle*pi/180
    R    = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
    digs = 8

    X = Vec3()
    for node in bl.nodes
        X = base + R*(node.coord-base)*conj(R)
        node.coord .= round.(X, digits=digs) 
    end

    # isinverted(bl) && flip!(bl)

    return bl

    # # unit vector
    # axis = axis/norm(axis)
    # a, b, c = axis
    # d = sqrt(b^2+c^2)

    # # unit vector for rotation
    # l = cos(angle*pi/180)
    # m = sin(angle*pi/180)

    # # Rotation matrices
    # if d != 0.0
    #     Rx  = [  1.0   0.0   0.0
    #              0.0   c/d  -b/d
    #              0.0   b/d   c/d ]

    #     Rxi = [  1.0   0.0   0.0
    #              0.0   c/d   b/d
    #              0.0  -b/d   c/d ]
    # end

    # Ry  = [   d    0.0  -a
    #          0.0    1.0  0.0
    #           a    0.0   d ]

    # Ryi = [   d    0.0   a
    #          0.0    1.0  0.0
    #          -a    0.0   d ]

    # Rz  = [   l   -m   0.0
    #           m    l   0.0
    #          0.0   0.0   1.0 ]

    # # all rotations matrix
    # if d != 0.0
    #     R = Rxi*Ryi*Rz*Ry*Rx
    # else
    #     R = Ryi*Rz*Ry
    # end

    # coords =getcoords(bl.nodes)

    # # equation: p2 = base + R*(p-base)
    # coords = ( base .+ R*(coords' .- base) )'

    # setcoords!(bl.nodes, coords)

    # return bl
end

"""
    $(SIGNATURES)

    Rotates `blocks` an `angle` (default 90 degrees) around an `axis` 
    that passes by a `base` point. The elements in `blocks` can be `Block`
    objects or even lists of `Block` objects.
"""
function LinearAlgebra.rotate!(blocks::Array; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )
    rotate!.(blocks, base=base, axis=axis, angle=angle)
    return blocks
    # for bl in blocks
        # rotate!(bl, base=base, axis=axis, angle=angle)
    # end
    # return blocks
end


"""
    $(SIGNATURES)

    Creates a `polar` array using copies of `block` by rotating it
    around `axis` along an `angle` domain. `n` representes 
    the number of cells in the polar direction.

    # Examples

    ```julia
    julia> block = Block([ 0 0; 1 0; 0.707 0.707; 0 1; 0.5 0; 0.924 0.383; 0.382 0.924; 0.354 0.354 ], nx=3, ny=3);

    julia> blocks = polar(block, base=[ 0, 0 ,0 ], axis=[ 0, 0, 1 ], angle=360, n=4);

    julia> length(blocks)
    4
    ```
"""
function polar(block::AbstractBlock; base=[0.0,0,0], axis=[0.0,0,1], angle=360.0, n=2 )
    blocks = [ block ]
    angle = angle/n
    for i=1:n-1
        bli = copy(block)
        rotate!(bli, base=base, axis=axis, angle=angle*i)
        push!(blocks, bli)
    end
    return blocks
end


"""
    $(SIGNATURES)

    Creates a `polar` array from `blocks` by rotating it
    around `axis` along an `angle` domain. `n` representes 
    the number of cells in the polar direction.
    The elements in `blocks` can be `Block` objects or even lists of 
    `Block` objects.
"""
function polar(blocks::Array; base=[0.0,0,0], axis=[0.0,0,1], angle=360, n=2 )
    # flat array

    # blocks = mapreduce(x->isa(x, Array) ? flat(x) : [x], append!, blocks, init=[])
    return polar.(blocks, base=base, axis=axis, angle=angle, n=n)
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
function LinearAlgebra.rotate!(mesh::Mesh; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )

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

#function rollaxes!(bl::AbstractBlock)
#    for p in bl.nodes
#        p.coord.x, p.coord.y, p.coord.z = p.coord.z, p.coord.x, p.coord.y
#    end
#    return nothing
#end
#
#rollaxes!(bls::Array{<:AbstractBlock,1}) = (rollaxes!.(bls); nothing)
#
#function rollaxes!(mesh::Mesh)
#    if mesh.ndim==2
#        for p in mesh.nodes
#            p.coord.x, p.coord.y = p.coord.y, p.coord.x
#        end
#    else
#        for p in mesh.nodes
#            p.coord.x, p.coord.y, p.coord.z = p.coord.z, p.coord.x, p.coord.y
#        end
#    end
#
#    if length(mesh.node_data)>0 || length(mesh.elem_data)>0
#        notify("rollaxes!: mesh associated data was not reordered according to new axes.")
#    end
#end


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

    fixup!(mesh, reorder=false)

    if length(mesh.node_data)>0 || length(mesh.elem_data)>0
        notify("changeaxes!: mesh associated data was not reordered according to new axes.")
    end
end


function isinverted(elem::AbstractCell)
    elem.shape.family==LINE_CELL && return false

    if elem.shape.ndim==2
        coords = getcoords(elem)
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
    elem.shape.family==LINE_CELL && return elem

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
