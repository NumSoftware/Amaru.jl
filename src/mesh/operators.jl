# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh


#Base.copy(p::Point) = Point(p.x, p.y, p.z, p.tag)
Base.copy(c::Cell)  = Cell(c.shape, c.points, tag=c.tag, ocell=c.ocell, nips=c.nips)

function Base.copy(bls::Array{<:AbstractBlock,1})
    return [ copy(obj) for obj in bls ]
end


"""
`move(block, [dx=0.0,] [dy=0.0,] [dz=0.0])`

Changes the coordinates of a `block`. Also returns a reference.
"""
function move!(bl::AbstractBlock; dx=0.0, dy=0.0, dz=0.0)
    for p in bl.points
        p.x += dx
        p.y += dy
        p.z += dz
    end
    return bl
end


"""
`move(blocks, [dx=0.0,] [dy=0.0,] [dz=0.0])`

Changes the coordinates of an array of blocks. Also returns a reference.
"""
function move!(blocks::Array; dx=0.0, dy=0.0, dz=0.0)
    for bl in blocks
        for p in bl.points
            p.x += dx
            p.y += dy
            p.z += dz
        end
    end
    return blocks
end


function scale!(bl::AbstractBlock; factor=1.0, base=[0.,0,0], axis=nothing)
    coords = getcoords(bl)
    coords = ( base .+ factor*(coords' .- base) )'

    for (i,p) in enumerate(bl.points)
        p.x, p.y, p.z = coords[i,:]
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
    coords = getcoords(bl.points)

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
    bl.points = [ Point(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]

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
    coords = getcoords(newmesh.points)
    distances = (coords .- p1')*normal       # d = n^.(xi - xp)
    coords    = coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

    # updating points
    for (i,p) in enumerate(newmesh.points)
        p.x = coords[i,1]
        p.y = coords[i,2]
        p.z = coords[i,3]
    end

    # fix connectivities
    for c in newmesh.cells
        if c.shape==HEX8
            idxs = [ 5:8; 1:4 ]
            c.points = c.points[idxs]
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
        Rx  = [  1.0    0.0    0.0
                 0.0   c/d  -b/d
                 0.0   b/d   c/d ]

        Rxi = [  1.0    0.0    0.0
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

    coords = getcoords(bl.points)

    # equation: p2 = base + R*(p-base)
    coords = ( base .+ R*(coords' .- base) )'

    setcoords!(bl.points, coords)

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
    for p in mesh.points
        p.x += dx
        p.y += dy
        p.z += dz
    end
    return mesh
end



function scale!(msh::Mesh; factor=1.0, base=[0.0,0,0])
    for p in msh.points
        p.x, p.y, p.z = base + ([p.x, p.y, p.z] - base)*factor
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
    #@show d
    d==0.0 && ( d=1.0 )
    #@show d

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

    for p in mesh.points
        p.x, p.y, p.z = base + R*([p.x, p.y, p.z] - base)
    end

    return mesh
end


# Roll Axes
# =========

function rollaxes!(bl::AbstractBlock)
    for p in bl.points
        p.x, p.y, p.z = p.z, p.x, p.y
    end
    return nothing
end

rollaxes!(bls::Array{<:AbstractBlock,1}) = (rollaxes!.(bls); nothing)

function rollaxes!(mesh::Mesh)
    if mesh.ndim==2
        for p in mesh.points
            p.x, p.y = p.y, p.x
        end
    else
        for p in mesh.points
            p.x, p.y, p.z = p.z, p.x, p.y
        end
    end
end


# TESTING
function get_surface_alt(cells::Array{Cell,1})
    # Actually slower....
    # Get all points
    pointsd = Dict{UInt64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    points = values(pointsd)

    # Get incidence matrix (shares) (fast)
    np = length(points)
    N = [ Cell[] for i=1:np]
    for cell in cells
        for pt in cell.points
            push!(N[pt.id], cell)
        end
    end

    @show 1
    # Get matrix of cells faces
    F = [ get_faces(cell) for cell in cells]
    nc = length(cells)
    #CF = Array(Array{Array{Int64,1},1}, nc)
    CF = Array(Array{UInt64,1}, nc)
    for cell in cells # fast
        #CF[cell.id] = [ sort([pt.id for pt in face.points]) for face in F[cell.id]]
        CF[cell.id] = [ hash(face) for face in F[cell.id]]
    end

    @show 2
    # Get cells boundary flag matrix
    CB = [ trues(length(CF[cell.id])) for cell in cells]
    for cell in cells
        for (i,fcon) in enumerate(CF[cell.id])
            for pid in fcon
                for cl in N[pid]
                    if cl.id == cell.id; continue end
                    if fcon in CF[cl.id]
                        CB[cell.id][i] = false
                    end
                end
            end
        end
    end

    @show 3
    # Get list of boundary faces (almost fast)
    facets = Cell[]
    for cell in cells
        for (i,face) in enumerate(F[cell.id])
            if CB[cell.id][i]
                push!(facets, face)
            end
        end
    end
    #return facets
end
