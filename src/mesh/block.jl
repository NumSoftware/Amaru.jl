# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function box_coords(C1::Array{<:Real,1}, C2::Array{<:Real,1}, ndim::Int)
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if ndim==2
        return [
                 x1     y1     0.0
                 x1+lx  y1     0.0
                 x1+lx  y1+ly  0.0
                 x1     y1+ly  0.0
               ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
                 x1     y1     z1
                 x1+lx  y1     z1
                 x1+lx  y1+ly  z1
                 x1     y1+ly  z1
                 x1     y1     z1+lz
                 x1+lx  y1     z1+lz
                 x1+lx  y1+ly  z1+lz
                 x1     y1+ly  z1+lz
                ]
    end
end


"""
    Block

A type that represents a segment, area or volume and is used to
aid the generation of structured meshes by subdivision.

# Fields
$(FIELDS)
"""
mutable struct Block <: AbstractBlock
    ndim::Int
    points::Vector{Point}
    shape::CellShape
    cellshape::CellShape
    nx::Int64
    ny::Int64
    nz::Int64
    rx::Float64
    ry::Float64
    rz::Float64
    tag::String
    id::Int64

    function Block(
        points::Vector{Point};
        nx::Int  = 0,
        ny::Int  = 0,
        nz::Int  = 0,
        n ::Int  = 0,
        rx::Real = 1.0,
        ry::Real = 1.0,
        rz::Real = 1.0,
        r ::Real = 0.0,
        cellshape = nothing,
        tag       = "",
        shape     = nothing,
        )

        if cellshape===nothing && shape!==nothing
            # notify("Block: argument shape was deprecated. Please use cellshape instead")
            cellshape = shape
        end

        shapes1d = (LIN2, LIN3, LIN4)
        shapes2d = (TRI3, TRI6, QUAD4, QUAD8, QUAD9, QUAD12)
        shapes3d = (TET4, TET10, HEX8, HEX20, HEX27, PYR5)

        # Get ndim
        sumy = sum(abs, [ p.coord.y for p in points ])
        sumz = sum(abs, [ p.coord.z for p in points ])

        ndim = 3
        n>0 && (nx=n)
        r>0 && (rx=r)
        sumz==0 && (ndim=2)
        sumy+sumz==0 && (ndim=1)

        # Check for surface or chord
        surface = ndim==3 && nz==0
        chord   = ndim>1 && ny==0 && nz==0
        chord

        nz==0 && ndim==3 && (nz=1)
        ny==0 && ndim>=2 && (ny=1)
        cellshape in shapes3d && (ndim==3 || error("Block: 3d points and nx, ny and nz are required for cell shape $(cellshape.name)"))

        npoints = length(points)

        if ndim==1 || chord
            npoints in (2, 3) || error("Block: invalid number of points ($npoints) for dimension $ndim or chord.")
            cellshape===nothing && (cellshape=LIN2)
            cellshape in shapes1d || error("Block: invalid cell type $(cellshape.name) for dimension $ndim.")
            shape = npoints==2 ? LIN2 : LIN3
        elseif ndim==2 || surface
            npoints in (4, 8) || error("Block: invalid number of points ($npoints) for dimension $ndim or surface.")
            cellshape===nothing && (cellshape=QUAD4)
            cellshape in shapes2d || error("Block: invalid cell type $(cellshape.name) for dimension $ndim or surface.")
            shape = npoints==4 ? QUAD4 : QUAD8
        else
            npoints in (8, 20) || error("Block: invalid number of points ($npoints) for dimension $ndim.")
            cellshape===nothing && (cellshape=HEX8)
            cellshape in shapes3d || error("Block: invalid cell type $(cellshape.name) for dimension $ndim.")
            shape = npoints==8 ? HEX8 : HEX20
        end

        for i in 1:length(points)
            points[i].id = i
        end

        return new(ndim, points, shape, cellshape, nx, ny, nz, rx, ry, rz, tag)
    end


    function Block(coords::Matrix{<:Real}; args...)
        nz = get(args, :nz, 0)
        ny = get(args, :ny, 0)

        ncoord, ncol = size(coords)
        ncol<=3 || error("Block: invalid coordinate matrix")

        # Get ndim
        sumy = ncol>=2 ? sum(abs, coords[:,2]) : 0.0
        sumz = ncol==3 ? sum(abs, coords[:,3]) : 0.0

        ndim = 3
        sumz==0 && (ndim=2)
        sumy+sumz==0 && (ndim=1)

        # Check for surface or chord
        surface = ndim==3 && nz==0
        chord   = ndim>1 && ny==0 && nz==0

        if ndim in (2,3) && ncoord==2 && !surface && !chord
            coords = box_coords(coords[1,:], coords[2,:], ndim)
            ncoord = size(coords,1)
        end
        points = [ Point(coords[i,:]) for i in 1:ncoord ]

        return Block(points; args...)
    end
end


"""
    $(TYPEDSIGNATURES)

Creates a copy of `block`.
"""
function Base.copy(block::Block)

    return Block(copy(getcoords(block.points)), nx=block.nx, ny=block.ny, nz=block.nz, cellshape=block.cellshape, tag=block.tag)
    # return Block(copy(block.points), nx=block.nx, ny=block.ny, nz=block.nz, cellshape=block.cellshape, tag=block.tag)
end


"""
    $(SIGNATURES)

Creates a copy of the array `blocks` containing `Block` objects.
"""
function Base.copy(blocks::Array{<:AbstractBlock,1})
    return [ copy(bl) for bl in blocks ]
end


function BlockGrid(
    X::Array{<:Real},            # list of x coordinates
    Y::Array{<:Real},            # list of y coordinates
    Z::Array{<:Real}=Float64[];  # list of z coordinates
    nx=[],  # list of divisions in the x direction
    ny=[],  # list of divisions in the y direction
    nz=[],  # list of divisions in the z direction
    rx=[],  # list of divisions ratios in the x direction
    ry=[],  # list of divisions ratios in the x direction
    rz=[],  # list of divisions ratios in the x direction
    cellshape=QUAD4, # element shape
    tag="",          # elements tag
    id=-1
    )

    length(rx)==0 && (rx = ones(length(nx)))
    length(ry)==0 && (ry = ones(length(ny)))
    length(rz)==0 && (rz = ones(length(nz)))
    blocks = Block[]
    for i in 1:length(nx)
        for j in 1:length(ny)
            coords = [
                X[i] Y[j]
                X[i+1] Y[j+1]
            ]
            bl = Block(coords, cellshape=cellshape, nx=nx[i], ny=ny[j], rx=rx[i], ry=ry[j])
            push!(blocks, bl)
        end
    end
    return blocks
end

