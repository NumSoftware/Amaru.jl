# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


function box_points(C1::Array{<:Real,1}, C2::Array{<:Real,1})
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if length(C1)==2
        return [
                 Node(x1   , y1   , 0.0),
                 Node(x1+lx, y1   , 0.0),
                 Node(x1+lx, y1+ly, 0.0),
                 Node(x1   , y1+ly, 0.0),
                ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
                 Node(x1   , y1   , z1 ),
                 Node(x1+lx, y1   , z1 ),
                 Node(x1+lx, y1+ly, z1 ),
                 Node(x1   , y1+ly, z1 ),
                 Node(x1   , y1   , z1+lz ),
                 Node(x1+lx, y1   , z1+lz ),
                 Node(x1+lx, y1+ly, z1+lz ),
                 Node(x1   , y1+ly, z1+lz ),
                ]
    end
end


"""
`Block(coords, [nx=1,] [ny=1,] [cellshape=QUAD4,] [tag=""] )`

Generates a block object for the mesh generation of 2D meshes.
`shape` can be TRI3, TRI6, QUAD4, QUAD8.
"""
mutable struct Block <: AbstractBlock
    nodes::Array{Node,1}
    shape::ShapeType
    cellshape::ShapeType
    nx::Int64
    ny::Int64
    nz::Int64
    tag::String
    nips::Int64
    id::Int64

    function Block(coords::Array{<:Real}; nx::Int=1, ny::Int=1, nz::Int=1, cellshape=nothing, tag="", id=-1, shape=nothing, nips=0)
        if shape != nothing
            @warn "Block: argument shape was deprecated. Please use cellshape instead"
            cellshape = shape
        end

        shapes1d = (LIN2, LIN3)
        shapes2d = (TRI3, TRI6, QUAD4, QUAD8, QUAD9, QUAD12)
        shapes3d = (TET4, TET10, HEX8, HEX20, HEX27)

        ncoord, ndim = size(coords)
        ndim<=3 || error("Block: invalid coordinate matrix")
        (ndim==3 && cellshape in shapes2d) && (ndim=2)
        cellshape in shapes3d && (ndim==3 || error("Block: 3d points are required for cell shape $cellshape"))

        if ndim==1
            ncoord==2 || error("Block: invalid coordinates matrix rows ($ncoord) for dimension $ndim")
            cellshape==nothing && (cellshape=LIN2)
            cellshape in shapes1d || error("Block: invalid cell type $cellshape for dimension $ndim")
            shape = LIN2
            nodes = box_points(coords[1,:], coords[2,:])
        elseif ndim==2
            ncoord in (2, 4, 8) || error("Block: invalid coordinates matrix rows ($ncoord) for dimension $ndim")
            cellshape==nothing && (cellshape=QUAD4)
            cellshape in shapes2d || error("Block: invalid cell type $cellshape for dimension $ndim")
            if ncoord==2
                nodes = box_points(coords[1,:], coords[2,:])
                shape = QUAD4
            else
                nodes = [ Node(coords[i,:]) for i=1:ncoord ]
                shape = length(nodes)==4 ? QUAD4 : QUAD8
            end
        else
            ncoord in (2, 8, 20) || error("Block: invalid coordinates matrix rows ($ncoord) for dimension $ndim")
            cellshape==nothing && (cellshape=HEX8)
            cellshape in shapes3d || error("Block: invalid cell type $cellshape for dimension $ndim")
            if ncoord==2
                nodes = box_points(coords[1,:], coords[2,:])
                shape = HEX8
            else
                nodes = [ Node(coords[i,:]) for i=1:ncoord ]
                shape = length(nodes)==8 ? HEX8 : HEX20
            end
        end

        for i=1:length(nodes)
            nodes[i].id = i
        end

        return new(nodes, shape, cellshape, nx, ny, nz, tag, nips, id)
    end
end

# For backward compatibility
Block2D = Block
Block3D = Block


function Base.copy(bl::Block; dx=0.0, dy=0.0, dz=0.0)
    newbl = Block3D(copy(get_coords(bl.nodes)), nx=bl.nx, ny=bl.ny, nz=bl.nz, cellshape=bl.cellshape, tag=bl.tag, nips=bl.nips)
end


# Splits a block
# TODO: replace msh::Mesh by points, pointdict and cells
# TODO: optimize matrix products
function split_block(bl::Block, msh::Mesh)
    #nx, ny = bl.nx, bl.ny
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    shape  = bl.shape # cell shape
    coords = get_coords(bl.nodes)
    cellshape = bl.cellshape

    if cellshape==QUAD4
        p_arr = Array{Node}(undef, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, digits=8)
                    p =get_node(msh._pointdict, C)
                    if p==nothing
                        p = Node(C); push!(msh.nodes, p)
                        msh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C);
                    push!(msh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                cell = Cell(cellshape, [p1, p2, p3, p4], tag=bl.tag, nips=bl.nips)
                push!(msh.elems, cell)
            end
        end
        return
    end

    if cellshape == QUAD8 || cellshape == QUAD9
        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                if cellshape==QUAD8 && iseven(i) && iseven(j) continue end

                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(msh._pointdict, C)
                    if p==nothing
                        p = Node(C); push!(msh.nodes, p)
                        msh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(msh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                if cellshape==QUAD8
                    cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag, nips=bl.nips)
                else
                    p9   = p_arr[i+1, j+1]
                    cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9], tag=bl.tag, nips=bl.nips)
                end
                push!(msh.elems, cell)
            end
        end
        return
    end

    if cellshape == QUAD12
        p_arr = Array{Node}(undef, 3*nx+1, 3*ny+1)
        for j = 1:3*ny+1
            for i = 1:3*nx+1
                if cellshape==QUAD12 && (i-1)%3>0 && (j-1)%3>0 continue end

                r = ((2/3)/nx)*(i-1) - 1.0
                s = ((2/3)/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 3*nx+1) || j in (1, 3*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(msh._pointdict, C)
                    if p==nothing
                        p = Node(C); push!(msh.nodes, p)
                        msh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(msh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:3:3*ny
            for i = 1:3:3*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+3, j  ]
                p3 = p_arr[i+3, j+3]
                p4 = p_arr[i  , j+3]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+3, j+1]
                p7 = p_arr[i+2, j+3]
                p8 = p_arr[i  , j+2]

                p9  = p_arr[i+2, j  ]
                p10 = p_arr[i+3, j+2]
                p11 = p_arr[i+1, j+3]
                p12 = p_arr[i  , j+1]

                cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], tag=bl.tag, nips=bl.nips)
                push!(msh.elems, cell)
            end
        end
        return
    end

    if cellshape == TRI3
        p_arr = Array{Node}(undef, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, digits=8)
                    p =get_node(msh._pointdict, C)
                    if p==nothing
                        p = Node(C); push!(msh.nodes, p)
                        msh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(msh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                cell1 = Cell(cellshape, [p1, p2, p3], tag=bl.tag, nips=bl.nips)
                cell2 = Cell(cellshape, [p4, p1, p3], tag=bl.tag, nips=bl.nips)
                push!(msh.elems, cell1)
                push!(msh.elems, cell2)
            end
        end
        return
    end

    if cellshape == TRI6

        #=   4       7       3
               @-----@-----@
               |         / |
               |       /   |
             8 @     @     @ 6
               |   /  9    |
               | /         |
               @-----@-----@
             1       5       2     =#

        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(msh._pointdict, C)
                    if p==nothing
                        p = Node(C); push!(msh.nodes, p)
                        msh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(msh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                p9   = p_arr[i+1, j+1]

                cell1 = Cell(cellshape, [p1, p2, p3, p5, p6, p9], tag=bl.tag, nips=bl.nips)
                cell2 = Cell(cellshape, [p4, p1, p3, p8, p9, p7], tag=bl.tag, nips=bl.nips)
                push!(msh.elems, cell1)
                push!(msh.elems, cell2)
            end
        end
        return
    end

    if cellshape==HEX8 || cellshape==TET4
        p_arr = Array{Node}(undef, nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        C = round.(C, digits=8)
                        p =get_node(msh._pointdict, C)
                        if p==nothing
                            p = Node(C); push!(msh.nodes, p)
                            msh._pointdict[hash(p)] = p
                        end
                    else
                        p = Node(C); push!(msh.nodes, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    p1 = p_arr[i  , j  , k  ]
                    p2 = p_arr[i+1, j  , k  ]
                    p3 = p_arr[i+1, j+1, k  ]
                    p4 = p_arr[i  , j+1, k  ]
                    p5 = p_arr[i  , j  , k+1]
                    p6 = p_arr[i+1, j  , k+1]
                    p7 = p_arr[i+1, j+1, k+1]
                    p8 = p_arr[i  , j+1, k+1]

                    if cellshape==HEX8
                        cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag, nips=bl.nips)
                        push!(msh.elems, cell)
                    end
                    if cellshape==TET4
                        push!( msh.elems, Cell(cellshape, [p2, p4, p1, p8], tag=bl.tag, nips=bl.nips) )
                        push!( msh.elems, Cell(cellshape, [p2, p1, p5, p8], tag=bl.tag, nips=bl.nips) )
                        push!( msh.elems, Cell(cellshape, [p2, p5, p6, p8], tag=bl.tag, nips=bl.nips) )
                        push!( msh.elems, Cell(cellshape, [p2, p6, p7, p8], tag=bl.tag, nips=bl.nips) )
                        push!( msh.elems, Cell(cellshape, [p2, p3, p4, p8], tag=bl.tag, nips=bl.nips) )
                        push!( msh.elems, Cell(cellshape, [p2, p7, p3, p8], tag=bl.tag, nips=bl.nips) )
                    end
                end
            end
        end
        return
    end

    #if cellshape == HEX20 || cellshape == TET10 # || cellshape == HEX27
    if cellshape in (TET10, HEX20, HEX27)
        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if cellshape==HEX20
                        if iseven(i) && iseven(j) continue end
                        if iseven(j) && iseven(k) continue end
                        if iseven(k) && iseven(i) continue end
                    end

                    r = (1.0/nx)*(i-1) - 1.0
                    s = (1.0/ny)*(j-1) - 1.0
                    t = (1.0/nz)*(k-1) - 1.0
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        C = round.(C, digits=8)
                        p =get_node(msh._pointdict, C)
                        if p==nothing
                            p = Node(C); push!(msh.nodes, p)
                            msh._pointdict[hash(p)] = p
                        end
                    else
                        p = Node(C); push!(msh.nodes, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:2:2*nz
            for j = 1:2:2*ny
                for i = 1:2:2*nx
                    p1  = p_arr[i  , j  , k  ]
                    p2  = p_arr[i+2, j  , k  ]
                    p3  = p_arr[i+2, j+2, k  ]
                    p4  = p_arr[i  , j+2, k  ]
                    p5  = p_arr[i  , j  , k+2]
                    p6  = p_arr[i+2, j  , k+2]
                    p7  = p_arr[i+2, j+2, k+2]
                    p8  = p_arr[i  , j+2, k+2]

                    p9  = p_arr[i+1, j  , k  ]
                    p10 = p_arr[i+2, j+1, k  ]
                    p11 = p_arr[i+1, j+2, k  ]
                    p12 = p_arr[i  , j+1, k  ]
                    p13 = p_arr[i+1, j  , k+2]
                    p14 = p_arr[i+2, j+1, k+2]
                    p15 = p_arr[i+1, j+2, k+2]
                    p16 = p_arr[i  , j+1, k+2]

                    p17 = p_arr[i  , j  , k+1]
                    p18 = p_arr[i+2, j  , k+1]
                    p19 = p_arr[i+2, j+2, k+1]
                    p20 = p_arr[i  , j+2, k+1]

                    #conn = [
                        #p_arr[i  , j  , k  ],
                        #p_arr[i+2, j  , k  ],
                        #p_arr[i+2, j+2, k  ],
                        #p_arr[i  , j+2, k  ],
                        #p_arr[i  , j  , k+2],
                        #p_arr[i+2, j  , k+2],
                        #p_arr[i+2, j+2, k+2],
                        #p_arr[i  , j+2, k+2],
#
                        #p_arr[i+1, j  , k  ],
                        #p_arr[i+2, j+1, k  ],
                        #p_arr[i+1, j+2, k  ],
                        #p_arr[i  , j+1, k  ],
                        #p_arr[i+1, j  , k+2],
                        #p_arr[i+2, j+1, k+2],
                        #p_arr[i+1, j+2, k+2],
                        #p_arr[i  , j+1, k+2],
#
                        #p_arr[i  , j  , k+1],
                        #p_arr[i+2, j  , k+1],
                        #p_arr[i+2, j+2, k+1],
                        #p_arr[i  , j+2, k+1]]

                    if cellshape == HEX20
                        cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20], tag=bl.tag, nips=bl.nips)
                        push!(msh.elems, cell)
                    end
                    if cellshape in (TET10, HEX27)

                        #p1  = p_arr[i  , j  , k  ]
                        #p2  = p_arr[i+2, j  , k  ]
                        #p3  = p_arr[i+2, j+2, k  ]
                        #p4  = p_arr[i  , j+2, k  ]
                        #p5  = p_arr[i  , j  , k+2]
                        #p6  = p_arr[i+2, j  , k+2]
                        #p7  = p_arr[i+2, j+2, k+2]
                        #p8  = p_arr[i  , j+2, k+2]
#
                        #p9  = p_arr[i+1, j  , k  ]
                        #p10 = p_arr[i+2, j+1, k  ]
                        #p11 = p_arr[i+1, j+2, k  ]
                        #p12 = p_arr[i  , j+1, k  ]
                        #p13 = p_arr[i+1, j  , k+2]
                        #p14 = p_arr[i+2, j+1, k+2]
                        #p15 = p_arr[i+1, j+2, k+2]
                        #p16 = p_arr[i  , j+1, k+2]
#
                        #p17 = p_arr[i  , j  , k+1]
                        #p18 = p_arr[i+2, j  , k+1]
                        #p19 = p_arr[i+2, j+2, k+1]
                        #p20 = p_arr[i  , j+2, k+1]

                        #p21 = p_arr[i+1, j+1, k  ]
                        #p22 = p_arr[i+1, j+1, k+2]
                        #p23 = p_arr[i+1, j  , k+1]
                        #p24 = p_arr[i+2, j+1, k+1]
                        #p25 = p_arr[i+1, j+2, k+1]
                        #p26 = p_arr[i  , j+1, k+1]
                        #p27 = p_arr[i+1, j+1, k+1]


                        p21 = p_arr[i  , j+1, k+1]
                        p22 = p_arr[i+2, j+1, k+1]
                        p23 = p_arr[i+1, j  , k+1]
                        p24 = p_arr[i+1, j+2, k+1]
                        p25 = p_arr[i+1, j+1, k  ]
                        p26 = p_arr[i+1, j+1, k+2]
                        p27 = p_arr[i+1, j+1, k+1]

                        if cellshape==TET10
                            push!( msh.elems, Cell(cellshape, [p2, p4, p1, p8, p25, p12, p9, p27, p20, p21], tag=bl.tag, nips=bl.nips) )
                            push!( msh.elems, Cell(cellshape, [p2, p1, p5, p8, p9, p17, p23, p27, p21, p16], tag=bl.tag, nips=bl.nips) )
                            push!( msh.elems, Cell(cellshape, [p2, p5, p6, p8, p23, p13, p18, p27, p16, p26], tag=bl.tag, nips=bl.nips) )
                            push!( msh.elems, Cell(cellshape, [p2, p6, p7, p8, p18, p14, p22, p27, p26, p15], tag=bl.tag, nips=bl.nips) )
                            push!( msh.elems, Cell(cellshape, [p2, p3, p4, p8, p10, p11, p25, p27, p24, p20], tag=bl.tag, nips=bl.nips) )
                            push!( msh.elems, Cell(cellshape, [p2, p7, p3, p8, p22, p19, p10, p27, p15, p24], tag=bl.tag, nips=bl.nips) )
                        else
                            cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27], tag=bl.tag, nips=bl.nips)
                            push!(msh.elems, cell)
                        end
                    end
                end
            end
        end
        return
    end
    error("block: Can not discretize using shape $(cellshape.name)")
end


mutable struct BlockCylinder <: AbstractBlock
    nodes::Array{Node,1}
    shape::ShapeType # LIN2
    cellshape::ShapeType # HEX8, HEX20
    r::Float64
    nr::Int64
    n::Int64
    tag::String
    nips::Int64
    id::Int64

    function BlockCylinder(coords::Array{<:Real}; r=1.0, nr=3, n=2, cellshape=HEX8, tag="", nips=0, id=-1)
        size(coords,1) != 2 && error("Invalid coordinates matrix for BlockCylinder")
        nr<2 && error("Invalid nr=$nr value for BlockCylinder")
        cellshape in (HEX8, HEX20) || error("BlockCylinder: cellshape must be HEX8 or HEX20")
        nodes = [ Node(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]
        return new(nodes, LIN2, cellshape, r, nr, n, tag, nips, id)
    end
end


function Base.copy(bl::BlockCylinder)
    newbl = BlockCylinder(copy(get_coords(bl.nodes)), r=bl.r, nr=bl.nr, n=bl.n, cellshape=bl.cellshape, tag=bl.tag, nips=bl.nips)
end


function split_block(bl::BlockCylinder, msh::Mesh)

    nx1 = round(Int, bl.nr/3)
    nx2 = bl.nr - nx1
    shape2D = bl.cellshape==HEX8 ? QUAD4 : QUAD8

    # constructing quadratic blocks
    coords = bl.r*[ 0 0; 1/3 0; 1/3 1/3; 0 1/3; 1/6 0; 1/3 1/6; 1/6 1/3; 0 1/6 ]
    bl1 = Block(coords, nx=nx1, ny=nx1, cellshape= shape2D, tag=bl.tag, nips=bl.nips)

    s45  = sin(45*pi/180)
    c45  = s45
    s225 = sin(22.5*pi/180)
    c225 = cos(22.5*pi/180)

    coords = bl.r*[ 1/3 0; 1 0; c45 s45; 1/3 1/3; 2/3 0; c225 s225; (c45+1/3)/2 (s45+1/3)/2; 1/3 1/6 ]
    bl2    = Block(coords, nx=nx2, ny=nx1, cellshape= shape2D, tag=bl.tag, nips=bl.nips)

    coords = bl.r*[ 0 1/3; 1/3 1/3; c45 s45; 0 1; 1/6 1/3; (c45+1/3)/2 (s45+1/3)/2; s225 c225; 0 2/3 ]
    bl3    = Block(coords, nx=nx1, ny=nx2, cellshape= shape2D, tag=bl.tag, nips=bl.nips)

    blocks = [bl1, bl2, bl3 ]

    # polar and move
    blocks = polar(blocks, n=4)
    coords =get_coords(bl.nodes)
    move!(blocks, dx=coords[1,1], dy=coords[1,2], dz=coords[1,3])

    # extrude
    len      = norm(coords[1,:] - coords[2,:])
    blocks3D = extrude(blocks, len=len, n=bl.n)

    # rotation
    zv    = [0.0, 0.0, 1.0]
    axis  = coords[2,:] - coords[1,:]
    angle = acos( dot(zv, axis)/(norm(zv)*norm(axis)) )*180/pi
    raxis = cross(zv, axis)
    norm(raxis)>1e-10 && rotate!(blocks3D, base=coords[1,:], axis=raxis, angle=angle)

    # split
    for bl in blocks3D
        split_block(bl, msh)
    end

end


mutable struct BlockGrid<: AbstractBlock
    points::Array{Node,1}
    shape::ShapeType # QUAD4
    cellshape::ShapeType # QUAD4, QUAD8
    r::Float64
    nr::Int64
    n::Int64
    tag::String
    nips::Int64
    id::Int64

    function BlockGrid(coords::Array{<:Real}; r=1.0, nr=3, n=2, cellshape=HEX8, tag="", nips=0, id=-1)
        # TODO
        #return new(points, LIN2, cellshape, r, nr, n, tag, nips, id)
    end
end

