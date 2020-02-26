# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh


export triangulate
import Base.repr
import Base.show

# Point class for triangulation
mutable struct TPoint
    x::Float64
    y::Float64
    id::Int64
    function TPoint(x,y)
        this = new(x,y)
        this.id = -1
        return this
    end
end

# Cell class for triangulation
mutable struct TCell
    points::Array{TPoint,1}
    adjacs::Array{Union{TCell,Nothing},1}

    function TCell(p0, p1, p2, t1=nothing, t2=nothing, t3=nothing)
        this = new()
        this.points = [p0, p1, p2]
        this.adjacs = [t1, t2, t3]
        return this
    end

end

function repr(cell::TCell)
    repr(cell.points)
end

function show(io::IO, cell::TCell)
    print(repr(cell))
end

# Function to increase the index of a triangle vertex
function inc(v::Int64, n::Int64)
    return v+n>3 ? (v+n)%3 : v+n
end

# Recursive algorithm to find a triangle cell
function find_cell(x::Float64, y::Float64, cell::TCell)
    for i=1:3
        pa = cell.points[i]
        pb = cell.points[inc(i,1)]
        if (pb.x-pa.x)*(y-pa.y) - (pb.y-pa.y)*(x-pa.x) < 0.0 # point lies on right side
            acell = cell.adjacs[i]
            if acell==nothing
                println("find_cell Error: Looking for a point outside the mesh")
            end
            return find_cell(x, y, acell)
        end
    end
    return cell
end

# Function to swap common edge of two triangles according to Delaunay triangulation
function swap_cells(cell0, cell1)
    # search for p0 and v0 v1
    local p0, v0, v1 ::TPoint
    for i=1:3
        p0 = cell0.points[i]
        if !(p0 in cell1.points)
            v0 = cell0.points[inc(i,2)]
            v1 = cell0.points[inc(i,1)]
            break
        end
    end

    pos_p0 = findfirst(cell0.points, p0)

    # Get p1
    pos_p1 = inc(findfirst(cell1.points, v1), 1)
    p1  = cell1.points[pos_p1]

    # Calculate if p0 from cell 0 is in circuncircle of cell1
    x01 = v0.x - p1.x  ;  y01 = v0.y - p1.y
    x11 = v1.x - p1.x  ;  y11 = v1.y - p1.y
    x00 = v0.x - p0.x  ;  y00 = v0.y - p0.y
    x10 = v1.x - p0.x  ;  y10 = v1.y - p0.y

    if (x01*x11 + y01*y11)*(x10*y00 - x00*y10) < (y01*x11 - x01*y11)*(x10*x00 + y00*y10)
        # related adjacent cells
        a = cell1.adjacs[inc(pos_p1,2)]
        b = cell1.adjacs[pos_p1]
        c = cell0.adjacs[inc(pos_p0,2)]
        d = cell0.adjacs[pos_p0]

        cell0.points = [p0, v1, p1]
        cell1.points = [p0, p1, v0]
        cell0.adjacs = [d, a, cell1]
        cell1.adjacs = [cell0, b, c]

        if a != nothing; a.adjacs[ findfirst(a.points, p1) ] = cell0 end
        if c != nothing; c.adjacs[ findfirst(c.points, p0) ] = cell1 end

        return true
    end

    return false
end

# Function to remove the adjacency of a cell
function remove_adjac(cell::TCell, acell::TCell)
    for (i,p) in enumerate(cell.points)
        if !(p in acell.points)
            cell.adjacs[inc(i,1)] = nothing
        end
    end
end

# Function that triangulates an array of points
function triangulate(coords::Array{Float64,2}; getedges=false)
    cells  = TCell[]
    points = TPoint[]
    n = size(coords,1)

    # input limits
    minx, miny = minimum(coords,1)
    maxx, maxy = maximum(coords,1)

    # auxiliar length to mount the super triangle
    Δ = 10.0*max(maxx-minx, maxy-miny, 1.0)

    # points of super triangle
    pp0 = TPoint(minx-Δ  , miny-Δ)
    pp1 = TPoint(maxx+3*Δ  , miny-Δ)
    pp2 = TPoint(minx-Δ  , maxy+3*Δ)

    # addition of super triangle
    push!(cells, TCell(pp0, pp1, pp2) )

    # loop for each input vertex
    for i=1:n
        # input point
        p = TPoint(coords[i,1:end]...)
        push!(points, p)

        # find cell where p falls
        cell = find_cell(p.x, p.y, cells[end])

        # get cell points and adjacent cells
        p0, p1, p2 = cell.points
        a0, a1, a2 = cell.adjacs

        # define new cells
        t0 = cell # replace current cell
        t0.points = [p, p0, p1]
        t1 = TCell(p, p1, p2); push!(cells, t1)
        t2 = TCell(p, p2, p0); push!(cells, t2)

        # new cells adjacency
        t0.adjacs = [t2, a0, t1]
        t1.adjacs = [t0, a1, t2]
        t2.adjacs = [t1, a2, t0]

        # update neighbors in adjacent cells
        if a0 != nothing; a0.adjacs[ findfirst(a0.points, p1) ] = t0 end
        if a1 != nothing; a1.adjacs[ findfirst(a1.points, p2) ] = t1 end
        if a2 != nothing; a2.adjacs[ findfirst(a2.points, p0) ] = t2 end

        # swap common edge in cells if necessary according to Delaunay triangulation
        toswap = [t0, t1, t2]
        while length(toswap)>0
            cell0 = pop!(toswap)
            cell1 = cell0.adjacs[2] # adjacent cell opposito to point P in cell0 cell
            if cell1==nothing
                continue
            end

            if swap_cells(cell0, cell1)
                push!(toswap, cell0)
                push!(toswap, cell1)
            end

        end

    end

    # filter cells and fix adjacency
    spts = [pp0, pp1, pp2] # point in super triangle
    n    = length(cells)
    idxs = Int64[]
    for (i,cell) in enumerate(cells)
        toremove = false
        for p in spts
            if p in cell.points
                toremove = true
                break
            end
        end

        if toremove
            for adj in cell.adjacs
                if adj == nothing; continue end
                remove_adjac(adj, cell)
            end
        else
            push!(idxs, i)
        end
    end
    cells = cells[idxs] # filter cells

    # numbering points
    for (i,point) in enumerate(points)
        point.id = i
    end

    # points coordinates matrix
    Co = Array{Float64}(undef, length(points), 2)
    for (i,point) in enumerate(points)
        Co[i, 1] = point.x
        Co[i, 2] = point.y
    end

    # generate cell connectivities
    n = length(cells)
    Con = Array{Int64}(n, 3)
    for (i,cell) in enumerate(cells)
        Con[i,1] = cell.points[1].id
        Con[i,2] = cell.points[2].id
        Con[i,3] = cell.points[3].id
    end
    Con = [ cell.points[j].id for cell in cells, j=1:3 ]

    # generate edge connectivities
    edges = Set{Tuple{Int64, Int64}}()
    lcon  = ((1,2), (2,3), (3,1)) # local connectivities
    for (i,cell) in enumerate(cells)
        for (pos1, pos2) in lcon
            id1, id2 = cell.points[pos1].id, cell.points[pos2].id
            push!(edges, (min(id1,id2), max(id1, id2)))
        end
    end
    E = [ edge[j] for edge in edges, j=1:2 ]

    if getedges
        return E
    else
        return Con
    end
end


if true
    #V = [ 1.0 2; 2 1; 1 3 ; 2 4; 3 2; 4 3 ]
    #V = [ 0.0 0; 1 0; 0 1 ]
    V = [ rand() for i=1:100, j=1:2]
    E = triangulate(V, getedges=true)

    using FemMesh
    bl = BlockTruss(V, E)
    mesh = Mesh(bl)
    save(mesh, "tri1.vtk")
end


