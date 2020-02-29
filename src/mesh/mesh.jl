# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

#include("delaunay.jl")

### Type Mesh
"""
    Mesh()

Constructs an unitialized mesh object to be used in finite element analyses.
It contains geometric fields as: points, cells, faces, edges, ndim, quality, etc.
"""
mutable struct Mesh
    ndim   ::Int
    points ::Array{Point,1}
    cells  ::Array{Cell,1}

    faces  ::Array{Cell,1}
    edges  ::Array{Cell,1}

    pointdict::Dict{UInt64,Point}
    cellpartition::CellPartition

    #quality::Float64
    #qmin   ::Float64

    # Data
    point_scalar_data::OrderedDict{String,Array}
    cell_scalar_data ::OrderedDict{String,Array}
    point_vector_data::OrderedDict{String,Array}

    function Mesh()
        this = new()
        this.points = []
        this.pointdict = Dict{UInt64, Point}()
        this.cells  = []
        this.faces  = []
        this.edges  = []
        this.ndim   = 0
        this.cellpartition = CellPartition()
        #this.quality = 0.0
        #this.qmin    = 0.0

        this.point_scalar_data = OrderedDict()
        this.cell_scalar_data  = OrderedDict()
        this.point_vector_data = OrderedDict()
        return this
    end
end


function Base.copy(mesh::Mesh)
    newmesh = Mesh()
    newmesh.ndim = mesh.ndim
    newmesh.points = copy.(mesh.points)
    newmesh.cells  = copy.(mesh.cells)
    for c in newmesh.cells
        ids = [ p.id for p in c.points ]
        c.points = newmesh.points[ids]
    end

    newmesh.pointdict = Dict( k => newmesh.points[p.id] for (k,p) in mesh.pointdict )
    fixup!(newmesh)
    return newmesh
end

function datafields(mesh::Mesh)
    return [
            collect(keys(mesh.point_scalar_data));
            collect(keys(mesh.point_vector_data));
            collect(keys(mesh.cell_scalar_data));
           ]
end


function get_surface(cells::Array{Cell,1})::Array{Cell,1}
    surf_dict = Dict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_faces(cell)
            #conns = [ pt.id for pt in face.points ]
            #hs = hash(conns)
            hs = hash(face)
            if haskey(surf_dict, hs)
                delete!(surf_dict, hs)
            else
                surf_dict[hs] = face
            end
        end
    end

    return [ face for face in values(surf_dict) ]
end

function get_edges(surf_cells::Array{Cell,1})::Array{Cell,1}
    edges_dict = Dict{UInt64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in get_faces(cell)
            edge.ocell = cell.ocell
            edges_dict[hash(edge)] = edge
        end
    end

    return [ edge for edge in values(edges_dict) ]
end

# Return a list of neighbors for each cell
function get_neighbors(cells::Array{Cell,1})::Array{Cell,1}
    faces_dict = Dict{UInt64, Cell}()
    neighbors = [ Cell[] for i=1:length(cells) ]

    # Get cell faces. If dup, original and dup are deleted but neigh info saved
    for cell in cells
        for face in get_faces(cell)
            hs = hash(face)
            other = get(faces_dict, hs, nothing)
            if other == nothing
                faces_dict[hs] = face
            else
                push!(neighbors[face.ocell.id], other.ocell)
                push!(neighbors[other.ocell.id], face.ocell)
                delete!(faces_dict, hs)
            end
        end
    end

    return neighbors

end

# Return the cell patch for each point
function get_patches(mesh::Mesh)
    patches = [ Cell[] for i=1:length(mesh.points) ]
    for cell in mesh.cells
        for pt in cell.points
            push!(patches[pt.id], cell)
        end
    end

    return mesh.points, patches
end

# Reverse Cuthill–McKee algorithm (RCM)
function reorder!(mesh::Mesh; sort_degrees=true, reversed=false)

    # Get all mesh edges
    all_edges = Dict{UInt64, Cell}()
    for cell in mesh.cells

        # adding cell edges
        if cell.shape.family == SOLID_SHAPE #is_solid(cell.shape)
            for edge in get_edges(cell)
                hs = hash(edge)
                all_edges[hs] = edge
            end

            #check for lagrangian elements
            if cell.shape==QUAD9
                edge = Cell(POLYV, [ cell.points[1], cell.points[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
            end
            continue
        end

        # joint1D cells (semi-embedded approach)
        if cell.shape in (JLINK2, JLINK3)
            edge = Cell(POLYV, [ cell.points[1], cell.points[end-1] ])
            hs = hash(edge)
            all_edges[hs] = edge
            continue
        end

        # embedded line cells
        if cell.shape.family == LINE_SHAPE && length(cell.linked_cells)>0
            edge1 = Cell(cell.shape, cell.points)
            edge2 = Cell(LIN2, [ cell.points[1], cell.linked_cells[1].points[1] ])
            all_edges[hash(edge1)] = edge1
            all_edges[hash(edge2)] = edge2
            continue
        end

        # all other cells
        edge = Cell(cell.shape, cell.points)
        hs = hash(edge)
        all_edges[hs] = edge

    end

    # Get neighbors ids
    npoints = length(mesh.points)
    neighs_ids = Array{Int64}[ [] for i in 1:npoints ]

    for edge in values(all_edges)
        points = edge.points
        np = length(points)
        for i=1:np-1
            for j=i+1:np
                push!(neighs_ids[points[i].id], points[j].id)
                push!(neighs_ids[points[j].id], points[i].id)
            end
        end
    end

    # removing duplicates and get neighbors
    neighs = Array{Point}[ mesh.points[unique(list)] for list in neighs_ids ]
    #neighs = Array{Point}[ mesh.points[unique(list)] for list in neighs_ids ]

    # list of degrees per point
    degrees = Int64[ length(list) for list in neighs]
    mindeg, idx  = findmin(degrees)

    if mindeg == 0
        # Case of overlapping elements where edges have at least one point with the same coordinates
        @warn "reorder!: Reordering nodes failed! Check for overlapping cells or non used points."
        return
    end

    N = [ mesh.points[idx] ] # new list of points
    L = Dict{Int64,Point}()  # last levelset. Use ids as keys instead of hash to avoid collisions of points with same coordinates
    L[idx] = mesh.points[idx]
    LL = Dict{Int64,Point}()  # levelset before the last one

    while length(N) < npoints
        # Generating current levelset A
        A = Dict{Int64,Point}()

        for p in values(L)
            for q in neighs[p.id]
                (haskey(L, q.id) || haskey(LL, q.id)) && continue
                A[q.id] = q
            end
        end
        if length(A)==0
            @error "reorder!: Reordering nodes failed! Possible error with cell connectivities."
            return
        end

        # Convert A into an array RA
        RA = collect(values(A))
        if sort_degrees
            D  = [ degrees[point.id] for point in RA ]
            RA = RA[sortperm(D)]
        end

        append!(N, RA)
        LL, L = L, A
    end

    # Reverse list of new nodes
    if reversed
        N = reverse(N)
    end

    # Renumerating
    for (i,p) in enumerate(N)
        p.id = i
    end

    mesh.points = N

    return nothing

end


function renumber!(mesh::Mesh)
    # Get ndim
    ndim = 1
    for point in mesh.points
        point.y != 0.0 && (ndim=2)
        point.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim

    # Numberig nodes
    for (i,p) in enumerate(mesh.points) p.id = i end

    # Numberig cells
    for (i,c) in enumerate(mesh.cells )
        c.id = i;
        c.ndim=ndim;
    end

    mesh.point_scalar_data["point-id"] = collect(1:length(mesh.points))
    mesh.cell_scalar_data["cell-id"]   = collect(1:length(mesh.cells))
    mesh.cell_scalar_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.cells ]
end

function set_faces_edges!(mesh::Mesh)
    # Facets
    if genfacets
        verbose && print("  finding facets...   \r")
        mesh.faces = get_surface(mesh.cells)
    end
    ndim==2 && (mesh.edges=mesh.faces)

    # Edges
    if genedges && ndim==3
        verbose && print("  finding edges...   \r")
        mesh.edges = get_edges(mesh.faces)
    end
end

function update_quality!(mesh::Mesh)
    # Quality
    Q = Float64[]
    for c in mesh.cells
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end
    mesh.cell_scalar_data["quality"]   = Q
end

# Updates numbering, faces and edges in a Mesh object
function fixup!(mesh::Mesh; verbose::Bool=false, genfacets::Bool=true, genedges::Bool=true, reorder::Bool=false)

    # Get ndim
    ndim = 1
    for point in mesh.points
        point.y != 0.0 && (ndim=2)
        point.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim

    # Numberig nodes
    for (i,p) in enumerate(mesh.points) p.id = i end

    # Numberig cells
    for (i,c) in enumerate(mesh.cells )
        c.id = i;
        c.ndim=ndim;
    end

    # Facets
    if genfacets
        verbose && print("  finding facets...   \r")
        mesh.faces = get_surface(mesh.cells)
    end
    ndim==2 && (mesh.edges=mesh.faces)

    # Edges
    if genedges && ndim==3
        verbose && print("  finding edges...   \r")
        mesh.edges = get_edges(mesh.faces)
    end

    # Quality
    Q = Float64[]
    for c in mesh.cells
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end

    # Ordering
    reorder && reorder!(mesh)

    # Reset data
    #mesh.point_scalar_data = OrderedDict()
    #mesh.cell_scalar_data  = OrderedDict()
    #mesh.point_vector_data = OrderedDict()
    mesh.point_scalar_data["point-id"] = collect(1:length(mesh.points))
    mesh.cell_scalar_data["quality"]   = Q
    mesh.cell_scalar_data["cell-id"]   = collect(1:length(mesh.cells))
    mesh.cell_scalar_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.cells ]

    return nothing
end

# Mesh quality
function quality!(mesh::Mesh)
    for c in mesh.cells
        c.quality = cell_quality(c)
    end
    Q = Float64[ c.quality for c in mesh.cells]
    #mesh.quality = sum(Q)/length(mesh.cells)
    #mesh.qmin    = minimum(Q)
    mesh.cell_scalar_data["quality"] = Q
    return nothing
end

# Adds m2 to m1
function join_mesh!(m1::Mesh, m2::Mesh)
    #m = Mesh()

    # copy m1 to m
    #m.points = copy(m1.points)
    #m.cells  = copy(m1.cells)
    #m.pointdict = copy(m1.pointdict)

    pid = length(m1.points)
    cid = length(m1.cells)

    #@show length(m1.points)
    #@show length(m2.points)

    #for (h,p) in m1.pointdict
        #@show (p.id, p.x, p.y, p.z)
    #end

    # Add new points from m2
    for p in m2.points
        hs = hash(p)
        if !haskey(m1.pointdict, hs)
            #@show (p.x, p.y, p.z)
            pid += 1
            p.id = pid
            push!(m1.points, p)
        end
    end

    #@show length(m1.points)

    # Fix m2 cells connectivities for points at m1 border
    for c in m2.cells
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            if haskey(m1.pointdict, hs)
                pp = m1.pointdict[hs]
                c.points[i] = pp
            else
                # update pointdict dict
                if haskey(m2.pointdict, hs)
                    m1.pointdict[hs] = p
                end
            end
        end

        cid += 1
        c.id = cid
        push!(m1.cells, c)
    end

    return nothing
end


"""
    Mesh(coords, conns, cellshapes=[])

# Arguments

`coords` : Matrix with point coordinates

`conns`  : Array of connectivities for each cell

`cellshapes=[]` : Array of cell shapes
"""
function Mesh(
              coords     :: Array{<:Real},
              conns      :: Array{Array{Int64,1},1},
              cellshapes :: Array{ShapeType,1}=ShapeType[];
              tag        :: String=""
             )
    n = size(coords, 1) # number of points
    m = size(conns , 1) # number of cells

    points = Point[]
    for i=1:n
        C = coords[i,:]
        push!(points, Point(C))
    end

    # Get ndim
    ndim = size(coords,2)

    cells = Cell[]
    for i=1:m
        pts = points[conns[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        cell = Cell(shape, pts, tag=tag)
        push!(cells, cell)
    end

    mesh = Mesh()
    mesh.points = points
    mesh.cells  = cells
    fixup!(mesh, reorder=false) # no node ordering
    mesh.ndim = ndim # force ndim

    return mesh
end


# flattens a list of nested lists
function flatten(x, y)
    ty = typeof(x)
    if ty <: Tuple || ty <: Array
        for item in x
            flatten(item, y)
        end
    else
        push!(y, x)
    end
    return y
end
flatten(x)=flatten(x, [])


"""
    Mesh(items, kwargs...)

Generates a mesh based on an array of geometrical objects.

# Arguments

`items`     : Array of objects used to generate a `Mesh` object.
These objects can be of type `Block` or `Mesh`.
Subarrays of these type of objects are also supported.

# Keyword arguments

`genfacets = true` : If true, generates facet cells

`genedges  = true` : If true, generates edge cells

`reorder   = true` : If true, reorder nodes numbering

`verbose   = false` : If true, prints extra information

`silent    = false` : If true, does not print anything
"""
function Mesh(
              items     :: Union{Mesh, AbstractBlock, Array{<:Union{AbstractBlock, Array},1}}...;
              genfacets :: Bool = true,
              genedges  :: Bool = true,
              reorder   :: Bool = true,
              verbose   :: Bool = false,
              silent    :: Bool = false,
             )

    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    # Flatten items list
    fitems = flatten(items)

    # Get list of blocks and update mesh objects
    blocks = []
    meshes = []
    for item in fitems
        if isa(item, AbstractBlock)
            push!(blocks, item)
        elseif isa(item, Mesh)
            push!(meshes, copy(item))
        else
            error("Mesh: item of type $(typeof(item)) cannot be used as Block or Mesh")
        end
    end

    nmeshes = length(meshes)
    nblocks = length(blocks)
    if verbose
        verbosity>0 && printstyled("Mesh generation:\n", bold=true, color=:cyan)
        nmeshes>0 && @printf "  %5d meshes\n" nmeshes
        @printf "  %5d blocks\n" nblocks
    end

    # New mesh object
    mesh = Mesh()

    # Join meshes
    for m in meshes
        join_mesh!(mesh, m)
    end

    # Split blocks: generates points and cells
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        verbosity>1 && print("  spliting block ", i, "...    \r")
    end

    # Updates numbering, quality, facets and edges
    fixup!(mesh, verbose=verbose, genfacets=genfacets, genedges=genedges, reorder=reorder)

    if verbosity>0
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" npoints
        @printf "  %5d cells\n" ncells
    end
    if verbosity>1
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        if genfacets
            @printf "  %5d faces\n" nfaces
        end
        if genedges
            @printf "  %5d surface edges\n" nedges
        end
        icount = 0
        for block in blocks
            if typeof(block) == BlockInset; icount += block.icount end
        end
        if icount>0
            @printf "  %5d intersections\n" icount
        end
    end

    return mesh
end


function Mesh(cells::Array{Cell,1})
    # New mesh object TODO: make a copy?
    mesh = Mesh()
    mesh.points = cells[:points]
    mesh.cells  = cells
    fixup!(mesh, reorder=false)

    return mesh
end


# Gets a part of a mesh filtering elements
function Base.getindex(mesh::Mesh, filter_ex::Expr)
    # TODO: make a copy?

    # filter cells
    cells  = mesh.cells[filter_ex]
    # get points
    points = get_points(cells)

    # ids from selected cells and poitns
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in points ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.points = points
    new_mesh.cells = cells

    # select relevant data
    for (key,vals) in mesh.point_scalar_data
        new_mesh.point_scalar_data[key] = vals[pids]
    end

    for (key,vals) in mesh.cell_scalar_data
        new_mesh.cell_scalar_data[key] = vals[cids]
    end

    for (key,vals) in mesh.point_vector_data
        new_mesh.point_vector_data[key] = vals[pids,:]
    end

    # update node numbering, facets and edges
    fixup!(new_mesh, reorder=false)

    return new_mesh

end


function threshold(mesh::Mesh, field::Union{Symbol,String}, minval::Float64, maxval::Float64)
    field = string(field)

    # check if field exists
    found = haskey(mesh.cell_scalar_data, field)
    if found
        vals = mesh.cell_scalar_data[field]
    else
        found = haskey(mesh.point_scalar_data, field)
        found || error("threshold: field $field not found")
        data  = mesh.point_scalar_data[field]
        vals = [ mean(data[i]) for i=1:ncells ]
    end

    # filter cells
    cells = Cell[]
    for (cell, val) in zip(mesh.cells, vals)
        if minval <= val <= maxval
            push!(cells, cell)
        end
    end

    # get points
    points = get_points(cells)

    # ids from selected cells and points
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in points ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.points = points
    new_mesh.cells = cells

    # select relevant data
    for (key,vals) in mesh.point_scalar_data
        new_mesh.point_scalar_data[key] = vals[pids]
    end

    for (key,vals) in mesh.cell_scalar_data
        new_mesh.cell_scalar_data[key] = vals[cids]
    end

    for (key,vals) in mesh.point_vector_data
        new_mesh.point_vector_data[key] = vals[pids,:]
    end

    # update node numbering, facets and edges
    fixup!(new_mesh, reorder=false)

    return new_mesh

end


export get_segment_data
function get_segment_data(msh::Mesh, X1::Array{<:Real,1}, X2::Array{<:Real,1}, filename::String=""; npoints=200)
    data = msh.point_scalar_data
    table = DTable(["s"; collect(keys(data))])
    X1 = [X1; 0.0][1:3]
    X2 = [X2; 0.0][1:3]
    Δ = (X2-X1)/(npoints-1)
    Δs = norm(Δ)
    s1 = 0.0

    for i=1:npoints
        X = X1 + Δ*(i-1)
        s = s1 + Δs*(i-1)
        cell = find_cell(X, msh.cells, msh.cellpartition, 1e-7, Cell[])
        coords = getcoords(cell)
        R = inverse_map(cell.shape, coords, X)
        N = cell.shape.func(R)
        map = [ p.id for p in cell.points ]
        vals = [ s ]
        for (k,V) in data
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
    end

    filename != "" && save(table, filename)

    return table
end


export randmesh
function randmesh(l::Real...)
    ndim = length(l)
    if ndim==2
        lx, ly = l
        nx, ny = rand(4:7, 2)
        cellshape = rand((TRI3, TRI6, QUAD4, QUAD8))
        m = Mesh(Block2D([0.0 0.0; lx ly], nx=nx, ny=ny, cellshape=cellshape), verbose=false)
    else
        lx, ly, lz = l
        nx, ny, nz = rand(4:7, 3)
        cellshape = rand((TET4, TET10, HEX8, HEX20))
        m = Mesh(Block3D([0.0 0.0 0.0; lx ly lz], nx=nx, ny=ny, nz=nz, cellshape=cellshape), verbose=false)
    end
end

function check_mesh(mesh::Mesh)
    cell_dict = Dict{UInt64, Cell}()

    # Get only unique elems. If dup, original and dup are deleted
    n = 0
    for cell in mesh.cells
        hs = hash(cell)
        if haskey(cell_dict, hs)
            n += 1
        else
            cell_dict[hs] = cell
        end
    end

    return n

end


