# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    Mesh()

Constructs an unitialized mesh object to be used in finite element analyses.
It contains geometric fields as: nodes, elems, faces, edges, ndim, quality, etc.
"""
mutable struct Mesh<:AbstractMesh
    ndim ::Int
    nodes::Array{Node,1}
    elems::Array{Cell,1}

    faces::Array{Cell,1}
    edges::Array{Cell,1}

    _pointdict::Dict{UInt64,Node}
    _elempartition::ElemPartition

    # Data
    node_data::OrderedDict{String,Array}
    elem_data ::OrderedDict{String,Array}

    function Mesh()
        this = new()
        this.nodes = []
        this._pointdict = Dict{UInt64, Node}()
        this.elems  = []
        this.faces  = []
        this.edges  = []
        this.ndim   = 0
        this._elempartition = ElemPartition()
        this.node_data = OrderedDict()
        this.elem_data = OrderedDict()
        return this
    end
end

function get_node(nodes::Dict{UInt64,Node}, C::AbstractArray{<:Real})
    hs = hash(Node(C))
    return get(nodes, hs, nothing)
end

function Base.copy(mesh::Mesh)
    newmesh = Mesh()
    newmesh.ndim = mesh.ndim
    newmesh.nodes = copy.(mesh.nodes)

    for elem in mesh.elems
        idxs = [ node.id for node in elem.nodes ]
        newelemnodes = newmesh.nodes[idxs]
        newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
        push!(newmesh.elems, newelem)
    end

    newmesh._pointdict = Dict( hash(node) => node for node in newmesh.nodes )

    fixup!(newmesh)
    return newmesh
end

function datafields(mesh::Mesh)
    return [
            collect(keys(mesh.node_data));
            collect(keys(mesh.elem_data));
           ]
end


function get_surface(cells::Array{<:AbstractCell,1})
    surf_dict = Dict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_faces(cell)
            #conns = [ pt.id for pt in face.nodes ]
            #hs = hash(conns)
            hs = hash(face)
            if haskey(surf_dict, hs)
                delete!(surf_dict, hs)
            else
                surf_dict[hs] = face
            end
        end
    end

    return Face[ face for face in values(surf_dict) ]
end

function get_edges(surf_cells::Array{<:AbstractCell,1})
    edges_dict = Dict{UInt64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in get_faces(cell)
            edge.oelem = cell.oelem
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
                push!(neighbors[face.oelem.id], other.oelem)
                push!(neighbors[other.oelem.id], face.oelem)
                delete!(faces_dict, hs)
            end
        end
    end

    return neighbors

end

# Return the cell patch for each point
function get_patches(mesh::Mesh)
    patches = [ Cell[] for i=1:length(mesh.nodes) ]
    for cell in mesh.elems
        for pt in cell.nodes
            push!(patches[pt.id], cell)
        end
    end

    return mesh.nodes, patches
end

# Reverse Cuthill–McKee algorithm (RCM)
function reorder!(mesh::Mesh; sort_degrees=true, reversed=false)

    # Get all mesh edges
    all_edges = Dict{UInt64, Cell}()
    for cell in mesh.elems

        # adding cell edges
        if cell.shape.family == SOLID_SHAPE #is_solid(cell.shape)
            for edge in get_edges(cell)
                hs = hash(edge)
                all_edges[hs] = edge
            end

            #check for lagrangian elements
            if cell.shape==QUAD9
                edge = Cell(POLYV, [ cell.nodes[1], cell.nodes[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
            end
            if cell.shape==HEX27
                edge = Cell(POLYV, [ cell.nodes[1], cell.nodes[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
                for i in (21,22,23,24,25,26)
                    edge = Cell(POLYV, [ cell.nodes[i], cell.nodes[end] ] )
                    hs   = hash(edge)
                    all_edges[hs] = edge
                end
            end
            continue
        end

        # joint1D cells (semi-embedded approach)
        if cell.shape in (JLINK2, JLINK3)
            npts = cell.shape.npoints
            edge = Cell(POLYV, [ cell.nodes[1]; cell.nodes[end-npts+1:end] ])
            hs = hash(edge)
            all_edges[hs] = edge
            continue
        end

        # embedded line cells
        if cell.shape.family == LINE_SHAPE && length(cell.linked_elems)>0
            edge1 = Cell(cell.shape, cell.nodes)
            edge2 = Cell(LIN2, [ cell.nodes[1], cell.linked_elems[1].nodes[1] ])
            all_edges[hash(edge1)] = edge1
            all_edges[hash(edge2)] = edge2
            continue
        end

        # all other cells
        edge = Cell(cell.shape, cell.nodes)
        hs = hash(edge)
        all_edges[hs] = edge

    end

    # Get neighbors ids
    npoints = length(mesh.nodes)
    neighs_ids = Array{Int64}[ [] for i in 1:npoints ]

    for edge in values(all_edges)
        nodes = edge.nodes
        np = length(nodes)
        for i=1:np-1
            for j=i+1:np
                push!(neighs_ids[nodes[i].id], nodes[j].id)
                push!(neighs_ids[nodes[j].id], nodes[i].id)
            end
        end
    end

    # removing duplicates and get neighbors
    neighs = Array{Node}[ mesh.nodes[unique(list)] for list in neighs_ids ]
    #neighs = Array{Node}[ mesh.nodes[unique(list)] for list in neighs_ids ]

    # list of degrees per point
    degrees = Int64[ length(list) for list in neighs]
    mindeg, idx  = findmin(degrees)

    if mindeg == 0
        # Case of overlapping elements where edges have at least one point with the same coordinates
        @warn "reorder!: Reordering nodes failed! Check for overlapping cells or non used nodes."
        return
    end

    N = [ mesh.nodes[idx] ] # new list of nodes
    L = Dict{Int64,Node}()  # last levelset. Use ids as keys instead of hash to avoid collisions of nodes with same coordinates
    L[idx] = mesh.nodes[idx]
    LL = Dict{Int64,Node}()  # levelset before the last one

    while length(N) < npoints
        # Generating current levelset A
        A = Dict{Int64,Node}()

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

    mesh.nodes = N

    return nothing

end


function update_quality!(mesh::Mesh)
    # Quality
    Q = Float64[]
    for c in mesh.elems
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end
    mesh.elem_data["quality"]   = Q
end


# Updates numbering, faces and edges in a Mesh object
function fixup!(mesh::Mesh; verbose::Bool=false, genfacets::Bool=true, genedges::Bool=true, reorder::Bool=false)

    # Get ndim
    ndim = 1
    for point in mesh.nodes
        point.coord.y != 0.0 && (ndim=2)
        point.coord.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim

    # Numberig nodes
    for (i,p) in enumerate(mesh.nodes) 
        p.id = i 
    end

    # Numberig cells
    for (i,c) in enumerate(mesh.elems )
        c.id = i;
    end

    # Facets
    if genfacets
        verbose && print("  finding facets...   \r")
        mesh.faces = get_surface(mesh.elems)
    end
    ndim==2 && (mesh.edges=mesh.faces)

    # Edges
    if genedges && ndim==3
        verbose && print("  finding edges...   \r")
        mesh.edges = get_edges(mesh.faces)
    end

    # Quality
    Q = Float64[]
    for c in mesh.elems
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end

    # Ordering
    reorder && reorder!(mesh)

    # Reset data
    mesh.node_data["node-id"] = collect(1:length(mesh.nodes))
    mesh.elem_data["quality"]   = Q
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))
    mesh.elem_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.elems ]

    return nothing
end

# Mesh quality
function quality!(mesh::Mesh)
    for c in mesh.elems
        c.quality = cell_quality(c)
    end
    Q = Float64[ c.quality for c in mesh.elems]
    mesh.elem_data["quality"] = Q
    return nothing
end


function join_mesh!(mesh::Mesh, m2::Mesh)

    for node in m2.nodes
        push!(mesh.nodes, node)
        mesh._pointdict[hash(node)] = node
    end

    for elem in m2.elems
        newelemnodes = Node[ mesh._pointdict[hash(node)] for node in elem.nodes ]
        newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
        push!(mesh.elems, newelem)
    end

    fixup!(mesh, reorder=false)

    return nothing
end


function join_meshes(m1::Mesh, m2::Mesh)
    mesh = copy(m1)
    join_mesh!(mesh, m2)
    return mesh
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
              tag        :: String="",
              verbose    :: Bool=true,
              silent     :: Bool=false
             )
    n = size(coords, 1) # number of nodes
    m = size(conns , 1) # number of cells

    nodes = Node[]
    for i=1:n
        C = coords[i,:]
        push!(nodes, Node(C))
    end

    # Get ndim
    ndim = size(coords,2)

    cells = Cell[]
    for i=1:m
        pts = nodes[conns[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        cell = Cell(shape, pts, tag=tag)
        push!(cells, cell)
    end

    mesh = Mesh()
    mesh.nodes = nodes
    mesh.elems  = cells
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
    if verbosity>0
        printstyled("Mesh generation:\n", bold=true, color=:cyan)
        nmeshes>0 && @printf "  %5d meshes\n" nmeshes
        @printf "  %5d blocks\n" nblocks
    end

    # New mesh object
    mesh = Mesh()

    # Join meshes
    for m in meshes
        join_mesh!(mesh, m)
    end

    # Split blocks: generates nodes and cells
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        verbosity>1 && print("  spliting block ", i, "...    \r")
    end

    # Updates numbering, quality, facets and edges
    fixup!(mesh, verbose=verbose, genfacets=genfacets, genedges=genedges, reorder=reorder)

    # Add field for embedded nodes
    if any( c.shape.family==JOINT1D_SHAPE for c in mesh.elems )
        ncells = length(mesh.elems)
        inset_data = zeros(Int, ncells, 3) # npoints, first link id, second link id
        for i=1:ncells
            cell = mesh.elems[i]
            if cell.shape.family==JOINT1D_SHAPE
                inset_data[i,1] = cell.shape.npoints
                inset_data[i,2] = cell.linked_elems[1].id
                inset_data[i,3] = cell.linked_elems[2].id
            end
        end
        mesh.elem_data["inset-data"] = inset_data
    end

    if verbosity>0
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d nodes\n" npoints
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


function Mesh(elems::Array{Cell,1})
    newmesh = Mesh()
    newmesh.nodes = copy.(get_nodes(elems))
    newmesh._pointdict = Dict{UInt64,Node}(hash(node)=>node for node in newmesh.nodes)

    for elem in elems
        newelemnodes = Node[ newmesh._pointdict[hash(node)] for node in elem.nodes ]
        newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
        push!(newmesh.elems, newelem)
    end

    fixup!(newmesh, reorder=false)

    return newmesh
end


function Mesh(mesh::Mesh, filter::Union{String,Expr,Symbol})
    elems = mesh.elems[filter]
    elemids = [ elem.id for elem in elems ]

    nodedict = OrderedDict{Int,Node}()
    for elem in elems
        for node in elem.nodes
            nodedict[node.id] = copy(node)
        end
    end
    nodeids  = collect(keys(nodedict))

    newmesh = Mesh()
    newmesh.nodes = collect(values(nodedict))

    for elem in elems
        newelemnodes = Node[ nodedict[node.id] for node in elem.nodes ]
        newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
        push!(newmesh.elems, newelem)
    end

    for (k,v) in mesh.elem_data
        newmesh.elem_data[k] = v[elemids]
    end

    for (k,v) in mesh.node_data
        newmesh.node_data[k] = v[nodeids]
    end

    fixup!(newmesh, reorder=false)

    return newmesh
end


function slice(mesh::Mesh, normal::Array{Float64,1})
    mesh = Mesh()
    return mesh
end


function clip(mesh::Mesh, normal::Array{Float64,1})
    mesh = Mesh()
    return mesh
end


function threshold(mesh::Mesh, field::Union{Symbol,String}, minval::Float64, maxval::Float64)
    field = string(field)

    # check if field exists
    found = haskey(mesh.elem_data, field)
    if found
        vals = mesh.elem_data[field]
    else
        found = haskey(mesh.node_data, field)
        found || error("threshold: field $field not found")
        data  = mesh.node_data[field]
        vals = [ mean(data[i]) for i=1:ncells ]
    end

    # filter cells
    cells = Cell[]
    for (cell, val) in zip(mesh.elems, vals)
        if minval <= val <= maxval
            push!(cells, cell)
        end
    end

    # get nodes
    nodes = get_nodes(cells)

    # ids from selected cells and nodes
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in nodes ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.nodes = nodes
    new_mesh.elems = cells

    # select relevant data
    for (key,vals) in mesh.node_data
        new_mesh.node_data[key] = vals[pids,:]
    end

    for (key,vals) in mesh.elem_data
        new_mesh.elem_data[key] = vals[cids]
    end

    # update node numbering, facets and edges
    fixup!(new_mesh, reorder=false)

    return new_mesh

end


export get_segment_data
function get_segment_data(msh::AbstractMesh, X1::Array{<:Real,1}, X2::Array{<:Real,1}, filename::String=""; n=200)
    data = msh.node_data
    table = DataTable(["s"; collect(keys(data))])
    X1 = [X1; 0.0][1:3]
    X2 = [X2; 0.0][1:3]
    Δ = (X2-X1)/(n-1)
    Δs = norm(Δ)
    s1 = 0.0

    for i=1:n
        X = X1 + Δ*(i-1)
        s = s1 + Δs*(i-1)
        cell = find_elem(X, msh.elems, msh._elempartition, 1e-7, Cell[])
        coords =get_coords(cell)
        R = inverse_map(cell.shape, coords, X)
        N = cell.shape.func(R)
        map = [ p.id for p in cell.nodes ]
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



export UMesh

function UMesh(
               polym::Union{Polygon,PolygonMesh},
               sizes::Array=[],
               embedded=Array=[];
               size::Real=NaN
             )
    if polym isa Polygon
        polym = PolygonMesh([polym])
    end

    coords = get_coords(polym)
    ndim = sum(abs, coords[:,end])==0.0 ? 2 : 3
    npoints = Base.size(coords,1)
    pdict = Dict{Point, Int}( p=>i for (i,p) in enumerate(polym.points) )
    ldict = Dict{UInt, Int}( hash(l)=>i for (i,l) in enumerate(polym.lines) )

    if isnan(size)
        bb = bounding_box(polym.points)
        size = maximum(diff(bb, dims=1))/3
    end

    # Set sizes
    _sizes = zeros(npoints)
    _sizes .= size
    for (filter,s) in sizes
        points = polym.points[filter]
        for p in points
            _sizes[pdict[p]] = s
        end
    end

    # Set embedded points in surfaces
    _embed = Dict{Int,Int}()
    polydict = Dict{UInt, Int}( hash(poly)=>i for (i,poly) in enumerate(polym.polygons) )
    for (filter,X) in embedded
        polygons = polym.polygons[filter]
        length(polygons)==0 && break
        poly = polygons[1]
        coords = [coords; X']
        push!(_sizes, size)
        @show coords
        #_embed[polydict[hash(poly)]] = Base.size(coords,2)
        _embed[ Base.size(coords,1) ] = polydict[hash(poly)]
    end

    # Find point indexes for lines
    pdict = Dict{Point, Int}( p=>i for (i,p) in enumerate(polym.points) )
    lineindexes = Array{Int,1}[]
    for line in polym.lines
        idx = Int[ pdict[p] for p in line.points ]
        push!(lineindexes, idx)
    end

    # Find line indexes for polygons
    polyindexes = Array{Int,1}[]
    for poly in polym.polygons
        idx = Int[ ldict[hash(l)] for l in poly.lines]
        push!(polyindexes, idx)
    end

    # Fix signal in line indexes for polygons
    for loop in polyindexes
        if !(lineindexes[loop[1]][end] in lineindexes[loop[2]])
            loop[1] = -loop[1]
        end

        for i in 2:length(loop)
            lineidx = loop[i]
            line = lineidx>0 ? lineindexes[lineidx] : reverse(lineindexes[abs(lineidx)])
            lastlineidx = loop[i-1]
            lastline = lastlineidx>0 ? lineindexes[lastlineidx] : reverse(lineindexes[abs(lastlineidx)])
            if line[1]!=lastline[end]
                loop[i] = -lineidx
            end
        end
    end

    # Mesh generation
    @eval begin
        coords = $coords
        sizes = $_sizes
        lines = $lineindexes
        surfaces = $polyindexes
        ndim = $ndim
        npoints = $npoints
        embedded = $_embed

        import Gmsh.gmsh

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("model01")

        npoints = size(coords,1)

        # adding points
        for i=1:npoints
            #gmsh.model.geo.addPoint(coords[i,1], coords[i,2], ndim==2 ? 0.0 : coords[i,3], sizes[i], i)
            gmsh.model.geo.addPoint(coords[i,1], coords[i,2], coords[i,3], sizes[i], i)
        end

        # adding lines
        nlines = length(lines)
        for (i,line) in enumerate(lines)
            if length(line)==2 # line
                gmsh.model.geo.addLine(line[1], line[2], i)
            else # circle arc
                gmsh.model.geo.addCircleArc(line[1], line[2], line[3], i)
            end
        end

        #if ndim==2 && length(surfaces)==0
            #surfaces = [ collect(1:nlines) ]
        #end

        # PhysicalGroup index
        ipg = 0

        # adding surfaces
        nsurfs = length(surfaces)
        for (i,loop) in enumerate(surfaces)
            gmsh.model.geo.addCurveLoop(loop, i)
            gmsh.model.geo.addPlaneSurface([i], i)
            if ndim==2
                global ipg += 1
                gmsh.model.addPhysicalGroup(2, [i], ipg) # ndim, entities, tag
            end
        end
        #if _2d
            #gmsh.model.geo.addCurveLoop(collect(1:nlines), 1)
            #gmsh.model.geo.addPlaneSurface([1], 1)
            #gmsh.model.addPhysicalGroup(2, [1], 1) # ndim, entities, tag
        #end

        if ndim==3
            gmsh.model.geo.addSurfaceLoop(collect(1:nsurfs), 1)
            gmsh.model.geo.addVolume([1], 1)
            global ipg += 1
            gmsh.model.addPhysicalGroup(3, [1], ipg) # ndim, entities, tag
        end

        gmsh.model.geo.synchronize()
        for (k,v) in embedded
            gmsh.model.mesh.embed(0,[k],2,v)
        end
        gmsh.model.mesh.generate(ndim)
        gmsh.model.mesh.smooth()

        tempfile = "_temp.vtk"
        gmsh.write(tempfile)
        gmsh.finalize()
    end

    mesh = Mesh(tempfile)
    rm(tempfile)
    return mesh

end
