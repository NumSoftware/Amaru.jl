# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    Mesh

A type that represents a finite element mesh. 
It contains geometric fields as: nodes, elems, faces, edges, ndim, quality, etc.

# Fields
$(FIELDS)
"""
mutable struct Mesh<:AbstractDomain
    " array of nodes"
    nodes::Array{Node,1} 
    " array of elements (cells)"
    elems::Array{Cell,1}
    "array of faces (cells)"
    faces::Array{Cell,1}
    "array of edges (cells)"
    edges::Array{Cell,1}
    "nodal data dictionary"
    node_data::OrderedDict{String,Array}
    "element data dictionary"
    elem_data ::OrderedDict{String,Array}
    "mesh environment"
    env::MeshEnv

    _pointdict::Dict{UInt64,Node}
    _elempartition::ElemPartition

    """
    Mesh()

    Constructs an unitialized mesh object to be used in finite element analyses.
    It contains geometric fields as: nodes, elems, faces, edges, ndim, quality, etc.
    """
    function Mesh(ndim=0)
        this = new()
        this.nodes = []
        this.elems  = []
        this.faces  = []
        this.edges  = []
        # this.ndim   = 0
        this.node_data = OrderedDict()
        this.elem_data = OrderedDict()
        this.env = MeshEnv(ndim)
        this._pointdict = Dict{UInt64, Node}()
        this._elempartition = ElemPartition()
        return this
    end
end


function get_node(nodes::Dict{UInt64,Node}, C::AbstractArray{<:Real})
    hs = hash(Node(C))
    return get(nodes, hs, nothing)
end


function Base.copy(mesh::AbstractDomain)
    ndim = mesh.env.ndim
    newmesh = Mesh(ndim)
    newmesh.nodes = copy.(mesh.nodes)

    for elem in mesh.elems
        idxs = [ node.id for node in elem.nodes ]
        newelemnodes = newmesh.nodes[idxs]
        newelem = Cell(elem.shape, newelemnodes, tag=elem.tag, id=elem.id, active=elem.active)
        push!(newmesh.elems, newelem)
    end

    compute_facets!(newmesh)
    
    newmesh._pointdict = Dict( hash(node) => node for node in newmesh.nodes )
    newmesh.node_data = copy(mesh.node_data)
    newmesh.elem_data = copy(mesh.elem_data)
    
    # fixing references for linked elements #todo check
    for elem in newmesh.elems
        if length(elem.linked_elems)>0
            idxs = [ e.id for e in elem.linked_elems ]
            elem.linked_elems = newmesh.elems[idxs]
        end
    end

    return newmesh
end


function get_outer_facets(cells::Array{<:AbstractCell,1})
    surf_dict = OrderedDict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in getfacets(cell)
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


function getedges(surf_cells::Array{<:AbstractCell,1})
    edges_dict = Dict{UInt64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in getedges(cell)
            edge.owner = cell.owner
            edges_dict[hash(edge)] = edge
        end
    end

    return collect(values(edges_dict))
end

# Return a list of neighbors for each cell
function get_neighbors(cells::Array{Cell,1})::Array{Cell,1}
    faces_dict = Dict{UInt64, Cell}()
    neighbors = [ Cell[] for i in 1:length(cells) ]

    # Get cell faces. If dup, original and dup are deleted but neigh info saved
    for cell in cells
        for face in getfacets(cell)
            hs = hash(face)
            other = get(faces_dict, hs, nothing)
            if other === nothing
                faces_dict[hs] = face
            else
                push!(neighbors[face.owner.id], other.owner)
                push!(neighbors[other.owner.id], face.owner)
                delete!(faces_dict, hs)
            end
        end
    end

    return neighbors

end

# Return the cell patch for each point
function get_patches(mesh::Mesh)
    patches = [ Cell[] for i in 1:length(mesh.nodes) ]
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
        if cell.shape.family == BULKCELL #is_solid(cell.shape)
            for edge in getedges(cell)
                hs = hash(edge)
                all_edges[hs] = edge
            end

            #check for lagrangian elements
            if cell.shape==QUAD9
                edge = Cell(POLYVERTEX, [ cell.nodes[1], cell.nodes[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
            end
            if cell.shape==HEX27
                edge = Cell(POLYVERTEX, [ cell.nodes[1], cell.nodes[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
                for i in (21,22,23,24,25,26)
                    edge = Cell(POLYVERTEX, [ cell.nodes[i], cell.nodes[end] ] )
                    hs   = hash(edge)
                    all_edges[hs] = edge
                end
            end
            continue
        end

        # joint1D cells (semi-embedded approach)
        if cell.shape in (JLINK2, JLINK3)
            npts = cell.shape.npoints
            edge = Cell(POLYVERTEX, [ cell.nodes[1]; cell.nodes[end-npts+1:end] ])
            hs = hash(edge)
            all_edges[hs] = edge
            continue
        end

        # embedded line cells
        if cell.shape.family == LINECELL && length(cell.linked_elems)>0
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
        for i in 1:np-1
            for j in i+1:np
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
        notify("reorder!: Reordering nodes failed. Possible causes: disconnected domain, non used nodes and overlapping cells.")
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
            #@error "reorder!: Reordering nodes failed! Possible error with cell connectivities."
            notify("reorder!: Reordering nodes failed. Possible causes: disconnected domain, non used nodes and overlapping cells.")
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
    mesh.elem_data["quality"] = Q
end


# Updates numbering, faces and edges in a Mesh object
function fixup!(mesh::Mesh; genfacets=true, reorder=false)
    @assert mesh.env.ndim!=0

    # Numberig nodes
    for (i,p) in enumerate(mesh.nodes) 
        p.id = i 
    end

    # Numberig cells and setting env
    for (i,elem) in enumerate(mesh.elems)
        elem.id = i
        elem.env = mesh.env
    end

    # Faces and edges
    if genfacets
        if mesh.env.ndim==2
            mesh.edges = get_outer_facets(mesh.elems)
            mesh.faces = mesh.edges
        elseif mesh.env.ndim==3
            solids = [ elem for elem in mesh.elems if elem.shape.ndim==3 ]
            mesh.faces = [ get_outer_facets(solids); [ elem for elem in mesh.elems if elem.shape.ndim==2 ] ]
            mesh.edges = getedges(mesh.faces)
        end
    end

    # Quality
    Q = Float64[]
    for c in mesh.elems
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end

    # tag
    tags = sort(unique([elem.tag for elem in mesh.elems]))
    tag_dict = Dict( tag=>i-1 for (i,tag) in enumerate(tags) )
    T = [ tag_dict[elem.tag] for elem in mesh.elems ] 

    # Ordering
    reorder && reorder!(mesh)

    # Reset data
    empty!(mesh.node_data)
    empty!(mesh.elem_data)
    mesh.node_data["node-id"]   = collect(1:length(mesh.nodes))
    mesh.elem_data["quality"]   = Q
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))
    mesh.elem_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.elems ]
    mesh.elem_data["tag"] = T


    return nothing
end

function compute_facets!(mesh::Mesh)
    if mesh.env.ndim==2
        mesh.edges = get_outer_facets(mesh.elems)
        mesh.faces = mesh.edges
    elseif mesh.env.ndim==3
        solids = [ elem for elem in mesh.elems if elem.shape.ndim==3 ]
        mesh.faces = [ get_outer_facets(solids); [ elem for elem in mesh.elems if elem.shape.ndim==2 ] ]
        mesh.edges = getedges(mesh.faces)
    end
end

function syncronize!(mesh::Mesh; reorder=false, cleandata=false)
    @assert mesh.env.ndim!=0

    # Numberig nodes
    for (i,p) in enumerate(mesh.nodes) 
        p.id = i 
    end

    # Numberig cells and setting env
    for (i,elem) in enumerate(mesh.elems)
        elem.id = i
        elem.env = mesh.env
    end

    # Faces and edges
    compute_facets!(mesh)

    # Quality
    Q = Float64[]
    for c in mesh.elems
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end

    # Ordering
    reorder && reorder!(mesh)

    # update data
    if cleandata
        empty!(mesh.node_data)
        empty!(mesh.elem_data)
    end

    # tag
    tags = sort(unique([elem.tag for elem in mesh.elems]))
    tag_dict = Dict( tag=>i-1 for (i,tag) in enumerate(tags) )
    
    T = [ tag_dict[elem.tag] for elem in mesh.elems ] 
    mesh.node_data["node-id"]   = collect(1:length(mesh.nodes))
    mesh.elem_data["quality"]   = Q
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))
    mesh.elem_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.elems ]
    mesh.elem_data["tag"] = T
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

    mesh.env.ndim = max(mesh.env.ndim, m2.env.ndim)

    pointdict = Dict{UInt, Node}()
    for m in (mesh, m2)
        for node in m.nodes
            pointdict[hash(node)] = node
        end    
    end
    nodes = collect(values(pointdict))

    elems = Cell[]
    for m in (mesh, m2)
        for elem in m.elems
            newelemnodes = Node[pointdict[hash(node)] for node in elem.nodes ]
            newelem = Cell(elem.shape, newelemnodes, tag=elem.tag)
            push!(elems, newelem)
        end
    end

    mesh.nodes = nodes
    mesh.elems = elems
    mesh._pointdict = pointdict

    fixup!(mesh, reorder=false)

    return nothing
end


function join_meshes(m1::Mesh, m2::Mesh)
    mesh = copy(m1)
    join_mesh!(mesh, m2)
    return mesh
end


"""
    $(SIGNATURES)

Creates a `Mesh` from a nodal `coordinates` matrix,
an array of `connectivities` and a list of `cellshapes`.
If `cellshapes` are not provided, they are computed based on the geometry.
A `tag` string for all generated cells can be provided optionally.

# Examples
    
```jldoctest
julia> using Amaru;
julia> C = [ 0 0; 1 0; 1 1; 0 1 ];
julia> L = [ [ 1, 2, 3 ], [ 1, 3, 4 ] ];
julia> S = [ TRI3, TRI3 ];
julia> Mesh(C, L, S, tag="triangle")

Mesh
  ndim: 2
  nodes: 4-element Vector{Node}:
    1: Node  id=1
    2: Node  id=2
    3: Node  id=3
    4: Node  id=4
  elems: 2-element Vector{Cell}:
    1: Cell  id=1  tag="triangle"
    2: Cell  id=2  tag="triangle"
  faces: 4-element Vector{Cell}:
    1: Cell  id=-1  tag="triangle"
    2: Cell  id=-1  tag="triangle"
    3: Cell  id=-1  tag="triangle"
    4: Cell  id=-1  tag="triangle"
  edges: 4-element Vector{Cell}:
    1: Cell  id=-1  tag="triangle"
    2: Cell  id=-1  tag="triangle"
    3: Cell  id=-1  tag="triangle"
    4: Cell  id=-1  tag="triangle"
  node_data: OrderedDict{String, Array} with 1 entry
    "node-id" => [1, 2, 3, 4]
  elem_data: OrderedDict{String, Array} with 3 entries
    "quality" => [0.8915188114208271, 0.8915188114208271]
    "elem-id" => [1, 2]
    "cell-type" => [5, 5]
```
"""
function Mesh(
              coordinates::Array{<:Real},
              conns      ::Array{Array{Int64,1},1},
              cellshapes ::Array{CellShape,1}=CellShape[];
              tag        ::String="",
              quiet      ::Bool=true,
             )
    

    n = size(coordinates, 1) # number of nodes
    m = size(conns , 1) # number of cells

    nodes = Node[]
    for i in 1:n
        C = coordinates[i,:]
        push!(nodes, Node(C))
    end

    # Get ndim
    ndim = size(coordinates,2)
    env  = MeshEnv(ndim)

    cells = Cell[]
    for i in 1:m
        pts = nodes[conns[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        cell = Cell(shape, pts, tag=tag, env=env)
        push!(cells, cell)
    end

    mesh = Mesh(ndim)
    mesh.nodes = nodes
    mesh.elems = cells
    fixup!(mesh, reorder=false) # no node ordering

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
    $(SIGNATURES)

Generates a `Mesh` based from a set of `items` that can be 
`Block` or `Mesh` objects, or a combination of both.
If `reorder` is `true` (default), the nodes are ordered using the Cuthill–McKee algorithm.
The mesh generation steps are printed in `stdout`.
The printed output can be set to `verbose` or `silent`.

# Examples
    
```@example
julia> using Amaru;
julia> B1 = Block([0 0; 1 1], nx=2, ny=2);
julia> B2 = Block([1 0; 2 1], nx=3, ny=2);
julia> Mesh(B1, B2)
Mesh
  ndim: 2
  nodes: 18-element Vector{Node}:
    1: Node  id=1
    2: Node  id=2
    3: Node  id=3
    4: Node  id=4
    ⋮
    15: Node  id=15
    16: Node  id=16
    17: Node  id=17
    18: Node  id=18
  elems: 10-element Vector{Cell}:
    1: Cell  id=1
    2: Cell  id=2
    3: Cell  id=3
    4: Cell  id=4
    ⋮
    7: Cell  id=7
    8: Cell  id=8
    9: Cell  id=9
    10: Cell  id=10
  faces: 14-element Vector{Cell}:
    1: Cell  id=-1
    2: Cell  id=-1
    3: Cell  id=-1
    4: Cell  id=-1
    ⋮
    11: Cell  id=-1
    12: Cell  id=-1
    13: Cell  id=-1
    14: Cell  id=-1
  edges: 14-element Vector{Cell}:
    1: Cell  id=-1
    2: Cell  id=-1
    3: Cell  id=-1
    4: Cell  id=-1
    ⋮
    11: Cell  id=-1
    12: Cell  id=-1
    13: Cell  id=-1
    14: Cell  id=-1
  node_data: OrderedDict{String, Array} with 1 entry
    "node-id" => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
  elem_data: OrderedDict{String, Array} with 3 entries
    "quality" => 10-element Vector{Float64}
    "elem-id" => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    "cell-type" => [9, 9, 9, 9, 9, 9, 9, 9, 9, 9]
```
"""
function Mesh(
    items     ::Union{Mesh, AbstractBlock, Array{<:Union{AbstractBlock, Array},1}}...;
    ndim      ::Int = 0,
    genfacets ::Bool = true,
    reorder   ::Bool = true,
    quiet     ::Bool = false,
)

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

    # Find ndim
    if length(blocks)>0
        blndim = maximum(bl.ndim for bl in blocks)
        ndim = max(ndim, blndim)
    end
    if length(meshes)>0
        msndim = maximum(mesh.env.ndim for mesh in meshes)
        ndim = max(ndim, msndim)
    end

    nmeshes = length(meshes)
    nblocks = length(blocks)
    if !quiet
        printstyled("Mesh generation:\n", bold=true, color=:cyan)
        nmeshes>0 && @printf "  %5d meshes\n" nmeshes
        @printf "  %5d blocks\n" nblocks
    end

    # New mesh object
    mesh = Mesh(ndim)

    # Join meshes
    for m in meshes
        join_mesh!(mesh, m)
    end

    # Split blocks: generates nodes and cells
    for (i,b) in enumerate(blocks)
        split_block(b, mesh)
        quiet || print("  spliting block ", i, "...    \r")
    end

    # Updates numbering, quality, facets and edges
    fixup!(mesh, genfacets=genfacets, reorder=reorder)

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        @printf "  %4dd mesh                             \n" mesh.env.ndim
        @printf "  %5d nodes\n" npoints
        @printf "  %5d cells\n" ncells
    end
    if !quiet
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        if genfacets
            @printf "  %5d faces\n" nfaces
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

export stats
function stats(mesh::Mesh)
    printstyled("Mesh stats:\n", bold=true, color=:cyan)

    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    @printf "  %3dd mesh                             \n" mesh.env.ndim
    @printf "  %4d nodes\n" npoints
    @printf "  %4d cells\n" ncells

    L = Float64[]
    for elem in mesh.elems.solids
        l = cell_extent(elem)^(1/mesh.env.ndim)
        push!(L, l)
    end
    lavg = mean(L)
    lmdn = quantile(L, 0.5)
    minl = minimum(L)
    maxl = maximum(L)

    bin = (maxl-minl)/10
    hist  = fit(Histogram, L, minl:bin:maxl, closed=:right).weights
    # @show hist
    lmod = (findmax(hist)[2]-1)*bin + bin/2

    @printf "  lmin = %7.5f\n" minl
    @printf "  lmax = %7.5f\n" maxl
    @printf "  lavg = %7.5f\n" lavg
    @printf "  lmdn = %7.5f\n" lmdn
    @printf "  lmod = %7.5f\n" lmod
end


function Mesh(elems::Array{Cell,1})
    length(elems)==0 && return Mesh()

    nodes = copy.(getnodes(elems))
    # digs = 8
    # for node in nodes
    #     node.coord = round.(node.coord, digits=digs)
    #     # round!(node.coord, digits=digs)
    # end

    newmesh = Mesh(getndim(nodes))
    newmesh.nodes = nodes

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
    length(elems)==0 && return Mesh()

    elemids = [ elem.id for elem in elems ]

    nodedict = OrderedDict{Int,Node}()
    for elem in elems
        for node in elem.nodes
            nodedict[node.id] = copy(node)
        end
    end
    nodeids  = collect(keys(nodedict))
    nodes = collect(values(nodedict))
    
    newmesh = Mesh(getndim(nodes))
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


# function slice(mesh::Mesh, normal::Array{Float64,1})
    # mesh = Mesh()
    # return mesh
# end


# function clip(mesh::Mesh, normal::Array{Float64,1})
    # mesh = Mesh()
    # return mesh
# end


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
        vals = [ mean(data[i]) for i in 1:ncells ]
    end

    # filter cells
    cells = Cell[]
    for (cell, val) in zip(mesh.elems, vals)
        if minval <= val <= maxval
            push!(cells, cell)
        end
    end

    # get nodes
    nodes = getnodes(cells)

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
function get_segment_data(msh::AbstractDomain, X1::Array{<:Real,1}, X2::Array{<:Real,1}, filename::String=""; n=200)
    data = msh.node_data
    table = DataTable(["s"; collect(keys(data))])
    X1 = [X1; 0.0][1:3]
    X2 = [X2; 0.0][1:3]
    Δ = (X2-X1)/(n-1)
    Δs = norm(Δ)
    s1 = 0.0

    for i in 1:n
        X = X1 + Δ*(i-1)
        s = s1 + Δs*(i-1)
        cell = find_elem(X, msh.elems, msh._elempartition, 1e-7, Cell[])
        coords =getcoords(cell)
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
function randmesh(n::Int...; cellshape=nothing)
    ndim = length(n)
    if ndim==2
        lx, ly = (1.0, 1.0)
        nx, ny = n
        isnothing(cellshape) && (cellshape=rand((TRI3, TRI6, QUAD4, QUAD8)))
        m = Mesh(Block([0.0 0.0; lx ly], nx=nx, ny=ny, cellshape=cellshape), quiet=true)
    else
        lx, ly, lz = (1.0, 1.0, 1.0)
        nx, ny, nz = n
        isnothing(cellshape) && (cellshape=rand((TET4, TET10, HEX8, HEX20)))
        m = Mesh(Block([0.0 0.0 0.0; lx ly lz], nx=nx, ny=ny, nz=nz, cellshape=cellshape), quiet=true)
    end
end


function Mesh(geo::GeoModel; recombine=false, size=0.1, quadratic=false, quiet=false, algorithm::Symbol=:delaunay)

    if !quiet
        printstyled("Unstructured mesh generation:\n", bold=true, color=:cyan)
        nsurfs = length(geo.surfaces)
        nvols  = length(geo.volumes)
        @printf "  %5d surfaces\n" nsurfs
        nvols>0 && @printf "  %5d volumes\n" nvols
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t1")

    isvolumemesh = length(geo.volumes)>0

    # add points
    # ndim = 2
    for p in geo.points
        sz = p.size==0 ? size : p.size
        gmsh.model.geo.addPoint(p.coord.x, p.coord.y, p.coord.z, sz, p.id)
        # p.coord.z != 0.0 && (ndim=3)
    end

    # add lines
    for l in geo.lines
        if l isa Line
            p1 = l.points[1].id
            p2 = l.points[2].id
            gmsh.model.geo.addLine(p1, p2, l.id)
        else #Arc
            p1 = l.points[1].id
            pc = l.points[2].id
            p2 = l.points[3].id
            gmsh.model.geo.addCircleArc(p1, pc, p2, l.id)
        end
    end

    function get_loop_idxs(loop::AbstractLoop)
        l_idxs = [ l.id for l in loop.lines ]
        l_conn = [ [ p.id for p in l.points] for l in loop.lines ]

        if !(l_conn[1][end] in l_conn[2])
            l_idxs[1] *= -1
        end


        for i in 2:length(l_conn)
            last_pidx = l_idxs[i-1]>0 ? l_conn[i-1][end] : l_conn[i-1][1]
            if l_conn[i][1]!=last_pidx
                l_idxs[i] *= -1
            end
        end

        return l_idxs
    end

    # add loops
    for loop in geo.loops
        line_idxs = get_loop_idxs(loop)
        gmsh.model.geo.addCurveLoop(line_idxs, loop.id)
    end

    # add surfaces
    for surf in geo.surfaces
        lo_idxs = [ lo.id for lo in surf.loops ]
        if surf isa PlaneSurface
            gmsh.model.geo.addPlaneSurface(lo_idxs, surf.id) # plane surface
        else
            gmsh.model.geo.addSurfaceFilling(lo_idxs, surf.id) # filling surf
        end
    end
    surf_idxs = [ surf.id for surf in geo.surfaces ]

    # add volumes
    for vol in geo.volumes
        surf_idxs = [ surf.id for surf in vol.surfaces ]

        gmsh.model.geo.addSurfaceLoop(surf_idxs, vol.id)
        gmsh.model.geo.addVolume([vol.id], vol.id) # not considering volume holes
    end
    vol_idxs = [ vol.id for vol in geo.volumes ]

    gmsh.model.geo.synchronize() # only after geometry entities are defined
    
    for l in geo.lines
        # transfinite
        if l.n>0
            gmsh.model.mesh.set_transfinite_curve(l.id, l.n+1)
        end
    end
    
    # generate mesh
    if algorithm==:delaunay
        gmsh.option.setNumber("Mesh.Algorithm", 5)
    elseif algorithm==:frontal
        gmsh.option.setNumber("Mesh.Algorithm", 6)
    else
        error("Mesh: Wrong algorithm")
    end

    for s in geo.surfaces
        s.recombine && gmsh.model.mesh.set_recombine(2, s.id)
        s.transfinite && gmsh.model.mesh.set_transfinite_surface(s.id)
    end

    if !isvolumemesh
        tagset = Set([ surf.tag for surf in geo.surfaces ])
        tagsdict = Dict( tag=>i for (i,tag) in enumerate(tagset) )
        for (tag, gidx) in tagsdict
            surf_idxs = [ surf.id for surf in geo.surfaces if surf.tag==tag ]
            gmsh.model.addPhysicalGroup(2, surf_idxs, gidx) # ndim, entities, group_id
        end

        # gmsh.model.addPhysicalGroup(2, surf_idxs) # ndim, entities
        # gmsh.model.mesh.generate(2)
    else
        tagset = Set([ vol.tag for vol in geo.volumes  ])
        tagsdict = Dict( tag=>i for (i,tag) in enumerate(tagset) )
        for (tag, gidx) in tagsdict
            vol_idxs = [ vol.id for vol in geo.volumes if vol.tag==tag]
            gmsh.model.addPhysicalGroup(3, vol_idxs, gidx) # ndim, entities, group_id
        end

        # tagsdict = Dict( tag=>i for (i,tag) in enumerate(tagset) )
        # for v in geo.entities
        #     v isa Volume || continue
        #     gmsh.model.addPhysicalGroup(3, [v.id], tagsdict[v.tag]) # ndim, entities
        # end

        # xx1 = gmsh.model.addPhysicalGroup(3, vol_idxs, 100) # ndim, entities
        # @show xx1
        # gmsh.model.mesh.generate(3)
    end

    # embed points
    for p in geo.points
        length(p.lines)==0 || continue
    
        # search surfaces
        for s in geo.surfaces
            if inside(p, s.loops[1])
                gmsh.model.mesh.embed(0,[p.id],2,s.id)
            end
        end

    end

    tempfile = "_temp.vtk"
    logfile = "_gmsh.log"
    try
        open(logfile, "w") do out
            redirect_stdout(out) do
                gmsh.model.mesh.generate(isvolumemesh ? 3 : 2)
                quadratic && gmsh.model.mesh.setOrder(2) # quadratic elements
                recombine && gmsh.model.mesh.recombine()
                gmsh.write(tempfile)
                gmsh.write("file.geo_unrolled")
            end
        end
    catch err
        error("Error generating unstructured mesh.")
    end
    
    gmsh.finalize()
    mesh = Mesh(tempfile)
    rm(tempfile, force=true)
    rm(logfile, force=true)

    # flip elements
    for elem in mesh.elems
        isinverted(elem) && flip!(elem)
    end

    # set tags for elements
    if !isvolumemesh
        invtagsdict = Dict( i=>tag for (tag,i) in tagsdict )
        for elem in mesh.elems
            elem.tag = invtagsdict[ mesh.elem_data["CellEntityIds"][elem.id] ]
        end
    else
        invtagsdict = Dict( i=>tag for (tag,i) in tagsdict )
        for elem in mesh.elems
            elem.tag = invtagsdict[ mesh.elem_data["CellEntityIds"][elem.id] ]
        end
    end

    # set tag for nodes
    ptagdict = Dict()
    for p in geo.points
        p.tag!="" || continue
        ptagdict[p.coord] = p.tag
    end
    
    for node in mesh.nodes
        tag = get(ptagdict, node.coord, "")
        tag != "" || continue
        node.tag = tag
    end

    fixup!(mesh, reorder=true)

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        @printf "  %4dd mesh                             \n" mesh.env.ndim
        @printf "  %5d nodes\n" npoints
        @printf "  %5d cells\n" ncells

        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        @printf "  %5d faces\n" nfaces
        @printf "  %5d surface edges\n" nedges
    end

    return mesh

end
