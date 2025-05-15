# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


mutable struct Mesh<:AbstractDomain
    nodes::Array{Node,1} 
    elems::Array{Cell,1}
    faces::Array{Cell,1}
    edges::Array{Cell,1}
    node_data::OrderedDict{String,Array}
    elem_data ::OrderedDict{String,Array}
    ctx::MeshEnv

    _pointdict::Dict{UInt64,Node}
    _elempartition::ElemPartition

    function Mesh(ndim=0)
        this = new()
        this.nodes = []
        this.elems  = []
        this.faces  = []
        this.edges  = []
        # this.ndim   = 0
        this.node_data = OrderedDict()
        this.elem_data = OrderedDict()
        this.ctx = MeshEnv(ndim)
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
    ndim = mesh.ctx.ndim
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
    face_d = OrderedDict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in getfacets(cell)
            hs = hash(face)
            if haskey(face_d, hs)
                delete!(face_d, hs)
            else
                face_d[hs] = face
            end
        end
    end

    return CellFace[ face for face in values(face_d) ]
end


function get_outer_facets_by_id(cells::Array{<:AbstractCell,1})
    # face_d = OrderedDict{UInt64, Cell}()
    face_d = Dict()
    hash1(edge) = sort([ n.id for n in edge.nodes ])
    @show "hiiiii"

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in getfacets(cell)
            hs = hash1(face)
            if haskey(face_d, hs)
                delete!(face_d, hs)
            else
                face_d[hs] = face
            end
        end
    end

    return collect(values(face_d))
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
function sortnodes!(mesh::Mesh; sort_degrees=true, reversed=true)

    # Get all mesh edges
    all_edges = Cell[]  # using array ( the use of dictionary got hash collisions for some embedded cells)
    for cell in mesh.elems

        # adding cell edges
        if cell.shape.family == BULKCELL #is_solid(cell.shape)
            for edge in getedges(cell)
                push!(all_edges, edge)
            end

            #check for lagrangian elements
            if cell.shape==QUAD9
                edge = Cell(POLYVERTEX, [ cell.nodes[1], cell.nodes[end] ] )
                push!(all_edges, edge)
            end
            if cell.shape==HEX27
                edge = Cell(POLYVERTEX, [ cell.nodes[1], cell.nodes[end] ] )
                push!(all_edges, edge)
                for i in (21,22,23,24,25,26)
                    edge = Cell(POLYVERTEX, [ cell.nodes[i], cell.nodes[end] ] )
                    push!(all_edges, edge)
                end
            end
            continue
        end

        # joint1D cells (semi-embedded approach)
        if cell.shape in (JLINK2, JLINK3)
            npts = cell.shape.npoints
            edge = Cell(POLYVERTEX, [ cell.nodes[1]; cell.nodes[end-npts+1:end] ])
            push!(all_edges, edge)
            continue
        end

        # embedded line cells
        if cell.shape.family == LINECELL && length(cell.linked_elems)>0
            edge1 = Cell(cell.shape, cell.nodes)
            edge2 = Cell(LIN2, [ cell.nodes[1], cell.linked_elems[1].nodes[1] ])
            push!(all_edges, edge1)
            push!(all_edges, edge2)
            continue
        end

        # all other cells
        edge = Cell(cell.shape, cell.nodes)
        push!(all_edges, edge)

    end

    # Get neighbors ids
    nnodes = length(mesh.nodes)
    neighs_ids = Array{Int64}[ [] for i in 1:nnodes ]

    # for edge in values(all_edges)
    for edge in all_edges
        nodes = edge.nodes
        np = length(nodes)

        for i in 1:np-1
            for j in i+1:np
                push!(neighs_ids[nodes[i].id], nodes[j].id)
                push!(neighs_ids[nodes[j].id], nodes[i].id)
            end
        end
    end

    # remove duplicates
    neighs_ids = [ unique(list) for list in neighs_ids ]

    # get neighbors
    neighs = Array{Node}[ mesh.nodes[list] for list in neighs_ids ]

    # list of degrees per point
    degrees = Int64[ length(list) for list in neighs]
    mindeg, idx  = findmin(degrees)

    if mindeg == 0
        # Case of overlapping elements where edges have at least one point with the same coordinates
        notify("sortnodes!: Reordering nodes failed. Possible causes: disconnected domain or non used nodes.")
        return
    end

    N = [ mesh.nodes[idx] ] # new list of nodes
    L = Dict{Int64,Node}()  # last levelset. Use ids as keys instead of hash to avoid collisions of nodes with same coordinates
    L[idx] = mesh.nodes[idx]
    LL = Dict{Int64,Node}()  # levelset before the last one

    while length(N) < nnodes
        # Generating current levelset A
        A = Dict{Int64,Node}()

        for p in values(L)
            for q in neighs[p.id]
                (haskey(L, q.id) || haskey(LL, q.id)) && continue
                A[q.id] = q
            end
        end
        if length(A)==0
            #@error "sortnodes!: Reordering nodes failed! Possible error with cell connectivities."
            notify("sortnodes!: Reordering nodes failed. Possible causes: disconnected domain, non used nodes and overlapping cells.")
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

    ids = [ node.id for node in N ]

    # update nodal fields
    for (key, data) in mesh.node_data
        key == "node-id" && continue
        sz = size(data)
        if length(sz)==1
            mesh.node_data[key] = data[ids]
        else
            mesh.node_data[key] = data[ids,:]
        end
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


function compute_facets!(mesh::Mesh)
    if mesh.ctx.ndim==2
        mesh.edges = get_outer_facets(mesh.elems)
        mesh.faces = mesh.edges
    elseif mesh.ctx.ndim==3
        solids = filter( elem->elem.shape.ndim==3, mesh.elems )
        planars = filter( elem->elem.shape.ndim==2 && elem.shape.family==BULKCELL, mesh.elems )
        for cell in planars # TODO: fix for BC in shells
            cell.owner = cell # set itself as owner
        end
        mesh.faces = [ get_outer_facets(solids); planars ]
        mesh.edges = getedges(mesh.faces)
    end
end


# Syncs the mesh data
function synchronize!(mesh::Mesh; sortnodes=false, cleandata=false)

    ndim = mesh.ctx.ndim
    if ndim!=3
        ndim = any( node.coord[3] != 0.0 for node in mesh.nodes ) ? 3 : 2
        if ndim == 2
            ndim = any( node.coord[2] != 0.0 for node in mesh.nodes ) ? 2 : 1
        end
    end

    mesh.ctx.ndim = max(ndim, mesh.ctx.ndim)

    # Numberig nodes
    for (i,p) in enumerate(mesh.nodes) 
        p.id = i 
    end

    # Numberig cells and setting ctx
    for (i,elem) in enumerate(mesh.elems)
        elem.id = i
        elem.ctx = mesh.ctx
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
    sortnodes && sortnodes!(mesh)

    # update data
    if cleandata
        empty!(mesh.node_data)
        empty!(mesh.elem_data)
    end

    # tag
    # tags = sort(unique([elem.tag for elem in mesh.elems]))
    # tag_dict = Dict( tag=>i-1 for (i,tag) in enumerate(tags) )
    
    mesh.node_data["node-id"]   = collect(1:length(mesh.nodes))
    mesh.elem_data["quality"]   = Q
    mesh.elem_data["elem-id"]   = collect(1:length(mesh.elems))
    mesh.elem_data["cell-type"] = [ Int(cell.shape.vtk_type) for cell in mesh.elems ]
    
    # tag
    # tags  = collect(Set(elem.tag for elem in mesh.elems))
    # tag_d = Dict( tag=>i for (i,tag) in enumerate(tags) )
    # T     = Int[ tag_d[elem.tag]-1 for elem in mesh.elems ]
    # mesh.elem_data["tag"] = T
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

    # mesh.ctx.ndim = max(mesh.ctx.ndim, m2.ctx.ndim)

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

    synchronize!(mesh, sortnodes=false)

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
If `cellshapes` are not provided, they are guessed based on the geometry.
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
    ctx  = MeshEnv(ndim)

    cells = Cell[]
    for i in 1:m
        pts = nodes[conns[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        cell = Cell(shape, pts, tag=tag, ctx=ctx)
        push!(cells, cell)
    end

    mesh = Mesh(ndim)
    mesh.nodes = nodes
    mesh.elems = cells
    synchronize!(mesh, sortnodes=false) # no node ordering

    return mesh
end



export stats
function stats(mesh::Mesh)
    printstyled("Mesh stats:\n", bold=true, color=:cyan)

    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    @printf "  %3dd mesh                             \n" mesh.ctx.ndim
    @printf "  %4d nodes\n" npoints
    @printf "  %4d cells\n" ncells

    L = Float64[]
    for elem in mesh.elems.solids
        l = cell_extent(elem)^(1/mesh.ctx.ndim)
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


# Construct a mesh given a set of elements
function Mesh(elems::Vector{Cell})
    length(elems)==0 && return Mesh()

    nodes = [ node for elem in elems for node in elem.nodes ]

    # check if all nodes have id
    all( node.id > 0 for node in nodes ) || error("Mesh: all input nodes must have a valid id")

    # newnodes
    nodes_d  = Dict{Int,Node}( node.id => copy(node) for node in nodes )
    newnodes = collect(values(nodes_d))
    node_ids = collect(keys(nodes_d))

    # copy elements
    newelems = Cell[]
    for elem in elems
        elem_nodes = Node[ nodes_d[node.id] for node in elem.nodes ]
        newelem    = Cell(elem.shape, elem_nodes, tag=elem.tag)
        push!(newelems, newelem)
    end

    ndim = any( node.coord[3] != 0.0 for node in newnodes ) ? 3 : 2

    newmesh       = Mesh(ndim)
    newmesh.nodes = newnodes
    newmesh.elems = newelems

    synchronize!(newmesh, sortnodes=true)
    return newmesh
end


function Base.getindex(mesh::Mesh, filter::Union{String,Expr,Symbol,Symbolic,CellFamily,CellShape})
    elems = mesh.elems[filter]
    length(elems)==0 && return Mesh()

    # newnodes
    nodes    = [ node for elem in elems for node in elem.nodes ]
    nodes_d  = Dict{Int,Node}( node.id => copy(node) for node in nodes )
    newnodes = collect(values(nodes_d))
    node_ids = collect(keys(nodes_d))
    ndim     = any( node.coord[3] != 0.0 for node in newnodes ) ? 3 : 2

    # copy elements
    newelems = Cell[]
    for elem in elems
        elem_nodes = Node[ nodes_d[node.id] for node in elem.nodes ]
        newelem    = Cell(elem.shape, elem_nodes, tag=elem.tag)
        push!(newelems, newelem)
    end
    
    # create new mesh
    newmesh       = Mesh(ndim)
    newmesh.nodes = newnodes
    newmesh.elems = newelems

    # update node fields
    for (key, data) in mesh.node_data
        sz = size(data)
        if length(sz)==1
            newmesh.node_data[key] = data[node_ids]
        else
            newmesh.node_data[key] = data[node_ids,:]
        end
    end
    
    # update elem fields
    elem_ids = [ elem.id for elem in elems ]
    for (key, data) in mesh.elem_data
        newmesh.elem_data[key] = data[elem_ids]
    end

    synchronize!(newmesh, sortnodes=true)

    return newmesh
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
    synchronize!(new_mesh, sortnodes=false)

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





