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


function compute_facets!(mesh::Mesh)
    if mesh.ctx.ndim==2
        mesh.edges = get_outer_facets(mesh.elems)
        mesh.faces = mesh.edges
    elseif mesh.ctx.ndim==3
        solids = filter( elem->elem.shape.ndim==3, mesh.elems )
        planars = filter( elem->elem.shape.ndim==2, mesh.elems )
        for cell in planars
            cell.owner = cell # set itself as owner
        end
        mesh.faces = [ get_outer_facets(solids); planars ]
        mesh.edges = getedges(mesh.faces)
    end
end


# Syncs the mesh data
function syncronize!(mesh::Mesh; reorder=false, cleandata=false)
    sumz = map( node->abs(node.coord.z), mesh.nodes ) |> sum
    if sumz==0
        sumy = map( node->abs(node.coord.y), mesh.nodes ) |> sum
        ndim = sumy==0 ? 1 : 2
    else
        ndim = 3
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

    syncronize!(mesh, reorder=false)

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
    syncronize!(mesh, reorder=false) # no node ordering

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

    syncronize!(newmesh, reorder=false)

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

    syncronize!(newmesh, reorder=false)

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
    syncronize!(new_mesh, reorder=false)

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





