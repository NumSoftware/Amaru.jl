# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Cell


"""
    Cell

A type that represents a mesh cell. It can be a 1D, 2D or 3D.

# Fields
$(FIELDS)
"""
mutable struct Cell<:AbstractCell
    id     ::Integer
    shape  ::CellShape
    nodes  ::Array{Node,1}
    tag    ::String
    active ::Bool
    quality::Float64              # quality index: surf/(reg_surf)
    embedded::Bool                # flag for embedded cells
    crossed::Bool                 # flag if cell crossed by linear inclusion
    owner  ::Union{AbstractCell,Nothing}  # owner cell if this cell is a face/edge
    linked_elems::Array{AbstractCell,1}   # neighbor cells in case of joint cell
    ctx::MeshEnv                 # mesh environment variables

    @doc """
        $(SIGNATURES)

    Constructs a `Cell` given a `shape` and an array of nodes.
    A `tag` string can be provided optionally.
    If the cell is a face or an edge, the `owner` element cell can be provided.

    # Examples
    
    ```jldoctest
    julia> using Amaru;
    julia> nodes = [ Node(0,0), Node(1,0), Node(1,1) ];
    julia> Cell(TRI3, nodes, tag="triangle");
    Cell
      id: -1
      shape: CellShape  name="TRI3"
      nodes: 3-element Vector{Node}:
          1: Node  id=-1
          2: Node  id=-1
          3: Node  id=-1
      tag: "triangle"
      quality: 0.0
      embedded: false
      crossed: false
      owner: nothing
      linked_elems: 0-element Vector{Amaru.AbstractCell}
    ```
    """
    function Cell(shape::CellShape, nodes::Array{Node,1}; ctx::MeshEnv=MeshEnv(0), tag::String="", owner=nothing, id::Int=-1, active=true)
        this = new()
        this.id = id
        this.shape = shape
        this.nodes = copy(nodes)
        this.tag = tag
        this.active  = active
        this.quality = 0.0
        this.embedded= false
        this.crossed = false
        this.owner   = owner
        this.linked_elems = []
        this.ctx     = ctx
        return this
    end
end

const CellFace=Cell
const CellEdge=Cell
const Facet=Cell


### Cell methods

Base.hash(cell::Cell) = sum(hash(node) for node in cell.nodes)
Base.isequal(c1::Cell, c2::Cell) = hash(c1)==hash(c2)

"""
    $(TYPEDSIGNATURES)

Creates a copy of `cell`.
"""
Base.copy(cell::Cell)  = Cell(cell.shape, cell.nodes, tag=cell.tag, owner=cell.owner)

"""
    $(SIGNATURES)

Returns a matrix with the nodal coordinates of `c`.
`ndim` (default 3) can be used to define the number axes.
"""
function getcoords(c::AbstractCell, ndim=3)
    n = length(c.nodes)
    C = Array{Float64}(undef, n, ndim)
    for (i,p) in enumerate(c.nodes)
        C[i,1] = p.coord.x
        if ndim>1 C[i,2] = p.coord.y end
        if ndim>2 C[i,3] = p.coord.z end
    end
    return C
end


# Return all nodes in cells
# function getnodes(cells::Array{<:AbstractCell,1})
#     return collect(Set(node for cell in cells for node in cell.nodes)) 
# end



function getnodes(elems::Vector{<:AbstractCell})
    # get all nodes using object id
    points_d = Dict{UInt64, Node}()
    for elem in elems
        for node in elem.nodes
            points_d[objectid(node)] = node
        end
    end

    return collect(values(points_d))
end


# function (cells::Array{<:AbstractCell,1})(filter=nothing; tag=nothing, shape=nothing, family=nothing, plane=nothing, invert=false)

#     filtered = collect(1:length(cells))
    
#     if isa(filter, Expr) && filter!=:()
#         indexes = copy(filtered)
#         filtered = Int[]
#         # cells id must be set
#         nodes = getnodes(cells)
#         pointmap = zeros(Int, maximum(node.id for node in nodes) ) # points and pointmap may have different sizes

#         T = Bool[]
#         for (i,node) in enumerate(nodes)
#             pointmap[node.id] = i
#             x, y, z = node.coord.x, node.coord.y, node.coord.z
#             push!(T, evaluate(filter, x=x, y=y, z=z))
#         end
 
#         filtered = Int[ i for i in indexes if all( T[pointmap[node.id]] for node in cells[i].nodes ) ]
#     end

#     if isa(filter, Symbol)
#         families = Dict(:solids=>BULKCELL, :lines=>LINECELL, :joints=>JOINTCELL, :linejoints=>LINEJOINTCELL, :tipjoints=>TIPJOINT)
#         family = families[filter]
#     end

#     if tag!==nothing
#         indexes = copy(filtered)
#         filtered = Int[ i for i in indexes if cells[i].tag==tag ]
#     end
    
#     if shape!==nothing
#         indexes = copy(filtered)
#         filtered = Int[ i for i in indexes if cells[i].shape==shape ]
#     end

#     if family!==nothing
#         indexes = copy(filtered)
#         filtered = Int[ i for i in indexes if cells[i].shape.family==family ]
#     end

#     if plane!==nothing
#         indexes = copy(filtered)
#         filtered = Int[]

#         # plane normal
#         tol = 1e-5
#         A = plane # coplanar points
#         I = ones(size(A,1))
#         n = normalize(pinv(A.+tol/100)*I) # best fit normal        

#         for i in indexes
#             cell = cells[i]
#             cell.shape.family == JOINTCELL || continue
#             ndim = cell.shape.ndim+1
#             @assert ndim==size(A,2)
#             C = getcoords(cells[i], ndim)
#             I = ones(size(C,1))
#             N = pinv(C.+tol/100)*I # best fit normal 
#             norm(C*N-I)<tol || continue # check if cell is coplanar
#             normalize!(N)
#             norm(N-n)<tol || continue # check if planes are parallel
#             dot(A[1,:]-C[1,:], n)<tol && push!(filtered, i)
#         end
#     end

#     if invert
#         filtered = setdiff(1:length(cells), filtered)
#     end

#     return cells[filtered]
# end


function Base.getproperty(c::AbstractCell, s::Symbol)
    s == :coords && return getcoords(c)
    s == :faces  && return getfaces(c)
    s == :edges  && return getedges(c)
    s == :extent && return cell_extent(c)
    return getfield(c, s)
end


function Base.getproperty(cells::Array{<:AbstractCell,1}, s::Symbol)
    s in (:all, :solids, :bulks, :lines, :joints, :joints1D, :linejoints, :tipjoints, :embeddeds) && return cells[s]
    s == :nodes  && return getnodes(cells)
    s == :filter && return cells
    s == :active && return filter(cell -> cell.active, cells)
    return getfield(cells, s)
end


function Base.getindex(cells::Array{<:AbstractCell,1}, filters::NTuple; kwargs...)
    return getindex(cells, filters...; kwargs...)
end


function Base.getindex(
    cells::Array{<:AbstractCell,1}, 
    filters::Union{Symbol,Expr,Symbolic,String,CellFamily,CellShape}...;
    normal = nothing,
    plane  = nothing,
    invert = false
    )
    
    filtered = collect(1:length(cells))

    for filter in filters

        if isa(filter, Symbol)
            if filter in (:solids, :bulks) 
                filter=BULKCELL
            elseif filter == :lines
                filter=LINECELL
            elseif filter == :joints
                filter=JOINTCELL
            elseif filter in (:linejoints, :joints1D, :joints1d) 
                filter=LINEJOINTCELL
            elseif filter == :tipjoints
                filter=TIPJOINT
            end
        end

        if isa(filter, Expr) || isa(filter, Symbolic) 
            # nodes = getnodes(cells) # it merge nodes
            nodes = [ node for cell in cells for node in cell.nodes ]
            max_id = maximum( n->n.id, nodes )
            pointmap = zeros(Int, max_id) # points and pointmap may have different sizes

            T = Bool[]
            for (i,node) in enumerate(nodes)
                pointmap[node.id] = i
                x, y, z = node.coord.x, node.coord.y, node.coord.z
                push!(T, evaluate(filter, x=x, y=y, z=z))
            end

            # @show filtered
            # @show pointmap
            # @show [ n.id for n in nodes ]
            # for c in cells
            #     @show [ n.id for n in c.nodes ]
            # end
    
            filtered = Int[ i for i in filtered if all( T[pointmap[node.id]] for node in cells[i].nodes ) ]
        elseif isa(filter, Symbol)
            if filter==:all
                # do nothing (don't filter)
            elseif filter==:active
                filtered = Int[ i for i in filtered if cells[i].active ]
            elseif filter in (:embedded, :embeddeds)
                filtered = Int[ i for i in filtered if cells[i].shape.family==LINECELL && length(cells[i].linked_elems)>0 ]
            elseif filter == :shells
                filtered = Int[ i for i in filtered if cells[i].shape.family==BULKCELL && cells[i].shape.ndim==2 ]
            else
                error("getindex: cannot filter array of Cell with symbol $(repr(filter))")
            end
        elseif isa(filter, String)
            filtered = Int[ i for i in filtered if cells[i].tag==filter ]
        elseif isa(filter, CellFamily)
            filtered = Int[ i for i in filtered if cells[i].shape.family==filter ]
        elseif isa(filter, CellShape)
            filtered = Int[ i for i in filtered if cells[i].shape==filter ]
        end

    end

    
    if !isnothing(normal) || !isnothing(plane)
        indexes = copy(filtered)
        filtered = Int[]
        tol = 1e-8

        if !isnothing(normal)
            n = normalize(normal) # best fit normal
            dim = length(n)
        else
            # filter is a plane: a matrix with at least three points
            A = plane # coplanar points
            I = ones(size(A,1))
            n = normalize(pinv(A.+tol/100)*I) # best fit normal   
            dim = size(A,2)
        end

        for i in indexes
            cell = cells[i]

            ndim = cell.shape.ndim+1
            @assert dim==ndim
            C = getcoords(cells[i], ndim)
            N = cellnormal(cell)
            isnothing(N) && continue
            isparallelto(N,n) || continue

            if !isnothing(plane)
                dot(A[1,:]-C[1,:], n)<tol || continue # check if faces are coplanar
            end
            push!(filtered, i)
        end
    end

    if invert
        filtered = setdiff(1:length(cells), filtered)
    end

    return cells[filtered]
end


function cellnormal(cell::AbstractCell)
    iscoplanar(cell) || return nothing
    
    tol = 1e-8
    ndim = cell.shape.ndim+1
    C = getcoords(cell, ndim)
    I = ones(size(C,1))
    N = pinv(C.+tol/100)*I # best fit normal 
    normalize!(N)
    return N
end


function isparallelto(A,B)
    tol = 1e-8

    dotAB = dot(A,B)
    normAB = norm(A)*norm(B)
    abs(dotAB-normAB) < tol && return true
    abs(dotAB+normAB) < tol && return true
    return false
end


function iscoplanar(cell::AbstractCell)
    tol = 1e-8

    coords = getcoords(cell, 3)
    
    # find a plane
    X1 = coords[1,:]
    X2 = coords[2,:]
    X1X2 = X2-X1

    # look for a non-collinear point
    local X, N
    for i in 3:length(cell.nodes)
        X = coords[i,:]
        X1X = X-X1
        N = cross(X1X2, X1X)
        norm(N) > tol && break
    end
    
    # test the plane at each point
    for i in 3:length(cell.nodes)
        X = coords[i,:]

        if dot(X-X1, N) > tol
            return false
        end
    end

    return true
end


function nearest(cells::Array{Cell,1}, coord)
    n = length(cells)
    D = zeros(n)
    X = vec(coord)

    for (i,cell) in enumerate(cells)
        C = vec(mean(getcoords(cell), dims=1))
        D[i] = norm(X-C)
    end

    return cells[sortperm(D)[1]]
end



# Gets the coordinates of a bounding box for an array of nodes
function bounding_box(nodes::Array{<:AbstractPoint,1})
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for node in nodes
        node.coord.x < minx && (minx = node.coord.x)
        node.coord.y < miny && (miny = node.coord.y)
        node.coord.z < minz && (minz = node.coord.z)
        node.coord.x > maxx && (maxx = node.coord.x)
        node.coord.y > maxy && (maxy = node.coord.y)
        node.coord.z > maxz && (maxz = node.coord.z)
    end
    return [ minx miny minz; maxx maxy maxz ]
end


# Gets the coordinates of a bounding box for a cell
function bounding_box(cell::AbstractCell)
    return bounding_box(cell.nodes)
end


# Gets the coordinates of a bounding box for an array of cells
function bounding_box(cells::Array{<:AbstractCell,1})
    nodes = unique( Node[ p for c in cells for p in c.nodes ] )
    return bounding_box(nodes)
end


# gets all facets of a cell
function getfacets(cell::AbstractCell)
    faces  = Cell[]
    all_facets_idxs = cell.shape.facet_idxs
    facet_shape    = cell.shape.facet_shape

    facet_shape==() && return faces

    sameshape = typeof(facet_shape) == CellShape # check if all facets have the same shape

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_facets_idxs)
        nodes = cell.nodes[face_idxs]
        shape  = sameshape ? facet_shape : facet_shape[i]
        face   = Cell(shape, nodes, tag=cell.tag, owner=cell)
        face.nodes = nodes # update nodes since Cell creates a copy
        push!(faces, face)
    end

    return faces
end

# gets all faces of a cell
function getfaces(cell::AbstractCell)
    # Return a cell with the same shape in case of a 2D cell in 3D space
    # cell.ctx.ndim==3 && 
    cell.shape.ndim==2 && return [ Cell(cell.shape, cell.nodes, tag=cell.tag, owner=cell) ]
    return getfacets(cell)
end


# gets all edges of a cell
function getedges(cell::AbstractCell)
    edges  = Cell[]
    all_edge_idxs = cell.shape.edge_idxs

    for edge_idx in all_edge_idxs
        nodes = cell.nodes[edge_idx]
        shape  = (LIN2, LIN3, LIN4)[length(nodes)-1]
        edge   = Cell(shape, nodes, tag=cell.tag, owner=cell)
        push!(edges, edge)
    end

    return edges
end


# Returns the volume/area/length of a cell
function cell_extent(c::AbstractCell)
    IP = get_ip_coords(c.shape)
    nip = size(IP,1)
    nldim = c.shape.ndim # cell basic dimension

    # get coordinates matrix
    C = getcoords(c)
    J = Array{Float64}(undef, size(C,2), nldim)
    
    # calc metric
    vol = 0.0
    for i in 1:nip
        R    = IP[i].coord
        w    = IP[i].w
        dNdR = c.shape.deriv(R)
        @mul J = C'*dNdR
        normJ = norm2(J)
        #if normJ<0
            #error("cell_extent: Negative Jacobian while calculating cell volume/area/length id=$(c.id) shape=$(c.shape.name) ")
        #end
        vol += normJ*w
    end
    return vol
end


# Returns the surface/perimeter of a regular element given the volume/area of a cell
function regular_surface(metric::Float64, shape::CellShape)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ]
        A = metric
        a = 2.0*√(A/√3.0)
        return 3*a
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ]
        A = metric
        a = √A
        return 4*a
    end
    if shape in [ PYR5, PYR13 ]
        V = metric
        a = ( 3.0*√2.0*V )^(1.0/3.0)
        return (1 + √3.0)*a^2
    end
    if shape in [ TET4, TET10 ]
        V = metric
        a = ( 6.0*√2.0*V )^(1.0/3.0)
        return √3.0*a^2
    end
    if shape in [ HEX8, HEX20, HEX27 ]
        V = metric
        a = V^(1.0/3.0)
        return 6.0*a^2.0
    end
    if shape in [ WED6, WED15 ]
        V = metric
        a2 = (16.0/3.0*V^2)^(1.0/3.0)
        return (3.0 + √3.0/2.0)*a2
    end
    error("No regular surface/perimeter value for shape $(shape.name)")
end


# Returns the area/volume of a regular element given the perimeter/surface
function regular_volume(metric::Float64, shape::CellShape)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ]
        p = metric
        a = p/3
        return a^2/4*√3.0
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ]
        p = metric
        a = p/4
        return a^2
    end
    #if shape in [ PYR5 ]
        #V = metric
        #a = ( 3.0*√2.0*V )^(1.0/3.0)
        #return (1 + √3.0)*a^2
    #end
    if shape in [ TET4, TET10 ]
        s = metric
        A = s/4
        a = 2.0*√(A/√3.0)
        return a^3/(6*√2.0)
    end
    if shape in [ HEX8, HEX20 ]
        s = metric
        A = s/6
        a = √A
        return a^3
    end
    #if shape in [ WED6, WED15 ]
        #s = metric
        #return (3.0 + √3.0/2.0)*a2
    #end
    error("No regular area/volume value for shape $(get_name(shape))")
end



# Returns the cell quality ratio as vol/reg_vol
function cell_quality_2(c::AbstractCell)::Float64
    # get faces
    faces = getfacets(c)
    length(faces)==0 && return 1.0

    # cell surface
    surf = sum( cell_extent(f) for f in faces )

    # quality calculation
    vol = cell_extent(c) # volume or area
    rvol = regular_volume(surf, c.shape)
    return min(vol/rvol, 1.0)
end


# Returns the cell quality ratio as reg_surf/surf
function cell_quality(c::AbstractCell)::Float64
    # get faces
    c.shape.family==JOINTCELL && return 1.0

    faces = getfacets(c)
    length(faces)==0 && return 1.0

    # cell surface
    surf = 0.0
    for f in faces
        surf += cell_extent(f)
    end

    # quality calculation
    extent = cell_extent(c) # volume or area
    rsurf  = regular_surface(abs(extent), c.shape)

    k = 2
    q = sign(extent)*(rsurf/surf)^k
    if q>1
        q = 2 - q
    end
    return q
end


function cell_aspect_ratio(c::AbstractCell)::Float64
    edges = getedges(c)
    L = [ cell_extent(e) for e in edges ]
    return maximum(L)/minimum(L)
end


# Get an array with shares for all nodes
function get_patches(cells::Array{<:AbstractCell,1})
    # get all nodes from cells if needed
    pointsd = Dict{UInt64, Node}()
    for cell in cells
        for node in cell.nodes
            pointsd[hash(node)] = node
        end
    end

    nodes = collect(values(pointsd))
    np     = length(nodes)

    # backup nodes ids
    bk_pt_id = [ pt.id for pt in nodes ]
    for i in 1:np
        nodes[i].id = i
    end

    # get incidence array
    patches  = [ AbstractCell[] for i in 1:np ]
    for cell in cells
        for pt in cell.nodes
            push!(patches[pt.id], cell)
        end
    end

    # restore nodes ids
    for i in 1:np
        nodes[i].id = bk_pt_id[i]
    end

    return nodes, patches
end


function inverse_map(cell::AbstractCell, X::AbstractArray{Float64,1}, tol=1.0e-7)
    return inverse_map(cell.shape, getcoords(cell), X, tol)
end


function get_point(s::Float64, coords::Array{Float64,2})
    #  Interpolates coordinates for s between 0 and 1
    #
    #  0               +1  -->s
    #  1---2---3---..---n  -->idx

    @assert 0<=s<=1
    n = size(coords,1)
    m = n - 1 # number of segments
    δ = 1/m   # parametric length per segment
    
    # segment index and parametric coordinates
    i  = floor(Int, s/δ) + 1 # index for current segment
    i  = min(i, m)
    si = (i-1)*δ     # global parametric coordinate at the beginning of segment i
    t  = (s - si)/δ  # local parametric coordiante

    return coords[i,:]*(1-t) + coords[i+1,:]*t
end

export select


function select(cells::Array{Cell,1}, polycoords::Array{Float64,2}, axis=[0.0, 0.0, 1.0])
    # Selects elements included in the projection of the give polygon and axis direction

    n = size(polycoords, 1)
    eps = 1e-8
    
    # project to plane
    V1 = polycoords[2,:] - polycoords[1,:]
    V2 = cross(axis, V1)
    V1 = cross(V2, axis)
    normalize!(V1)
    normalize!(V2)
    N = normalize(axis)

    R = [V1'; V2'; N']
    XY = (polycoords*R')[:,1:2]
    XY = round.(XY, digits=8)

    selected = Cell[]

    for cell in cells
        allin = true

        for node in cell.nodes
            X = R*node.coord
            x, y = X[1], X[2]
        
            # find if point is inside
            ints = 0

            for i in 1:n
                x1, y1 = XY[i,1], XY[i,2]
        
                # check if point is equal to vertex
                if abs(x1-x)<eps && abs(y1-y)<eps
                    ints = 1
                    break
                end

                # get second point of segment
                if i!=n
                    x2, y2 = XY[i+1,1], XY[i+1,2]
                else
                    x2, y2 = XY[1,1], XY[1,2]
                end

                if y1==y2
                    continue
                end

                if y>=min(y1,y2) && y<=max(y1,y2)
                    # check if point is contained in line
                    if abs((x2-x1)/(y2-y1) - (x2-x)/(y2-y)) < eps
                        ints = 1
                        break
                    end
                    
                    xi = x1 + (x2-x1)/(y2-y1)*(y-y1)
        
                    if xi > x 
                        ints += 1
                    end
                end
        
            end

            if ints%2==0
                allin = false
                break
            end
        end

        allin && push!(selected, cell)
    end

    return selected

end
