# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

#abstract type AbstractCell end

abstract type AbstractBlock<:AbstractCell end


# Cell
# ====

"""
A geometry type that represent a definite line, area or volume.
"""
mutable struct Cell<:AbstractCell
    id     ::Integer
    shape  ::ShapeType
    nodes  ::Array{Node,1}
    tag    ::String
    quality::Float64              # quality index: surf/(reg_surf)
    embedded::Bool                # flag for embedded cells
    crossed::Bool                 # flag if cell crossed by linear inclusion
    oelem  ::Union{AbstractCell,Nothing}  # owner cell if this cell is a face/edge
    linked_elems::Array{AbstractCell,1}   # neighbor cells in case of joint cell
    function Cell(shape::ShapeType, nodes::Array{Node,1}; tag::String="", oelem=nothing)
        this = new()
        this.id = -1
        this.shape = shape
        this.nodes = nodes
        this.tag = tag
        this.quality = 0.0
        this.embedded= false
        this.crossed = false
        this.oelem   = oelem
        this.linked_elems = []
        return this
    end
end

const Face=Cell
const Edge=Cell
const Facet=Cell


### Cell methods

Base.hash(c::Cell) = sum(hash(p) for p in c.nodes)

function get_coords(c::AbstractCell, ndim=3)::Array{Float64,2}
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
function get_nodes(cells::Array{<:AbstractCell,1})::Array{Node,1}
    nodes = Set{Node}()
    for cell in cells
        for node in cell.nodes
            push!(nodes, node)
        end
    end
    R = Node[node for node in nodes]

    return R
end


function Base.getproperty(c::AbstractCell, s::Symbol)
    s == :coords && return get_coords(c)
    s == :faces  && return get_faces(c)
    s == :edges  && return get_edges(c)
    s == :extent && return cell_extent(c)
    s == :points && return c.nodes
    return getfield(c, s)
end

function Base.getproperty(cells::Array{<:AbstractCell,1}, s::Symbol)
    s == :solids   && return filter(cell -> cell.shape.family==SOLID_SHAPE, cells)
    s == :lines    && return filter(cell -> cell.shape.family==LINE_SHAPE, cells)
    s == :joints   && return filter(cell -> cell.shape.family==JOINT_SHAPE, cells)
    s in (:joints1D, :joints1d) && return filter(cell -> cell.shape.family==JOINT1D_SHAPE, cells)
    s == :nodes && return get_nodes(cells)

    error("type $(typeof(cells)) has no property $s")
end


# Index operator for a collection of elements. This function is not type stable
function Base.getindex(cells::Array{<:AbstractCell,1}, s::Symbol)
    s == :all      && return cells
    s == :solids   && return filter(cell -> cell.shape.family==SOLID_SHAPE, cells)
    s == :lines    && return filter(cell -> cell.shape.family==LINE_SHAPE, cells)
    s == :joints   && return filter(cell -> cell.shape.family==JOINT_SHAPE, cells)
    s == :joints1D && return filter(cell -> cell.shape.family==JOINT1D_SHAPE, cells)
    s == :nodes   && return get_nodes(cells)
    error("getindex: $(typeof(cells)) has no property $s")
end


function Base.getindex(cells::Array{Cell,1}, tag::String)
    return filter(cell -> cell.tag==tag, cells)
end


function Base.getindex(cells::Array{Cell,1}, filter_ex::Expr)
    length(cells)==0 && return Cell[]

    # cells id must be set
    nodes = unique(p->p.id, p for c in cells for p in c.nodes )
    sort!(nodes, by=p->p.coord.x+p.coord.y+p.coord.z)
    pointmap = zeros(Int, maximum(node.id for node in nodes) )
    # points and pointmap may have different sizes

    T = Bool[]
    for (i,node) in enumerate(nodes)
        pointmap[node.id] = i
        x, y, z = node.coord.x, node.coord.y, node.coord.z
        push!(T, eval_arith_expr(filter_ex, x=x, y=y, z=z))
    end

    R = Cell[]
    for cell in cells
        all( T[pointmap[node.id]] for node in cell.nodes ) && push!(R, cell)
    end
    return R
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
function get_faces(cell::AbstractCell)
    faces  = Cell[]

    all_faces_idxs = cell.shape.facet_idxs
    facet_shape    = cell.shape.facet_shape

    if facet_shape==() return faces end

    sameshape  = typeof(facet_shape) == ShapeType # check if all facets have the same shape

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_faces_idxs)
        nodes = cell.nodes[face_idxs]
        shape  = sameshape ? facet_shape : facet_shape[i]
        face   = Cell(shape, nodes, tag=cell.tag, oelem=cell)
        push!(faces, face)
    end

    return faces
end


# gets all edges of a cell
function get_edges(cell::AbstractCell)
    if cell.shape.ndim==2 return get_faces(cell) end

    edges  = Cell[]
    all_edge_idxs = cell.shape.edge_idxs

    for edge_idx in all_edge_idxs
        nodes = cell.nodes[edge_idx]
        shape  = (LIN2, LIN3, LIN4)[length(nodes)-1]
        edge   = Cell(shape, nodes, tag=cell.tag, oelem=cell)
        push!(edges, edge)
    end

    return edges
end


# Pseudo determinant of non-square matrices
function norm2(J::Array{Float64,2})

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    if r==1; return norm(J) end
    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
    end
    if r==c; return det(J) end
    error("No rule to calculate norm2 of a $r x $c matrix")
end


# Returns the volume/area/length of a cell
function cell_extent(c::AbstractCell)
    IP = get_ip_coords(c.shape)
    nip = size(IP,1)
    nldim = c.shape.ndim # cell basic dimension

    # get coordinates matrix
    C =get_coords(c)
    J = Array{Float64}(undef, nldim, size(C,2))

    # calc metric
    vol = 0.0
    for i=1:nip
        R    = vec(IP[i,1:3])
        dNdR = c.shape.deriv(R)

        @gemm J = dNdR*C
        w    = IP[i,4]
        normJ = norm2(J)
        #if normJ<0
            #error("cell_extent: Negative Jacobian while calculating cell volume/area/length id=$(c.id) shape=$(c.shape.name) ")
        #end
        vol += normJ*w
    end
    return vol
end


# Returns the surface/perimeter of a regular element given the volume/area of a cell
function regular_surface(metric::Float64, shape::ShapeType)
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
function regular_volume(metric::Float64, shape::ShapeType)
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


#= Returns the cell aspect ratio
function cell_aspect_ratio(c::AbstractCell)
    # get faces
    faces = get_faces(c)
    if length(faces)==0
        return 1.0
    end

    # cell surface
    fmetrics = [ cell_extent(f) for f in faces ]
    surf = sum(fmetrics)
    ar = minimum(fmetrics)/maximum(fmetrics)

    # quality calculation
    metric = cell_extent(c) # volume or area
    rsurf  = regular_surface(metric, c.shape)
    return rsurf/surf
end =#



# Returns the cell quality ratio as vol/reg_vol
function cell_quality_2(c::AbstractCell)::Float64
    # get faces
    faces = get_faces(c)
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
    faces = get_faces(c)
    if length(faces)==0
        return 1.0
    end

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
    faces = get_faces(c)
    len = [ cell_extent(f) for f in faces ]
    c.quality = minimum(len)/maximum(len)
    return c.quality
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
    for i=1:np
        nodes[i].id = i
    end

    # get incidence array
    patches  = [ AbstractCell[] for i=1:np ]
    for cell in cells
        for pt in cell.nodes
            push!(patches[pt.id], cell)
        end
    end

    # restore nodes ids
    for i=1:np
        nodes[i].id = bk_pt_id[i]
    end

    return nodes, patches
end


function inverse_map(cell::AbstractCell, X::AbstractArray{Float64,1}, tol=1.0e-7)
    return inverse_map(cell.shape,get_coords(cell), X, tol)
end

