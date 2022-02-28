# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Cell
# ====

"""
    Cell

A type that represents a mesh cell. It can be a 1D, 2D or 3D.

# Fields
$(FIELDS)
"""
mutable struct Cell<:AbstractCell
    "identification number"
    id     ::Integer
    "geometric shape"
    shape  ::CellShape
    "array of nodes"
    nodes  ::Array{Node,1}
    "string tag used to group cells"
    tag    ::String
    "cell quality in the range from 0 to 1"
    quality::Float64              # quality index: surf/(reg_surf)
    "defines if it is an embedded cell"
    embedded::Bool                # flag for embedded cells
    "defines if it is a crossed cell"
    crossed::Bool                 # flag if cell crossed by linear inclusion
    "owner cell if the cell is a face or edge"
    owner  ::Union{AbstractCell,Nothing}  # owner cell if this cell is a face/edge
    "array of coupled cells"
    linked_elems::Array{AbstractCell,1}   # neighbor cells in case of joint cell

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
    function Cell(shape::CellShape, nodes::Array{Node,1}; tag::String="", owner=nothing, id::Int=-1)
        this = new()
        this.id = id
        this.shape = shape
        this.nodes = copy(nodes)
        this.tag = tag
        this.quality = 0.0
        this.embedded= false
        this.crossed = false
        this.owner   = owner
        this.linked_elems = []
        return this
    end
end

const Face=Cell
const Edge=Cell
const Facet=Cell


### Cell methods

Base.hash(c::Cell) = sum(hash(p) for p in c.nodes)

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
function getcoords(c::AbstractCell, ndim=3)::Array{Float64,2}
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
function getnodes(cells::Array{<:AbstractCell,1})::Array{Node,1}
    return collect(Set(node for cell in cells for node in cell.nodes)) 
end


function (cells::Array{<:AbstractCell,1})(filter=nothing; tag=nothing, shape=nothing, family=nothing, plane=nothing, invert=false)

    filtered = collect(1:length(cells))
    
    if isa(filter, Expr) && filter!=:()
        indexes = copy(filtered)
        filtered = Int[]
        # cells id must be set
        nodes = getnodes(cells)
        pointmap = zeros(Int, maximum(node.id for node in nodes) ) # points and pointmap may have different sizes

        T = Bool[]
        for (i,node) in enumerate(nodes)
            pointmap[node.id] = i
            x, y, z = node.coord.x, node.coord.y, node.coord.z
            push!(T, eval_arith_expr(filter, x=x, y=y, z=z))
        end
 
        filtered = Int[ i for i in indexes if all( T[pointmap[node.id]] for node in cells[i].nodes ) ]
    end

    if isa(filter, Symbol)
        families = Dict(:solids=>SOLID_CELL, :lines=>LINE_CELL, :joints=>JOINT_CELL, :linejoints=>LINEJOINT_CELL, :tipjoints=>TIPJOINT_CELL)
        family = families[filter]
    end

    if tag!==nothing
        indexes = copy(filtered)
        filtered = Int[ i for i in indexes if cells[i].tag==tag ]
    end
    
    if shape!==nothing
        indexes = copy(filtered)
        filtered = Int[ i for i in indexes if cells[i].shape==shape ]
    end

    if family!==nothing
        indexes = copy(filtered)
        filtered = Int[ i for i in indexes if cells[i].shape.family==family ]
    end

    if plane!==nothing
        indexes = copy(filtered)
        filtered = Int[]

        # plane normal
        tol = 1e-5
        A = plane # coplanar points
        I = ones(size(A,1))
        n = normalize(pinv(A.+tol/100)*I) # best fit normal        

        for i in indexes
            cell = cells[i]
            cell.shape.family == JOINT_CELL || continue
            ndim = cell.shape.ndim+1
            @assert ndim==size(A,2)
            C = getcoords(cells[i], ndim)
            I = ones(size(C,1))
            N = pinv(C.+tol/100)*I # best fit normal 
            norm(C*N-I)<tol || continue # check if cell is coplanar
            normalize!(N)
            norm(N-n)<tol || continue # check if planes are parallel
            dot(A[1,:]-C[1,:], n)<tol && push!(filtered, i)
        end
    end

    if invert
        filtered = setdiff(1:length(cells), filtered)
    end

    return cells[filtered]
end

function Base.getproperty(c::AbstractCell, s::Symbol)
    s == :coords && return getcoords(c)
    s == :faces  && return getfaces(c)
    s == :edges  && return getedges(c)
    s == :extent && return cell_extent(c)
    s == :points && return c.nodes
    return getfield(c, s)
end

function Base.getproperty(cells::Array{<:AbstractCell,1}, s::Symbol)
    s == :solids   && return filter(cell -> cell.shape.family==SOLID_CELL, cells)
    s == :lines    && return filter(cell -> cell.shape.family==LINE_CELL, cells)
    s == :joints   && return filter(cell -> cell.shape.family==JOINT_CELL, cells)
    s in (:linejoints, :joints1D, :joints1d) && return filter(cell -> cell.shape.family==LINEJOINT_CELL, cells)
    s == :tipjoints   && return filter(cell -> cell.shape.family==TIPJOINT_CELL, cells)
    s == :embedded && return filter(cell -> cell.shape.family==LINE_CELL && length(cell.linked_elems)>0, cells)
    s == :nodes && return getnodes(cells)
    s == :filter && return cells

    error("type $(typeof(cells)) has no property $s")
end


# Index operator for a collection of elements. This function is not type stable
function Base.getindex(cells::Array{<:AbstractCell,1}, s::Symbol)
    s == :all      && return cells
    s == :solids   && return filter(cell -> cell.shape.family==SOLID_CELL, cells)
    s == :lines    && return filter(cell -> cell.shape.family==LINE_CELL, cells)
    s == :joints   && return filter(cell -> cell.shape.family==JOINT_CELL, cells)
    s in (:linejoints, :joints1D, :joints1d) && return filter(cell -> cell.shape.family==LINEJOINT_CELL, cells)
    s == :tipjoints   && return filter(cell -> cell.shape.family==TIPJOINT_CELL, cells)
    s == :nodes   && return getnodes(cells)
    error("getindex: $(typeof(cells)) has no property $s")
end


function Base.getindex(cells::Array{<:AbstractCell,1}, tag::String)
    return filter(cell -> cell.tag==tag, cells)
end


function Base.getindex(cells::Array{Ty,1}, filter_ex::Expr) where Ty<:AbstractCell
    length(cells)==0 && return Ty[]

    # cells id must be set
    nodes = unique(p->p.id, p for c in cells for p in c.nodes )
    sort!(nodes, by=p->p.coord.x+p.coord.y+p.coord.z)
    pointmap = zeros(Int, maximum(node.id for node in nodes) )
    # points and pointmap may have different sizes

    T = Bool[]
    for (i,node) in enumerate(nodes)
        pointmap[node.id] = i
        x, y, z = node.coord
        push!(T, eval_arith_expr(filter_ex, x=x, y=y, z=z))
    end

    R = Ty[]
    for cell in cells
        all( T[pointmap[node.id]] for node in cell.nodes ) && push!(R, cell)
    end
    return R
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
function getfaces(cell::AbstractCell)
    faces  = Cell[]

    all_faces_idxs = cell.shape.facet_idxs
    facet_shape    = cell.shape.facet_shape

    if facet_shape==() return faces end

    sameshape  = typeof(facet_shape) == CellShape # check if all facets have the same shape

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_faces_idxs)
        nodes = cell.nodes[face_idxs]
        shape  = sameshape ? facet_shape : facet_shape[i]
        face   = Cell(shape, nodes, tag=cell.tag, owner=cell)
        push!(faces, face)
    end

    return faces
end


# gets all edges of a cell
function getedges(cell::AbstractCell)
    if cell.shape.ndim==2 return getfaces(cell) end

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


# # Pseudo determinant of non-square matrices
# function norm2(J::Array{Float64,2})

#     if ndims(J)==1; return norm(J) end

#     r, c = size(J)
#     if r==1; return norm(J) end
#     if r==2 && c==3
#         j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
#         j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
#         j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
#         return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
#     end
#     if r==c; return det(J) end
#     error("No rule to calculate norm2 of a $r x $c matrix")
# end


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
    for i=1:nip
        R    = vec(IP[i,1:3])
        dNdR = c.shape.deriv(R)
        @gemm J = C'*dNdR
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


#= Returns the cell aspect ratio
function cell_aspect_ratio(c::AbstractCell)
    # get faces
    faces = getfaces(c)
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
    faces = getfaces(c)
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
    faces = getfaces(c)
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
    return inverse_map(cell.shape, getcoords(cell), X, tol)
end


function get_point(s::Float64, coords::Array{Float64,2})
    #  Interpolates coordinates for s between 0 and 1
    #
    #  0               +1  -->s
    #  1---2---3---..---n  -->idx

    # @show s
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

function select(cells::Array{Cell,1}, coords::Array{Float64,2}, normal=[0.0, 0.0, 1.0])
    n = size(coords, 1)
    N = normal
    eps = 1e-8
    
    # project to plane
    v1 = normalize(coords[2,:] - coords[1,:])
    v2 = normalize(cross(N, v1))
    v1 = cross(v2, N)

    R = [v1'; v2'; N']
    XY = (coords*R')[:,1:2]
    XY = round.(XY, digits=8)
    # display(XY)

    # @show R

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

                # if node.id==46
                #     @show "hi1"
                #     @show y1
                #     @show y2
                # end
        
                # y1==y2 && continue
                if y1==y2
                    # if y==y1 && x<max(x1,x2)
                        # ints += 1
                    # end
                    continue
                end

                # if node.id==46
                #     @show "hi2"
                # end
        
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

            # if cell.id==19
            #     @show node.coord
            #     @show node.id
            #     @show X
            #     @show ints
            # end

            if ints%2==0
                allin = false
                break
            end
        end

        allin && push!(selected, cell)
    end

    return selected

end

# function select(msh, coords::Array{Float64,2})
#     # Tolerances
#     ε  = 1e-6
#     εn = 1e-4
#     εc = 1e-9
#     # λ  = 1.0

#     itcount = 0

#     # Initial conditions
#     bdist = 0.0  # boundary function initial value
#     len   = 1.0

#     # # Defining required vectors
#     # X1 = vec(coords[  1,:])
#     # Xn = vec(coords[end,:])

#     # Find the initial and final element
    
#     s0 = get_point(εn, coords)
#     ccell = find_elem(s0, msh.elems, msh._elempartition, εc, exclude=Cell[]) # The first tresspased cell
    
#     ccell === nothing && error("Point $(s0) outside the mesh")
    
#     # Initializing more variables
#     ccell = ccell
#     crossedcells = [ ccell ]
#     end_reached  = false
#     s  = 0.0
#     sp = 0.0
#     # nits = round(Int, 1.0/λ)

#     # Splitting inset
#     k = 0
#     while true
#         k +=1
#         ccell_coords = getcoords(ccell)
#         # Default step
#         step  = 0.5*(1.0-s)

#         # Finding step
#         # st = s     # trial point
#         # for i in 1:nits
#         #     st += λ
#         #     if st>1.0; break end
#         #     X = get_point(st, coords, curvetype)
#         #     is_in = is_inside(ccell.shape, ccell_coords, X, ε)
#         #     if !is_in
#         #         step  = 0.5*(st-s)
#         #         break
#         #     end
#         # end



#         s += step
#         @show s

#         X  = get_point(s, coords)
#         @show X
#         @show ccell.id
#         n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

#         itcount+=n ##

#         for i in 1:n

#             step *= 0.5
#             if is_inside(ccell.shape, ccell_coords, X, εc)
#                 s += step
#             else
#                 s -= step
#             end
#             @show s

#             X = get_point(s, coords)

#             R     = inverse_map(ccell.shape, ccell_coords, X)
#             bdist = bdistance(ccell.shape, R)
#         end

#         # Check if end was reached
#         if s > len - εn
#             end_reached = true
#             # X = Xn
#         end

#         end_reached && break

#         # Preparing for the next iteration
#         ncell  = find_elem(get_point(s + εn, coords), msh.elems, msh._elempartition, εc, exclude=[ccell])
#         ncell === nothing && error("Hole found while searching for next crossed cell")

#         ccell = ncell
#         push!(crossedcells, ccell)
#         sp = s
#         s = s + εn
#     end

#     [ c.id for c in crossedcells ]
# end

