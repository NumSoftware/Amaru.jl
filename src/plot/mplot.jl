export mplot, mplotcolorbar

const MOVETO = 1
const LINETO = 2
const CURVE3 = 3
const CURVE4 = 4
const CLOSEPOLY = 79

# function plot_data_for_cell2d(points::Array{Array{Float64,1},1}, shape::CellShape)
function plot_data_for_cell2d(points::Array{Vec3,1}, shape::CellShape)

    points = [ p[1:2] for p in points ]

    if shape==LIN2
        verts = points
        codes = [ MOVETO, LINETO ]
    elseif shape == LIN3
        p1, p2, p3 = points
        cp    = 2*p3 - 0.5*p1 - 0.5*p2
        verts = [ p1, cp, p2 ]
        codes = [ MOVETO, CURVE3, CURVE3]
    elseif shape in (TRI3, QUAD4)
        n = shape==TRI3 ? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i in 1:n
            p2 = i<n ? points[i+1] : points[1]
            push!(verts, p2)
            push!(codes, LINETO)
        end
    elseif shape in (TRI6, QUAD8, QUAD9)
        n = shape==TRI6 ? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i in 1:n
            p1 = points[i]
            p2 = i<n ? points[i+1] : points[1]
            p3 = points[i+n]
            cp = 2*p3 - 0.5*p1 - 0.5*p2
            append!(verts, [cp, p2])
            append!(codes, [CURVE3, CURVE3])
        end
    elseif shape in (QUAD12, QUAD16)
        n = 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i in 1:n
            p1 = points[i]
            p2 = i<n ? points[i+1] : points[1]
            p3 = points[2*i+3]
            p4 = points[2*i+4]
            cp2 = 1/6*(-5*p1+18*p2-9*p3+2*p4)
            cp3 = 1/6*( 2*p1-9*p2+18*p3-5*p4)
            append!(verts, [cp2, cp3, p2])
            append!(codes, [CURVE4, CURVE4, CURVE4])
        end
    elseif shape==JLIN2
        verts = points[1:2]
        codes = [ MOVETO, LINETO ]
    elseif shape == JLIN3
        p1, p2, p3 = points[1:3]
        cp    = 2*p3 - 0.5*p1 - 0.5*p2
        verts = [ p1, cp, p2 ]
        codes = [ MOVETO, CURVE3, CURVE3]
    else
        error("plot_data_for_cell2d: Not implemented for ", shape.name)
    end

    return verts, codes
end

#function plot_data_for_cell3d(points::Array{Array{Float64,1},1}, shape::CellShape)
function plot_data_for_cell3d(points::Array{Vec3,1}, shape::CellShape, V::Vec3=Vec3(1,0,0), width::Float64=0.0)
    if shape == LIN2
        
        W = normalize(cross(V, points[2]-points[1]))
        verts = [
            points[1] - width/2*W,
            points[2] - width/2*W,
            points[2] + width/2*W,
            points[1] + width/2*W,
        ]
    elseif shape == LIN3
        n = 7
        X = [ p[i] for p in points, i in 1:3 ]
        verts = Array{Amaru.Vec3}(undef, 2*n)
        for (i,r) in enumerate(range(-1.0, 1.0, length=n))
            R = [r]
            N = shape.func(R)
            D = shape.deriv(R)
            # J = Vec3(D*X)
            J = Vec3(X'*D)
            W = normalize(cross(V, J))
            Xi = Vec3(N'*X)
            verts[i] = Xi + width/2*W
            verts[2*n-i+1] = Xi - width/2*W
        end
    elseif shape in (TRI3, QUAD4)
        verts = points
    elseif shape == TRI6
        verts = points[[1,4,2,5,3,6]]
    elseif shape in (QUAD8, QUAD9)
        # verts = points[[1,5,2,6,3,7,4,8]]
        n       = 7
        X       = [ p[i] for p in points, i in 1:3 ]
        verts   = Vec3[]
        corners = [[-1.0,-1.0], [1.0,-1.0], [1.0,1.0], [-1.0,1.0], [-1.0,-1.0]]
        
        for k in 1:4
            for r in range(0.0, 1.0, length=n)[1:n-1]
                R = corners[k] + r*(corners[k+1]-corners[k])
                N = shape.func(R)
                push!(verts, Vec3(N'*X))
            end
        end
    else
        error("plot_data_for_cell3d: Not implemented for ", shape.name)
    end
    return verts
end

function plot_data_for_marker3d(point::Vec3, V2, V3, d::Float64)
    n = 12
    verts = Vec3[]
    for θ in range(0, 2*π, length=n)
        vert = point + d/2*cos(θ)*V2 + d/2*sin(θ)*V3
        push!(verts, vert)
    end
    return verts
end


"""
    mplot(blocks, filename="", kwargs...)

Plots an array of blocks using `PyPlot` backend

# Arguments

`blocks` : An array of `Block` objects. Subarrays are also supported.

`filename` = ""` : If provided, a file with the output is saved

# See also

See documentation of `mplot(mesh, filename="", kwargs...)` for details about keyword arguments.
"""
function mplot(items::Union{Block, Array}, filename::String=""; args...)
    # Get list of blocks and check type
    blocks = unfold(items)

    for item in blocks
        isa(item, Block) || error("mplot: Block object expected")
    end

    # Using Nodes and Cell types
    nodes = Array{Node,1}()
    cells  = Array{Cell,1}()

    for bl in blocks
        append!(nodes, bl.nodes)

        if bl.shape.family==BULKCELL
            cell = Cell(bl.shape, bl.nodes)
            push!(cells, cell)
        elseif bl.shape.family==LINECELL
            lines = [ Cell(LIN2, bl.nodes[i-1:i]) for i in 2:length(bl.nodes)]
            append!(cells, lines)
        else
            continue
        end

    end

    mesh = Mesh(cells)
    mplot(mesh, filename; args...)
end


function get_main_edges(cells::Array{<:AbstractCell,1}, angle=120)
    edge_dict  = Dict{UInt64,Cell}()
    faces_dict = Dict{UInt64,Int}( hash(f)=>i for (i,f) in enumerate(cells) )
    main_edges = Cell[]
    # Get faces normals
    normals = [ get_facet_normal(f) for f in cells ]

    # Get edges with non-coplanar adjacent faces
    for face in cells
        face.shape.family == BULKCELL || continue # only surface cells
        face_idx = faces_dict[hash(face)]
        for edge in getedges(face)
            hs = hash(edge)
            edge0 = get(edge_dict, hs, nothing)
            if edge0===nothing
                edge_dict[hs] = edge
            else
                delete!(edge_dict, hs)
                n1 = normals[face_idx] # normal from face
                face0_idx = faces_dict[hash(edge0.owner)]
                n2 = normals[face0_idx] # normal from edge0's parent
                α = 180 - acos( abs(clamp(dot(n1,n2),-1,1)) )*180/pi
                α = round(α, digits=2)
                α<=angle && push!(main_edges, edge)
            end
        end
    end

    return main_edges
end


# import PyCall: PyObject, pyimport, @pydef # required
import PyCall: pyimport
import PyPlot: plt, matplotlib, figure, art3D, Axes3D, ColorMap, gcf

"""
    mplot(mesh, filename="", kwargs...)

Plots a `mesh` using `PyPlot` backend.

# Arguments

`mesh` : A finite element mesh

`filename = ""` : If provided, a file with the output is saved

# Keyword arguments

`axis          = true` : If true, show axes

`lw            = 0.5` : Line width

`markers  = false` : If true, shows node markers

`nodelabels   = false` : If true, shows node labels

`celllabels    = false` : If true, shows cell labels

`opacity       = 1.0`   : Opacity,

`field         = nothing` : If provided, plots corresponding field

`fieldmult    = 1.0` : Factor multiplied to `field` values

`fieldlims     = ()` : Tuple `(min, max)` with field limits

`vectorfield   = nothing` : If provided, plots corresponding vector field

`arrowscale    = 0.0` : Factor multiplied to `vectorfield` values

`colormap      = "coolwarm"` : Colormap according to PyPlot

`colormaplims  = (0.0, 1.0)` : Colormap range to be used

`shrinkcolors  = false` : If true, shrinks the color scale of the colormap

`darkcolors    = false` : If true, makes colormap colors darker

`lightcolors   = false` : If true, makes colormap colors lighter

`vividcolors   = false` : If true, makes colormap colors more vivid

`divergingcolors = false` : If true, makes colormap centralized at zero

`colorbarscale = 0.9` : Scale of the colorbar

`colorbarlabel = ""` : Label of the colorbar

`colorbarlocation = ""` : Location of colorbar (top, bottom, left and right)

`colorbarpad   = 0.0` : Separation of colorbar from the plot

`warpscale     = 0.0` : Factor multiplied to "U" field when available

`hicells       = 0` : Cell number to be highlighted

`elev          = 30.0` : 3D plot elevation

`azim          = 45.0` : 3D plot azimute

`dist          = 10.0` : 3D plot distance from observer

`outline       = true` : Highlight main edges of 3D meshes in the saved output

`figsize       = (3,3.0)` : Figure size

`leaveopen     = false` : If true, leaves the plot open so other drawings can be added
"""
function mplot(
               mesh    ::AbstractDomain,
               filename::String    = "";
               axis                = false,
               lw                  = 0.4,
               rodlw               = 1.0,
               rodcolor            = "indianred",
               facecolor           = "aliceblue",
               markers             = false,
               markerscolor        = "black",
               ms                  = 1.5,
               rodmarkers          = false,
               rodms               = 1.5,
               nodelabels          = false,
               celllabels          = false,
               field               = nothing,
               fieldmult           = 1.0,
               fieldlims           = (),
               vectorfield         = nothing,
               arrowscale          = 0.0,
               shrink              = 1.0,
               opacity             = 1.0,
               lightvector         = nothing,
               colormap            = "coolwarm",
               colormaplims        = (0.0,1.0),
               shrinkcolors        = false,
               darkcolors          = false,
               ligth               = true,
               lightcolors         = false,
               vividcolors         = false,
               divergingcolors     = false,
               shiftcolors         = nothing,
               colorbar            = true,
               colorbarscale       = 0.9,
               colorbarlabel       = "",
               colorbarlocation    = "right",
               colorbarorientation = "vertical",
               colorbarpad         = 0.0,
               colorbarmin         = false,
               colorbarmax         = false,
               warpscale           = 0.0,
               hicells             = 0,
               hicolor             = "ivory",
               elev                = 30.0,
               azim                = 45.0,
               dist                = 10.0,
               outline             = true,
               outlineangle        = 100,
               figsize             = (3,3.0),
               leaveopen           = false,
               crop                = false,
               quiet               = false,
               copypath            = ""
              )

    quiet || headline("Mesh plotting")
    quiet || isempty(filename) && message("generating plot to file $filename")

    if !quiet
        hint("Optional arguments:", level=2)
        options = "axis, lw, markers, nodelabels, celllabels, opacity, field,
                   fieldmult, fieldlims, vectorfield, arrowscale, colormap, colorbarscale,
                   colorbarlabel, colorbarlocation, colorbarorientation, colorbarpad, 
                   warpscale, hicells, elev, azim, dist, outline, outlineangle,
                   figsize, leaveopen, quiet"
        hint(options, level=3)
        hint("Available node fields:", level=2)
        hint(join(keys(mesh.node_data), ", "), level=3)
        hint("Available element fields:", level=2)
        hint(join(keys(mesh.elem_data), ", "), level=3)
    end

    if hicells isa Int && hicells!=0
        hicells  = [ hicells ]
    end

    field !== nothing && (field=string(field))

    if length(mesh.elems.active)==0
        error("mplot: No active elements")
    end

    mmesh = copy(mesh) # copy of original mmesh

    if shrink < 1.0
        nodes = Node[]
        node_ids = Int[]

        # Detach bulk elements by duplicating nodes
        for cell in mmesh.elems
            for (i,node) in enumerate(cell.nodes)
                push!(node_ids, node.id)
                if cell.shape.family == BULKCELL
                    newnode = Node(node.coord)
                    cell.nodes[i] = newnode
                    push!(nodes, newnode)
                else
                    push!(nodes, node)
                end
            end
        end
        mmesh.nodes = nodes

        # update nodal data
        for (field, data) in mmesh.node_data
            mmesh.node_data[field] = data[node_ids,:]
        end

        # update node ids
        for (i,node) in enumerate(mmesh.nodes)
            node.id = i
        end
    end

    # Get initial info from mmesh
    ndim = mmesh.env.ndim
    if ndim==2
        areacells = [ elem for elem in mmesh.elems.active if elem.shape.family==BULKCELL ]
        linecells = [ cell for cell in mmesh.elems.active if cell.shape.family==LINECELL]

        mmesh.elems = [ areacells; linecells ]
        cl_ids = [ [c.id for c in areacells]; [c.id for c in linecells] ]

        mmesh.nodes = getnodes(mmesh.elems.active)
        pt_ids = [ p.id for p in mmesh.nodes ]

    elseif ndim==3
        # get surface cells and update
        volcells  = [ elem for elem in mmesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==3 ]
        areacells = [ elem for elem in mmesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==2 ]
        surfcells = get_surface(volcells)
        linecells = [ cell for cell in mmesh.elems.active if cell.shape.family==LINECELL]
        outlinecells = outline ? get_outline_edges(surfcells) : Cell[]

        mmesh.elems = [ surfcells; areacells; linecells ]
        cl_ids = [ [c.owner.id for c in surfcells]; [c.id for c in linecells]; [c.id for c in areacells] ]

        mmesh.nodes = getnodes(mmesh.elems.active)
        pt_ids = [ p.id for p in mmesh.nodes ]

        # observer and light vectors
        V = Vec3( cosd(elev)*cosd(azim), cosd(elev)*sind(azim), sind(elev) )

        lightvector===nothing && (lightvector=V) 
        if lightvector isa AbstractArray
            L = lightvector
        else
            error("mplot: lightvector must be a vector.")
        end
    end

    # update data
    for (field, data) in mmesh.node_data
        mmesh.node_data[field] = data[pt_ids,:]
    end
    for (field, data) in mmesh.elem_data
        mmesh.elem_data[field] = data[cl_ids]
    end

    # update node ids
    for (i,node) in enumerate(mmesh.nodes)
        node.id = i
    end
    
    # update cell ids
    for (i,cell) in enumerate(mmesh.elems)
        cell.id = i
    end
    
    connect = [ [ node.id for node in cell.nodes ] for cell in mmesh.elems  ]

    ncells = length(mmesh.elems)
    nnodes = length(mmesh.nodes)
    node_data = mmesh.node_data
    elem_data = mmesh.elem_data
  
    # Change coords if warping
    if warpscale>0.0
        if haskey(node_data, "U")
            U = node_data["U"]
            for (i,node) in enumerate(mmesh.nodes)
                node.coord = node.coord + warpscale*U[i,:]  
            end
        else
            alert("mplot: Vector field U not found for warping.")
        end
    end

    # Change coords if shrink
    if shrink<1.0
        for cell in mmesh.elems
            shape = cell.shape
            if shape.family==BULKCELL
                X = getcoords(cell)
                C = vec(mean(X, dims=1))
                for node in cell.nodes
                    node.coord = (node.coord - C)*shrink + C
                end
            end
        end
    end
    

    # Data limits
    limX = collect(extrema( node.coord[1] for node in mmesh.nodes) )
    limY = collect(extrema( node.coord[2] for node in mmesh.nodes) )
    limZ = collect(extrema( node.coord[3] for node in mmesh.nodes) )

    ll = max(diff(limX)[1], diff(limY)[1], diff(limZ)[1])

    isinter = plt.isinteractive()
    filename!="" && plt.ioff()

    if facecolor isa String
        facecolor = matplotlib.colors.to_rgba(facecolor)
    end

    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=6)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", lw=0.5)
    plt.rc("legend", fontsize=6)
    plt.rc("figure", figsize=figsize) # suggested size (4.5,3)

    # Configure plot
    if ndim==3
        ax = Axes3D(figure())

        # Set limits
        meanX = mean(limX)
        meanY = mean(limY)
        meanZ = mean(limZ)
        limX = [meanX-ll/2, meanX+ll/2]
        limY = [meanY-ll/2, meanY+ll/2]
        limZ = [meanZ-ll/2, meanZ+ll/2]
        ax.set_xlim( meanX-ll/2, meanX+ll/2)
        ax.set_ylim( meanY-ll/2, meanY+ll/2)
        ax.set_zlim( meanZ-ll/2, meanZ+ll/2)
        ax.set_box_aspect((diff(limX)[1], diff(limY)[1], diff(limZ)[1]))  # instead of ax.set_aspect("equal")

        # Labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        axis == false && plt.axis("off")
    else
        ax = plt.axes()
        ax.set_aspect("equal", "datalim")

        # Set limits
        limX .+= [ -0.01*ll, 0.01*ll]
        limY .+= [ -0.01*ll, 0.01*ll]
        ax.set_xlim(limX...)
        ax.set_ylim(limY...)

        # Labels
        ax.set_xlabel.("x")
        ax.set_ylabel.("y")
        axis == false && plt.axis("off")
    end

    has_field = field !== nothing
    if has_field
        colorbarlabel = colorbarlabel=="" ? field : colorbarlabel
        field = string(field)
        found = haskey(elem_data, field)
        if found
            fvals = elem_data[field]
        else
            found = haskey(node_data, field)
            found || error("mplot: field $field not found")
            data  = node_data[field]
            fvals = [ mean(data[connect[i]]) for i in 1:ncells ]
        end
        fvals *= fieldmult
        fieldlims==() && (fieldlims = extrema(fvals))

        if !(fieldlims[1]<0 && fieldlims[2]>0)
            divergingcolors = false
        end

        if colormap isa String
            # colormap may be "coolwarm", "bone", "plasma", "inferno", etc.
            colormaps = matplotlib.pyplot.colormaps()
            if colormap in colormaps
                cmap = matplotlib.cm.get_cmap(colormap)
            else
                error("mplot: Invalid colormap $colormap \n", 
                      "colormap should be one of:\n", colormaps)
            end
        elseif !(colormap isa ColorMap)
            cdict = Dict("red"   => [(0.0,  0.8, 0.8), (0.5, 0.7, 0.7), (1.0, 0.0, 0.0)],
                         "green" => [(0.0,  0.2, 0.2), (0.5, 0.7, 0.7), (1.0, 0.2, 0.2)],
                         "blue"  => [(0.0,  0.0, 0.0), (0.5, 0.7, 0.7), (1.0, 0.6, 0.6)])

            cmap = matplotlib.colors.LinearSegmentedColormap("custom_colormap", cdict, 256)
        else
            cmap = colormap
        end

        if colormaplims!=(0.0,1.0)
            allcolors = [ cmap(p) for p in range(colormaplims[1], colormaplims[2], length=21) ]
            Q = range(0,1,length=21)

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
        end

        if shrinkcolors
            allcolors = [ cmap(p) for p in range(0,1,length=21) ]
            Q = [ 0.5+0.5*(abs((p-0.5)/0.5))^3.0*sign(p-0.5) for p in range(0,1,length=21) ]

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
        end

        if shiftcolors!==nothing
            src, tgt = shiftcolors
            allcolors = [ cmap(p) for p in range(0,1,length=21) ]
            Q = [ ( p<src ? p*tgt/src : tgt + (p-src)*(1-tgt)/(1-src) ) |> x->round(x,digits=3)  for p in range(0,1,length=21) ]

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
        end

        if darkcolors || lightcolors || vividcolors
            if darkcolors
                allcolors = [ 0.9.*cmap(p) for p in range(0,1,length=21) ]
            elseif lightcolors
                function lighten(c)
                    return c .+ 0.1.*(1.0.-c)
                end
                allcolors = [ lighten(cmap(p)) for p in range(0,1,length=21) ]
            else
                colorsys = pyimport("colorsys")
                function vivid(c)
                    h, l, s = colorsys.rgb_to_hls(c[1], c[2], c[3])
                    return colorsys.hls_to_rgb(h, l*0.9, s + 0.1*(1-s))
                end
                allcolors = [ vivid(cmap(p)) for p in range(0,1,length=21) ]
            end
            Q = range(0,1,length=21)

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
        end

        if divergingcolors
            q0 = -fieldlims[1]/(fieldlims[2]-fieldlims[1])
            # function to recalculate colors positions
            function q(p::Float64) 
                if p<0.5
                    return 2*p*q0
                else
                    return q0 + 2*(p-0.5)*(1-q0)
                end
            end

            P = range(0,1,length=21)
            allcolors = [ cmap(p) for p in P ]
            Q = q.(P)

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
        end
    end

    # Check for line field
    has_line_field = false
    if has_field
        for (i,cell) in enumerate(mmesh.elems)
            cell.shape.family == LINECELL || continue
            if fvals[i]!=0.0
                has_line_field = true
                break
            end
        end
    end

    # Plot cells
    if ndim==3
        # Plot cells
        all_verts = []
        edgecolor = []
        facecolors = []
        lineweight = []

        for (i,cell) in enumerate(mmesh.elems)
            shape = cell.shape
            shape.family in (BULKCELL, LINECELL) || continue
            
            points = [ node.coord for node in cell.nodes ]

            if shape.family==BULKCELL
                verts = plot_data_for_cell3d(points, shape)

                has_line_field && continue

                if !has_field || has_line_field
                    # op = has_line_field ? 0.0 : opacity
                    fc = (0.94, 0.97, 1.0, opacity)
                    # ec = (0.4, 0.4, 0.4, 1-0.75*(1-opacity))
                    # ec = (0.4, 0.4, 0.4, opacity)
                    ec = (0.5, 0.5, 0.5, opacity)
                else
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    fc = cmap(v)
                    ec = Tuple( (fc .+ [0.2, 0.2, 0.2, opacity])./2 )
                end

                ew = lw

                N = get_facet_normal(cell)
                R = 2*N*dot(L,N) - L
                
                f = 1.0
                if ligth
                    #f = 0.7+0.3*abs(dot(L,N)*(1-dot(V,R))/2)
                    #f = 0.6+0.3*abs(dot(L,N)) + 0.2*(1-dot(V,R))/2
                    # f = 0.6+0.2*abs(dot(L,N)) + 0.2*(1+dot(V,R))/2
                    f = 0.8 + 0.1*abs(dot(L,N)) + 0.1*(1+dot(V,R))/2
                    f = min(f, 1.0)

                    # f = 0.9+0.05*abs(dot(L,N)) + 0.05*(1+dot(V,R))/2
                end
                fc = (f*fc[1], f*fc[2], f*fc[3], opacity)
            elseif shape.family==LINECELL
                verts = plot_data_for_cell3d(points, shape, V, 0.0075*ll*rodlw)

                if has_line_field
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    ec = cmap(v)
                else
                    ec = rodcolor 
                end
                # fc = (0.0, 0.0, 0.0, 0.0)
                fc = ec
                ew = rodlw  #! lineweight is not working in Poly3DCollection
                ew = lw
            end

            push!(all_verts, verts)
            push!(edgecolor, ec)
            push!(facecolors, fc)
            push!(lineweight, ew)
        end

        for cell in outlinecells
            ΔX = 0.01*ll*V
            points = [ node.coord+ΔX for node in cell.nodes ]
            verts = plot_data_for_cell3d(points, cell.shape, V, 0.001)

            # ec = "black"
            ec = "dimgray"
            fc = (0.0, 0.0, 0.0, 0.0)

            push!(all_verts, verts)
            push!(lineweight, lw*1.1)
            push!(edgecolor, ec)
            push!(facecolors, fc)
        end
        
        if rodmarkers && filename!=""
            V2 = normalize(cross(V, rand(3)))
            V3 = normalize(cross(V, V2))
            ΔX = 0.01*ll*V
            for line in linecells, node in line.nodes
                verts = plot_data_for_marker3d(node.coord+4*ΔX, V2, V3,  0.015*ll*rodlw)
                ec = (0.0, 0.0, 0.0, 0.0)
                fc = "black"
                push!(all_verts, verts)
                push!(edgecolor, ec)
                push!(facecolors, fc)
            end
        end

        cltn = art3D.Poly3DCollection(all_verts, facecolor=facecolors, edgecolor=edgecolor, lw=lineweight, alpha=opacity) # ! lineweight is not working in Poly3DCollection
        ax.add_collection3d(cltn)

        if has_field && colorbar
            normalized = matplotlib.colors.Normalize(fieldlims[1], fieldlims[2])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalized)
            sm.set_array([])

            cbar = plt.colorbar(sm, label=colorbarlabel, shrink=colorbarscale, aspect=10*colorbarscale*figsize[2], 
                                pad=colorbarpad, location=colorbarlocation)
        end


    elseif ndim==2

        all_patches = []
        edgecolor   = []
        facecolors  = []
        lineweight  = []

        
        for (i,cell) in enumerate(mmesh.elems)
            shape = cell.shape
            shape.family in (BULKCELL, LINECELL) || continue
            points = [ node.coord for node in cell.nodes ]
            verts, codes = plot_data_for_cell2d(points, shape)
            path  = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path)


            if shape.family==BULKCELL
                if !has_field || has_line_field
                    fc = facecolor
                    ec = Tuple( (fc .+ [0.3, 0.3, 0.3, 1.0])./2 )

                    if cell.id in hicells      
                        ec = "black"
                        fc = hicolor
                    end
                else
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    fc = cmap(v)
                    ec = Tuple( (fc .+ [0.3, 0.3, 0.3, 1.0])./2 )
                end
                push!(lineweight, lw)
            elseif shape.family==LINECELL
                if has_line_field
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    ec = cmap(v)
                else
                    ec = rodcolor 
                end
                fc = (0.0, 0.0, 0.0, 0.0)
                push!(lineweight, rodlw )
            end
            push!(edgecolor, ec)
            push!(facecolors, fc)
            push!(all_patches, patch)
        end

        cltn = matplotlib.collections.PatchCollection(all_patches, edgecolor=edgecolor, facecolor=facecolors, lw=lineweight)
        ax.add_collection(cltn)
        
        if has_field && colorbar
            normalized = matplotlib.colors.Normalize(fieldlims[1], fieldlims[2])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalized)
            sm.set_array([])

            h = colorbarorientation=="vertical" ? figsize[2] : figsize[1]
            h = norm(figsize)
            cbar = plt.colorbar(sm, label=colorbarlabel, shrink=colorbarscale, aspect=4*colorbarscale*h, 
                                # format="%.2f", 
                                pad=colorbarpad+0.05, orientation=colorbarorientation)
            # cbar.ax.tick_params(labelsize=6)
            # cbar.outline.set_linewidth(0.0)
            # cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
            # cbar.update_ticks()
            # cbar.solids.set_alpha(1)
        end
    end

    if has_field && colorbar
        cbar.ax.tick_params(labelsize=6)
        cbar.outline.set_linewidth(0.0)
        cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
        cbar.update_ticks()
        cbar.solids.set_alpha(1)

        if colorbarmax
            ticks = [ cbar.get_ticks(); fieldlims[2] ]
            # @show ticks
            label = @sprintf("%g", round(fieldlims[2], sigdigits=3))
            if colorbarorientation=="vertical"
                labels = cbar.ax.get_yticklabels()
            else
                labels = cbar.ax.get_xticklabels()
            end
            labels = [ labels; label ]
            # @show labels
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(labels)
        end
        if colorbarmin
            ticks = [ fieldlims[1]; cbar.get_ticks()]
            label = @sprintf("%g", round(fieldlims[1], sigdigits=3))
            if colorbarorientation=="vertical"
                labels = cbar.ax.get_yticklabels()
            else
                labels = cbar.ax.get_xticklabels()
            end
            labels = [ label; labels ]
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(labels)
        end
    end

    # Draw nodes
    if markers
        if ndim==3
            coords = [ node.coord for cell in mmesh.elems.solids for node in cell.nodes ]

            # X = [ node.coord.x for node in mmesh.nodes ]
            # Y = [ node.coord.y for node in mmesh.nodes ]
            # Z = [ node.coord.z for node in mmesh.nodes ]
            X = getindex.(coords,1)
            Y = getindex.(coords,2)
            Z = getindex.(coords,3)
            ax.scatter(X, Y, Z, color=markerscolor, marker="o", s=ms)
        else
            subcoords = [ node.coord for cell in mmesh.elems.solids for node in cell.nodes ]

            X = getindex.(subcoords,1)
            Y = getindex.(subcoords,2)
            plt.plot(X, Y, color="black", marker="o", mfc=markerscolor, fillstyle="full", markersize=ms, mew=0, lw=0, clip_on=false)
        end
    end

   # Node markers on line cells
    if rodmarkers && ndim==2
        subcoords = [ node.coord for cell in mmesh.elems.lines for node in cell.nodes ]
        X = getindex.(subcoords,1)
        Y = getindex.(subcoords,2)

        plt.plot(X, Y, color="black", marker="o", markersize=rodms, lw=0, mew=0, clip_on=false)

        # if ndim==3
        #     Z = [ node.coord.z for line in linecells for node in line.nodes ]
        #     ax.scatter(X, Y, Z, marker="o", s=ms, facecolor="black", edgecolor="none", depthshade=false)
        # else
            # plt.plot(X, Y, color="black", marker="o", markersize=ms, lw=0)
        # end
    end

  #  # Draw arrows
  #  if vectorfield!=nothing && ndim==2
  #      data = node_data[vectorfield]
  #      color = "blue"
  #      if arrowscale==0
  #          plt.quiver(X, Y, data[:,1], data[:,2], color=color)
  #      else
  #          plt.quiver(X, Y, data[:,1], data[:,2], color=color, scale=1.0/arrowscale)
  #      end
  #  end

    # Draw node numbers
    if nodelabels
        # nnodes = length(X)
        nnodes = length(mesh.nodes)
        for i in 1:nnodes
            x, y, z = mesh.nodes[i].coord .+ [0.015, -0.03, -0.02]*ll
            # @show  i
            # @show  mesh.nodes[i].coord 
            # x = nodes[i].coords.x
            # x = nodes[i].coords.x
            # x = X[i] + 0.01*L
            # y = Y[i] - 0.01*L
            # z = Z[i] - 0.01*L
            if ndim==3
                ax.text(x, y, z, mesh.nodes[i].id, va="center", ha="center", backgroundcolor="none")
            else
                ax.text(x, y, mesh.nodes[i].id, va="top", ha="left", backgroundcolor="none")
            end
        end
    end

    # Draw cell numbers
    if celllabels && ndim==2
        ncells = length(mesh.elems)
        for i in 1:ncells
            coo = getcoords(mesh.elems[i])
            x = mean(coo[:,1]) + 0.015*ll
            y = mean(coo[:,2]) - 0.035*ll
            # ax.text(x, y, "($i)", va="top", ha="left", color="blue", backgroundcolor="none")
            ax.text(x, y, "($(mesh.elems[i].id))", va="top", ha="left", backgroundcolor="none")
        end
    end

    if ndim==3
        ax.view_init(elev=elev, azim=azim)
        ax.dist = dist
    end

    if filename==""
        return gcf()
    else
        _, format = splitext(filename)
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.0, format=format[2:end])

        if crop && format==".pdf"
            if Sys.islinux()
                cmd = `pdfcrop $filename $filename`
                out = Pipe()
                err = Pipe()
                run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
                close(out.in)
                close(err.in)
            else
                alert("crop option is not available on $(Sys.KERNEL)")
            end
        end

        quiet || info("file $filename saved")

        if copypath!=""
            if isdir(copypath)
                copyfile = joinpath(copypath, basename(filename))
            else
                copyfile = copypath
            end
            cp(filename, copyfile, force=true)
            quiet || info("file $copyfile saved")
        end

        isinter && plt.ion()
    end

    # Do not close if in IJulia
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return nothing
    end

    leaveopen || plt.close("all")

    return nothing
end


function mplotcolorbar(
    filename::String = "";
    colormap         = "coolwarm",
    fieldlims        = (0.0,1.0),
    nbins            = 0,
    digits           = -5,
    discrete         = false,
    fontsize         = 6.5,
    scale            = 0.9,
    label            = "",
    orientation      = "vertical",
    figsize          = (3,3.0),
    crop             = false,
    aspect           = 20,
    colormaplims      = (0.0,1.0),
    #shrinkcolor       = false,
    #darkcolor         = false,
    #lightcolor        = false,
    #vividcolor        = false,
    #divergingcolor    = false,
)    

    # import PyPlot:plt, matplotlib, figure, art3D, Axes3D, ioff, ColorMap, gca, gcf
    # @eval ioff()

    # fix PyPlot
    # @eval import PyPlot:getproperty, LazyPyModule
    # if ! @eval hasmethod(getproperty, (LazyPyModule, AbstractString))
    #     @eval Base.getproperty(lm::LazyPyModule, s::AbstractString) = getproperty(PyObject(lm), s)
    # end

    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", lw=0.5)
    plt.rc("legend", fontsize=fontsize)
    plt.rc("figure", figsize=figsize) # suggested size (4.5,3)

    fig   = plt.figure()
    scale = 0.5
    axes  = gcf().add_axes([0, 0, scale/aspect, scale])

    cmap = matplotlib.cm.get_cmap(colormap)

    if colormaplims!=(0.0,1.0)
        allcolors = [ cmap(p) for p in range(colormaplims[1], colormaplims[2], length=21) ]
        Q = range(0,1,length=21)

        cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                     "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                     "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

        cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap", cdict, 256)
    end

    if discrete
        cmap = matplotlib.cm.get_cmap(colormap, nbins)
    end

    if nbins>0
        ticks = collect(range(fieldlims[1], fieldlims[2], length=nbins+1))
        if digits>-5
            ticks = round.(ticks, digits=digits)
        end
    end

    format = nothing
    if digits>-5
        format = "%0."*string(digits)*"f"
    end

    cbar = matplotlib.colorbar.ColorbarBase(
        axes,
        label       = label,
        norm        = matplotlib.colors.Normalize(fieldlims[1], fieldlims[2]),
        drawedges   = false,
        orientation = orientation,
        cmap        = cmap,
        format      = format,
        #values=fieldlims,
        #boundaries=fieldlims,
        # ticks=ticks,
    )

    cbar.ax.tick_params(labelsize=fontsize)
    cbar.outline.set_linewidth(0.0)
    # cbar.locator = matplotlib.ticker.MaxNLocator(nbins=nbins)
    # cbar.set_ticks(ticks) # forces ticks
    # cbar.update_ticks()
    cbar.solids.set_alpha(1)

    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight")

        if crop
            if Sys.islinux()
                cmd = `pdfcrop $filename $filename`
                out = Pipe()
                err = Pipe()
                run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
                close(out.in)
                close(err.in)
            else
                alert("crop option is not available on $(Sys.KERNEL)")
            end
        end
    end
end



function round_for_scale(x)
    x = round(x, sigdigits=1)
    ex = floor(log10(x)) # exponent
    n = round(x/10^ex, sigdigits=1) # first significant digit
    if n>=3 
        if n==3
            n=3.0
        elseif n<=7
            n=5
        else
            n=10
        end
    end
    return round(n*10^ex, sigdigits=1)
end

using LaTeXStrings
export mplot_linear

"""
    mplot(elems, filename="", kwargs...)

Plots linear `elems` using `PyPlot` backend. If `filename` is provided it writes a pdf file containing the plot.

# Arguments

`mesh` : A finite element mesh

`filename = ""` : If provided, a pdf file with the output is saved

# Keyword arguments

`axis          = true` : If true, show axes

"""

function mplot(
    elems   ::Array{Element,1},
    filename::String = "";
    field            = nothing,
    fieldmult        = 1.0,
    fieldunits       = "",
    barscale         = 1.0,
    axis             = false,
    xlabel           = "",
    ylabel           = "",
    xlim             = nothing,
    ylim             = nothing,
    legendlabels     = [],
    showscale        = true,
    scalepos         = (0.55, 0.2),
    fullscale        = false,
    fontsize         = 6.5,
    copypath         = "",
    figsize          = (3.0, 3.0),
    minmax           = false
)
    
    # @eval import PyPlot:plt, matplotlib, figure, gca, gcf
    filename != "" && ioff()

    # @assert barscale>0

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("figure", figsize=figsize)
    maxv = 0.0

    lines = [ elem.shape.family==LINECELL ? elem : elem.linked_elems[2] for elem in elems ]
    coords = getcoords(getnodes(lines))
    #sumx = maximum(abs, coords[:,1])
    #sumy = maximum(abs, coords[:,2])
    #sumz = maximum(abs, coords[:,3])

    n = size(coords, 1)
    sumx = sum(coords[:,1])
    sumy = sum(coords[:,2])
    sumz = sum(coords[:,3])
    avgx = sumx/n
    avgy = sumy/n
    avgz = sumz/n
    devx = sum((coords[:,1].-avgx).^2)/n
    devy = sum((coords[:,2].-avgy).^2)/n
    devz = sum((coords[:,3].-avgz).^2)/n
    tol = 1e-10

    if devy<tol
        xidx = 1
        yidx = 3
        xlabel=="" && (xlabel=raw"$x$")
        ylabel=="" && (ylabel=raw"$z$")
    elseif devx<tol
        xidx = 2
        yidx = 3
        xlabel=="" && (xlabel=raw"$y$")
        ylabel=="" && (ylabel=raw"$z$")
    elseif devz<tol
        xidx = 1
        yidx = 2
        xlabel=="" && (xlabel=raw"$x$")
        ylabel=="" && (ylabel=raw"$y$")
    else
        xidx = 1
        yidx = 3
        xlabel=="" && (xlabel=raw"$x$")
        ylabel=="" && (ylabel=raw"$z$")
    end
    
    # Get maximum values
    vmax = 0.0
    vmin = 0.0
    xmin = 0.0
    xmax = 0.0
    for elem in elems
        for ip in elem.ips
            v = ip_state_vals(elem.matparams, ip.state)[Symbol(field)]
            vmin = min(vmin, v)
            vmax = max(vmax, v)
            xmin = min(xmin, ip.coord[xidx])
            xmax = max(xmax, ip.coord[xidx])
        end
    end

    if field === nothing 
        printstyled("mplot_linear:\n", color=:cyan )
    else
        printstyled("mplot_linear: plotting field $field\n", color=:cyan )
        printstyled("  (min,max): ($(vmin*abs(fieldmult)), $(vmax*abs(fieldmult)))\n", color=:light_black)
    end

    maxv = max(abs(vmin), abs(vmax))
    xwidth = 1.3*(xmax-xmin) # estimative of xwidth
    fieldtrans = 0.2*xwidth/maxv
    dn = 0.005*xwidth

    for elem in elems
        line = elem.shape.family==LINECELL ? elem : elem.linked_elems[2]
        X = [ node.coord[xidx] for node in line.nodes ]
        Y = [ node.coord[yidx] for node in line.nodes ]
        plt.plot(X, Y, "tab:gray", lw=2, solid_capstyle="round", zorder=2) # plot line

        field===nothing && continue

        V = [ X[2]-X[1], Y[2]-Y[1] ]
        normalize!(V)
        N = [ -V[2], V[1] ]

        ips = elem.ips
        for ip in ips
            x = ip.coord[xidx]
            y = ip.coord[yidx]
            v = ip_state_vals(elem.matparams, ip.state)[Symbol(field)]*sign(fieldmult)
            X = [ x, x+N[1]*v*fieldtrans*barscale ]
            Y = [ y, y+N[2]*v*fieldtrans*barscale ]
            X .+= dn*sign(v)*N[1]
            Y .+= dn*sign(v)*N[2]

            color = v>0 ? "tab:red" : "tab:blue"
            plt.plot(X,Y, color, ls="-", lw=1, solid_capstyle="round", zorder=1) # normal bar
            # plt.plot(x,y,"k+", ms=3, mew=0.5, zorder=3) # cross symbol
            plt.plot(x,y,"k.", ms=2, mew=0.5, zorder=3) # cross symbol
        end
    end

    ca = @eval gca()

    if showscale
        xmin, xmax = ca.get_xlim()
        xwidth = xmax-xmin # recalculate xwidth according to plot

        refval = round_for_scale(maxv)
        scalelen = refval*fieldtrans*abs(barscale)*1.0/xwidth

        xpos, ypos = scalepos

        scaletitle = "Scale"
        if fullscale
            scaletitle *= fieldunits!="" ? " [$fieldunits]" : "" 
        end
        plt.text(xpos, ypos, scaletitle, transform = ca.transAxes, ha="left", va="top")
        
        ypos -= 0.06
        ypos0 = ypos
        if vmin<0
            plt.plot([xpos, xpos+scalelen], [ypos, ypos], "tab:blue", ls="-", lw=1, solid_capstyle="butt", zorder=1, transform = ca.transAxes)
            ypos -= 0.02
        end
        if vmax>0
            plt.plot([xpos, xpos+scalelen], [ypos, ypos], "tab:red", ls="-", lw=1, solid_capstyle="butt", zorder=1, transform = ca.transAxes)
            ypos -= 0.02
        end
        

        if fullscale
            ticks = (0.0, 0.5, 1.0)
        else
            ticks = (0.0, 1.0)
        end
        for f in ticks
            plt.plot([xpos+f*scalelen, xpos+f*scalelen], [ypos+0.01, ypos], "black", ls="-", lw=0.5, zorder=1, transform = ca.transAxes)
        end

        ypos -= 0.01
        if fullscale
            for f in ticks
                plt.text(xpos+f*scalelen, ypos, "$(f*refval*abs(fieldmult))", transform=ca.transAxes, ha="center", va="top")
            end
        else
            plt.text(xpos+0.5*scalelen, ypos, "$(refval*abs(fieldmult)) $fieldunits", transform=ca.transAxes, ha="center", va="top")
        end
    end

    #ca.set_aspect("equal", "datalim")
    ca.set_aspect("equal", "box")

    # Set axis limits
    xlim!==nothing && plt.xlim(xlim)
    ylim!==nothing && plt.ylim(ylim)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    if length(legendlabels)>0
        lines = []
        idxs = []
        if vmin<0
            push!(lines, matplotlib.lines.Line2D([0],[0], lw=1, color="tab:blue", label="-"))
            push!(idxs, 1)
        end

        if vmax>0
            push!(lines, matplotlib.lines.Line2D([0],[0], lw=1, color="tab:red", label="+"))
            push!(idxs, 2)
        end

        ncol = 1
        if length(legendlabels)==2
            legendlabels = legendlabels[idxs]
            ncol = 2
        end

        legend = plt.legend(lines, legendlabels,
            loc            = "lower left",
            bbox_to_anchor = (-0.02, 1.01, 1.04, 0.2),
            edgecolor      = "k",
            ncol           = ncol,
        )
    end

    if minmax
        vrmin = round(vmin*fieldmult, sigdigits=4)
        vrmax = round(vmax*fieldmult, sigdigits=4)
        if fieldmult<0
            vrmin, vrmax = vrmax, vrmin
        end

        tmin  = "\$$(vrmin>0 ? "+" : "" )$vrmin\$ $fieldunits"
        tmax  = "\$$(vrmax>0 ? "+" : "" )$vrmax\$ $fieldunits"
        plt.text(0.5, 0.98, "min: $tmin   max: $tmax", transform = ca.transAxes, ha="center", va="top")
    end
    
    ax = plt.axes()
    ax.xaxis.set_tick_params(width=0.3)
    ax.yaxis.set_tick_params(width=0.3)
    ax.xaxis.set_tick_params(size=3.0)
    ax.yaxis.set_tick_params(size=3.0)
    #if ticksinside
        ax.tick_params(which="minor", axis="x", direction="in")
        ax.tick_params(which="minor", axis="y", direction="in")
        ax.tick_params(which="major", axis="x", direction="in")
        ax.tick_params(which="major", axis="y", direction="in")
    #end


    # show or save plot
    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.01, format="pdf")
        plt.close("all")
        printstyled("  file $filename saved\n", color=:cyan)
        if copypath!=""
            if isdir(copypath)
                copyfile = joinpath(copypath, basename(filename))
            else
                copyfile = copypath
            end
            cp(filename, copyfile, force=true)
            printstyled("  file $copyfile saved\n", color=:cyan)
        end
    end

    return @eval gcf()
end


function mplot(geo::GeoModel, filename::String;
    axis=false,
    quiet=false,
    azim=45,
    elev=30, 
    dist=10.0,
    crop=false,
    copypath="",
    )
    # check dimension
    ndim = 2
    X = Float64[]
    Y = Float64[]
    Z = Float64[]
    for p in geo.entities
        p isa Point || continue
        push!(X, p.coord[1])
        push!(Y, p.coord[2])
        push!(Z, p.coord[3])
        if p.coord[3] != 0.0
            ndim = 3
        end
    end

    # Data limits
    limX = collect(extrema(X))
    limY = collect(extrema(Y))
    limZ = collect(extrema(Z))
    ll = max(diff(limX)[1], diff(limY)[1], diff(limZ)[1])

    plt.close("all")

    # Configure plot
    if ndim==3
        ax = Axes3D(figure())

        # Set limits
        meanX = mean(limX)
        meanY = mean(limY)
        meanZ = mean(limZ)
        limX = [meanX-ll/2, meanX+ll/2]
        limY = [meanY-ll/2, meanY+ll/2]
        limZ = [meanZ-ll/2, meanZ+ll/2]
        ax.set_xlim( meanX-ll/2, meanX+ll/2)
        ax.set_ylim( meanY-ll/2, meanY+ll/2)
        ax.set_zlim( meanZ-ll/2, meanZ+ll/2)
        ax.set_box_aspect((diff(limX)[1], diff(limY)[1], diff(limZ)[1]))  # instead of ax.set_aspect("equal")

        # Labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        axis == false && plt.axis("off")

        ax.view_init(elev=elev, azim=azim)
        ax.dist = dist
    else
        ax = plt.axes()
        ax.set_aspect("equal", "datalim")

        # Set limits
        ax.set_xlim(limX...)
        ax.set_ylim(limY...)

        # Labels
        ax.set_xlabel.("x")
        ax.set_ylabel.("y")
        axis == false && plt.axis("off")
    end


    # draw surfaces
    if ndim==2 
        # plot all lines
        all_patches = []
        edgecolor   = []
        facecolors  = []
        lineweight  = []
        for s in geo.entities
            s isa Surface || continue
            points = getpoints(s.loops[1])
            npoints = length(points)
            # verts = [ p.coord[1:2] for p in points ]
            # push!(verts, verts[1])
            # verts = [ verts; [ verts[1] ] ]
            # codes = [ MOVETO; repeat([LINETO], npoints) ]
            # @show verts
            # @show codes


            # verts = [ s.loops[1].curves[1].points[1].coord[1:2] ]
            verts = []
            codes = []
            for (i,c) in enumerate(s.loops[1].curves)
                p1 = points[i]
                p2 = i<npoints ? points[i+1] : points[1]
                c = getcurve(geo, p1, p2)

                if i==1 
                    push!(verts, p1.coord[1:2])
                    push!(codes, MOVETO)
                end

                if c isa Line
                    push!(verts, p2.coord[1:2])
                    push!(codes, LINETO)

                else
                    p3 = p2
                    p2 = c.points[2]
                    P1, P2, P3 = p1.coord, p2.coord, p3.coord
                    r = norm(P3-P2)
                    α = acos(clamp(dot(P1-P2, P3-P2)/r^2, -1, 1))
                    θ = α/2
                    N = normalize(cross(P1-P2, P3-P2))
                    # N = normalize(N)
                    R = Quaternion(cos(θ/2), N[1]*sin(θ/2), N[2]*sin(θ/2), N[3]*sin(θ/2))
                    P = P2 + R*(P1-P2)*conj(R)

                    CP = 2*P - 0.5*P1 - 0.5*P3

                    append!(verts, [CP[1:2], P3[1:2]])
                    append!(codes, [CURVE3, CURVE3])
                end
            end

            path  = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path)

            ec = (0.0, 0.0, 0.0, 0.0)
            ec = (0.4, 0.4, 0.4, 1.0)

            fc = "aliceblue"
            # fc = (0.0, 0.0, 0.0, 0.0)

            push!(edgecolor, ec)
            push!(facecolors, fc)
            push!(all_patches, patch)
        end
        cltn = matplotlib.collections.PatchCollection(all_patches, edgecolor=edgecolor, facecolor=facecolors, lw=lineweight)
        ax.add_collection(cltn) 
    
        # plot all lines
        all_patches = []
        edgecolor   = []
        facecolors  = []
        lineweight  = []
    
        # draw curves
        for c in geo.entities
            c isa Curve || continue

            if c isa Line
                verts = [ p.coord[1:2] for p in c.points ]
                codes = [ MOVETO, LINETO ]
            else
                p1 = c.points[1]
                p2 = c.points[2]
                p3 = c.points[3]

                P1, P2, P3 = p1.coord, p2.coord, p3.coord
                r = norm(P3-P2)
                α = acos(clamp(dot(P1-P2, P3-P2)/r^2, -1, 1))
                θ = α/2
                N = normalize(cross(P1-P2, P3-P2))
                R = Quaternion(cos(θ/2), N[1]*sin(θ/2), N[2]*sin(θ/2), N[3]*sin(θ/2))
                P = P2 + R*(P1-P2)*conj(R)

                CP = 2*P - 0.5*P1 - 0.5*P3

                verts = [P1[1:2], CP[1:2], P3[1:2]]
                codes = [MOVETO, CURVE3, CURVE3]

            end
            path  = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path)
            ec = (0.4, 0.4, 0.4, 1.0)
            # ec = (0.0, 0.0, 0.0, 0.0)
            fc = (0.0, 0.0, 0.0, 0.0)
                
            push!(edgecolor, ec)
            push!(facecolors, fc)
            push!(all_patches, patch)
    
        end
    
        cltn = matplotlib.collections.PatchCollection(all_patches, edgecolor=edgecolor, facecolor=facecolors, lw=0)
        ax.add_collection(cltn) 
    elseif ndim==3
        all_verts = []
        edgecolors = []
        facecolors = []
        lineweight = []

        for s in geo.entities
            s isa Surface || continue
            
            points = getpoints(s.loops[1])
            npoints = length(points)
            verts = [ p.coord for p in points ]
            ec = (0, 0, 0, 1.0)
            fc = "aliceblue"
            
            push!(edgecolors, ec)
            push!(facecolors, fc)
            push!(all_verts, verts)
        end

        cltn = art3D.Poly3DCollection(all_verts, facecolor=facecolors, edgecolor=edgecolors, lw=0.5, alpha=1) # ! lineweight is not working in Poly3DCollection
        ax.add_collection3d(cltn)
    end

    # Draw node numbers
    for o in geo.entities
        if o isa Point
            x = o.coord[1] + 0.01*ll
            y = o.coord[2] - 0.01*ll
            z = o.coord[3] - 0.01*ll
            ndim==2 && ax.text(x, y, o.id, va="top", ha="left", backgroundcolor="none")
            ndim==3 && ax.text(x, y, z, o.id, fontsize=4, va="center", ha="center", backgroundcolor="none")
        elseif o isa Line
            verts = [ p.coord for p in o.points ]
            C = sum(verts)/length(verts)
            x = C[1] + 0.01*ll
            y = C[2] - 0.01*ll
            z = C[3] - 0.01*ll
            ndim==2 && ax.text(x, y, o.id, color="blue", va="top", ha="left", backgroundcolor="none")
            ndim==3 && ax.text(x, y, z, o.id, fontsize=4, color="blue", va="center", ha="center", backgroundcolor="none")
        elseif o isa Curve
            p1 = o.points[1]
            p2 = o.points[2]
            p3 = o.points[3]

            P1, P2, P3 = p1.coord, p2.coord, p3.coord
            r = norm(P3-P2)
            α = acos(clamp(dot(P1-P2, P3-P2)/r^2, -1, 1))
            θ = α/2
            N = normalize(cross(P1-P2, P3-P2))
            # N = normalize(N)
            R = Quaternion(cos(θ/2), N[1]*sin(θ/2), N[2]*sin(θ/2), N[3]*sin(θ/2))
            C = P2 + R*(P1-P2)*conj(R)
            # C = P2 + 1.3*r*normalize(P-P2)
            # @show C

            x = C[1] + 0.01*ll
            y = C[2] - 0.01*ll
            z = C[3] - 0.01*ll
            ndim==2 && ax.text(x, y, o.id, color="blue", va="top", ha="left", backgroundcolor="none")
            ndim==3 && ax.text(x, y, z, o.id, fontsize=4, color="blue", va="center", ha="center", backgroundcolor="none")
        elseif o isa Surface
            points = getpoints(o.loops[1])
            verts = [ p.coord for p in points ]
            C = sum(verts)/length(verts)
            x = C[1]
            y = C[2]
            z = C[3]
            ndim==2 && ax.text(x, y, o.id, color="red", va="center", ha="center", backgroundcolor="none")
            ndim==3 && ax.text(x, y, z, o.id, fontsize=4, color="red", va="center", ha="center", backgroundcolor="none")
        end
    end

    # if nodelabels
    #     nnodes = length(X)
    #     for i in 1:nnodes
    #         x = X[i] + 0.01*L
    #         y = Y[i] - 0.01*L
    #         z = Z[i] - 0.01*L
    #         if ndim==3
    #             ax.text(x, y, z, i, va="center", ha="center", backgroundcolor="none")
    #         else
    #             ax.text(x, y, i, va="top", ha="left", backgroundcolor="none")
    #         end
    #     end
    # end

  #  # Draw cell numbers
  #  if celllabels && ndim==2
  #      for i in 1:ncells
  #          coo = getcoords(cells[i])
  #          x = mean(coo[:,1])
  #          y = mean(coo[:,2])
  #          ax.text(x, y, i, va="top", ha="left", color="blue", backgroundcolor="none", size=8)
  #      end
  #  end


    if filename!=""
        _, format = splitext(filename)
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.00, format=format[2:end])

        if crop && format==".pdf"
            if Sys.islinux()
                cmd = `pdfcrop $filename $filename`
                out = Pipe()
                err = Pipe()
                run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
                close(out.in)
                close(err.in)
            else
                alert("crop option is not available on $(Sys.KERNEL)")
            end
        end

        quiet || info("file $filename saved")

        if copypath!=""
            if isdir(copypath)
                copyfile = joinpath(copypath, basename(filename))
            else
                copyfile = copypath
            end
            cp(filename, copyfile, force=true)
            quiet || info("file $copyfile saved")
        end

        # isinter && plt.ion()
    end

    return gcf()
end