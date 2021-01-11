export mplot, mplotcolorbar

const MOVETO = 1
const LINETO = 2
const CURVE3 = 3
const CURVE4 = 4
const CLOSEPOLY = 79

function plot_data_for_cell2d(points::Array{Array{Float64,1},1}, shape::ShapeType)

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
        for i=1:n
            p2 = i<n ? points[i+1] : points[1]
            push!(verts, p2)
            push!(codes, LINETO)
        end
    elseif shape in (TRI6, QUAD8, QUAD9)
        n = shape==TRI6 ? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
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
        for i=1:n
            p1 = points[i]
            p2 = i<n ? points[i+1] : points[1]
            p3 = points[2*i+3]
            p4 = points[2*i+4]
            cp2 = 1/6*(-5*p1+18*p2-9*p3+2*p4)
            cp3 = 1/6*( 2*p1-9*p2+18*p3-5*p4)
            append!(verts, [cp2, cp3, p2])
            append!(codes, [CURVE4, CURVE4, CURVE4])
        end
    else
        error("plot_data_for_cell2d: Not implemented for ", shape.name)
    end

    return verts, codes
end

#function plot_data_for_cell3d(points::Array{Array{Float64,1},1}, shape::ShapeType)
function plot_data_for_cell3d(points::Array{Vec3,1}, shape::ShapeType, V::Vec3=Vec3(1,0,0), width::Float64=0.0)
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
            J = Vec3(D*X)
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
        
        for k=1:4
            for r in range(0.0, 1.0, length=n)[1:n-1]
                R = corners[k] + r*(corners[k+1]-corners[k])
                N = shape.func(R)
                push!(verts, Vec3(N'*X))
            end
        end
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

Plots an array of blocks using `PyPlot` backend. If filename is provided it saves the output in a pdf file.

# Arguments

`blocks` : An array of `Block` objects. Subarrays are also supported.

`filename` = ""` : If provided, a pdf file with the output is saved

# See also

See documentation of `mplot(mesh, filename="", kwargs...)` for details about keyword arguments.
"""
function mplot(items::Union{Block, Array}, filename::String=""; args...)
    # Get list of blocks and check type
    blocks = unfold(items) # only close if not saving to file

    for item in blocks
        isa(item, Block) || error("mplot: Block object expected")
    end

    # Using Nodes and Cell types
    nodes = Array{Node,1}()
    cells  = Array{Cell,1}()

    for bl in blocks
        append!(nodes, bl.nodes)

        if bl.shape.family==SOLID_SHAPE
            cell = Cell(bl.shape, bl.nodes)
            push!(cells, cell)
        elseif bl.shape.family==LINE_SHAPE
            lines = [ Cell(LIN2, bl.nodes[i-1:i]) for i=2:length(bl.nodes)]
            append!(cells, lines)
        else
            continue
        end

    end

    # Get ndim
    ndim = 1
    for node in nodes
        node.coord.y != 0.0 && (ndim=2)
        node.coord.z != 0.0 && (ndim=3; break)
    end

    mesh = Mesh()
    mesh.ndim = ndim
    mesh.nodes = nodes
    mesh.elems = cells
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
        face.shape.family == SOLID_SHAPE || continue # only surface cells
        face_idx = faces_dict[hash(face)]
        for edge in get_edges(face)
            hs = hash(edge)
            edge0 = get(edge_dict, hs, nothing)
            if edge0===nothing
                edge_dict[hs] = edge
            else
                delete!(edge_dict, hs)
                n1 = normals[face_idx] # normal from face
                face0_idx = faces_dict[hash(edge0.oelem)]
                n2 = normals[face0_idx] # normal from edge0's parent
                α = 180 - acos( abs(clamp(dot(n1,n2),-1,1)) )*180/pi
                α = round(α, digits=2)
                α<=angle && push!(main_edges, edge)
            end
        end
    end

    return main_edges
end


import PyCall: PyObject, pyimport, @pydef # required

"""
    mplot(mesh, filename="", kwargs...)

Plots a `mesh` using `PyPlot` backend. If `filename` is provided it writes a pdf file containing the plot.

# Arguments

`mesh` : A finite element mesh

`filename = ""` : If provided, a pdf file with the output is saved

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

`colormap      = nothing` : Colormap according to PyPlot

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

`outline       = true` : Highlight main edges of 3D meshes in the pdf output

`outlineangle  = 100` : Limit angle to identify main edges

`figsize       = (3,3.0)` : Figure size

`leaveopen     = false` : If true, leaves the plot open so other drawings can be added
"""
function mplot(
               mesh    ::AbstractMesh,
               filename::String    = "";
               axis                = true,
               lw                  = 0.4,
               rodlw               = 1.0,
               rodcolor            = "indianred",
               markers             = false,
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
               opacity             = 1.0,
               lightvector         = nothing,
               colormap            = nothing,
               colormaplims        = (0.0,1.0),
               shrinkcolors        = false,
               darkcolors          = false,
               lightcolors         = false,
               vividcolors         = false,
               divergingcolors     = false,
               colorbarscale       = 0.9,
               colorbarlabel       = "",
               colorbarlocation    = "right",
               colorbarorientation = "vertical",
               colorbarpad         = 0.0,
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
               verbose             = true,
               silent              = false,
               copypath            = ""
              )

    headline("Mesh plotting")
    if filename!=""
        message("generating plot to file $filename")
    end

    if verbose
        hint("Optional arguments:", level=2)
        options = "axis, lw, markers, nodelabels, celllabels, opacity, field,
                   fieldmult, fieldlims, vectorfield, arrowscale, colormap, colorbarscale,
                   colorbarlabel, colorbarlocation, colorbarorientation, colorbarpad, 
                   warpscale, hicells      , elev, azim, dist, outline, outlineangle,
                   figsize, leaveopen, verbose"
        hint(options, level=3)
        hint("Available node fields:", level=2)
        hint(join(keys(mesh.node_data), ", "), level=3)
        hint("Available element fields:", level=2)
        hint(join(keys(mesh.elem_data), ", "), level=3)
    end

    # Get a copy
    mesh = copy(mesh)

    if hicells       isa Int 
        hicells       = [ hicells       ]
    end

    # Get initial info from mesh
    ndim = mesh.ndim
    if ndim==2
        node_data = mesh.node_data
        elem_data = mesh.elem_data

        # filter bulk and line elements
        areacells  = [ elem for elem in mesh.elems if elem.shape.ndim==2 ]
        linecells  = [ cell for cell in mesh.elems if cell.shape.family==LINE_SHAPE]
        newcells   = [ areacells; linecells ]
        c_ids      = [ [c.id for c in areacells]; [c.id for c in linecells] ]
        newnodes   = [ p for c in newcells for p in c.nodes ]
        pt_ids     = [ p.id for p in newnodes ]

        # update data
        for (field, data) in mesh.node_data
            node_data[field] = data[pt_ids,:]
        end
        for (field, data) in mesh.elem_data
            elem_data[field] = data[c_ids]
        end

        # nodes and cells
        nodes = newnodes
        cells = newcells
        
        # connectivities
        id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(nodes) )
        connect = [ [ id_dict[p.id] for p in c.nodes ] for c in cells ] # Do not use type specifier inside comprehension to avoid problem with Revise


        
        # connect = [ [ p.id for p in c.nodes ] for c in cells ] # Do not use type specifier inside comprehension to avoid problem with Revise
        # id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(nodes) )
    else
        node_data = OrderedDict{String,Array}()
        elem_data  = OrderedDict{String,Array}()

        # get surface cells and update
        volume_cells = [ elem for elem in mesh.elems if elem.shape.ndim==3 ]
        areacells    = [ elem for elem in mesh.elems if elem.shape.ndim==2 ]
        scells       = get_surface(volume_cells)
        linecells    = [ cell for cell in mesh.elems if cell.shape.family==LINE_SHAPE]
        outlinecells = outline ? get_outline_edges(scells) : Cell[]

        # linecells = []

        newcells = [ scells; areacells; linecells ]
        oc_ids = [ [c.oelem.id for c in scells]; [c.id for c in linecells]; [c.id for c in areacells] ]



        newnodes = [ p for c in newcells for p in c.nodes ]
        pt_ids = [ p.id for p in newnodes ]

        # update data
        for (field, data) in mesh.node_data
            node_data[field] = data[pt_ids,:]
        end
        for (field, data) in mesh.elem_data
            elem_data[field] = data[oc_ids]
        end

        # nodes and cells
        nodes = newnodes
        cells = newcells

        # connectivities
        id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(nodes) )
        connect = [ [ id_dict[p.id] for p in c.nodes ] for c in cells ]

        # observer and light vectors
        V = Vec3( cosd(elev)*cosd(azim), cosd(elev)*sind(azim), sind(elev) )

        lightvector===nothing && (lightvector=V) 
        if lightvector isa AbstractArray
            L = lightvector
        else
            error("mplot: lightvector must be a vector.")
        end


    end

    ncells = length(cells)
    nnodes = length(nodes)
    #pts = [ [p.coord.x, p.coord.y, p.coord.z] for p in nodes ]
    #XYZ = [ pts[i][j] for i=1:nnodes, j=1:3]


    # All nodes coordinates
    if warpscale>0 
        if haskey(node_data, "U")
            U = node_data["U"]
            coords = [ node.coord for node in nodes ]
            #display(U)
            #@show size(U)
            #@show length(nodes)
            for (i,node) in enumerate(nodes)
                #@show size(node.coord)
                #@show size(U[i,:]')
                #@show node.coord 
                node.coord = coords[i] + warpscale*U[i,:]
                #@show warpscale.*U[i,:]
                #@show node.coord 
                #error()
            end
            #XYZ .+= warpscale.*node_data["U"]
        else
            alert("mplot: Vector field U not found for warp.")
        end
    end

    limX = collect(extrema( node.coord.x for node in nodes ))
    limY = collect(extrema( node.coord.y for node in nodes ))
    limZ = collect(extrema( node.coord.z for node in nodes ))
    #limX = limX + 0.05*[-1, 1]*norm(limX)
    #limY = limY + 0.05*[-1, 1]*norm(limY)
    #limZ = limZ + 0.05*[-1, 1]*norm(limZ)

    #X = XYZ[:,1]
    #Y = XYZ[:,2]
    #Z = XYZ[:,3]

    #limX = collect(extrema(X))
    #limY = collect(extrema(Y))
    #limZ = collect(extrema(Z))
    #limX = limX + 0.05*[-1, 1]*norm(limX)
    #limY = limY + 0.05*[-1, 1]*norm(limY)
    #limZ = limZ + 0.05*[-1, 1]*norm(limZ)
    ll = max(norm(limX), norm(limY), norm(limZ))

    # if ndim==3 && outline # move outline cells towards observer
    #     θ, γ = (azim+0)*pi/180, elev*pi/180
    #     ΔX = [ cos(θ)*cos(γ), sin(θ)*cos(γ), sin(γ) ]*0.015*ll

    #     for edge in outlinecells
    #         for node in edge.nodes
    #             node = Node(edge.nodes[1].coord + ΔX)
    #         end
    #         # edge.nodes[2] = Node(edge.nodes[2].coord + ΔX)
    #     end
    #  end


    # Lazy import of PyPlot
    @eval import PyPlot:plt, matplotlib, figure, art3D, Axes3D, ioff, ColorMap
    @eval ioff()

    # fix PyPlot
    @eval import PyPlot:getproperty, LazyPyModule
    if ! @eval hasmethod(getproperty, (LazyPyModule, AbstractString))
        @eval Base.getproperty(lm::LazyPyModule, s::AbstractString) = getproperty(PyObject(lm), s)
    end

    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=6)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", lw=0.5)
    plt.rc("legend", fontsize=6)
    plt.rc("figure", figsize=figsize) # suggested size (4.5,3)



    # Configure plot
    if ndim==3
        ax = @eval Axes3D(figure())
        try
            ax.set_aspect("equal")
        catch err
            alert("mplot: Could not set aspect ratio to equal")
        end

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

        # Labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        axis == false && plt.axis("off")
    else
        ax = plt.axes()
        ax.set_aspect("equal", "datalim")

        # Set limits
        # limX .+= [ -0.01*ll, 0.1*ll]
        # limY .+= [ -0.01*ll, 0.1*ll]
        ax.set_xlim(limX...)
        ax.set_ylim(limY...)

        # Labels
        ax.set_xlabel.("x")
        ax.set_ylabel.("y")
        axis == false && plt.axis("off")
    end

    has_field = field != nothing
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
            fvals = [ mean(data[connect[i]]) for i=1:ncells ]
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
        for (i,cell) in enumerate(cells)
            cell.shape.family == LINE_SHAPE || continue
            if fvals[i]!=0.0
                has_line_field = true
                break
            end
        end
    end

    # has_line_field = false

    # Plot cells
    if ndim==3
        # Plot cells
        all_verts = []
        edgecolor = []
        facecolor = []
        lineweight = []

        for (i,cell) in enumerate(cells)
            shape = cell.shape
            points = [ node.coord for node in cell.nodes ]

            if shape.family==SOLID_SHAPE
                verts = plot_data_for_cell3d(points, shape)

                if !has_field || has_line_field
                    fc = (0.94, 0.97, 1.0, opacity)
                    ec = (0.4, 0.4, 0.4, 1-0.75*(1-opacity))
                else
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    fc = cmap(v)
                    ec = Tuple( (fc .+ [0.2, 0.2, 0.2, opacity])./2 )
                end
                # push!(lineweight, lw)
                ew = lw

                N = get_facet_normal(cell)
                R = 2*N*dot(L,N) - L
                
                #@show 0
                #@show dot(L,N)
                #@show (1-dot(V,R))/2
                
                #f = 0.7+0.3*abs(dot(L,N)*(1-dot(V,R))/2)
                #f = 0.6+0.3*abs(dot(L,N)) + 0.2*(1-dot(V,R))/2
                
                # f = 0.6+0.2*abs(dot(L,N)) + 0.2*(1+dot(V,R))/2
                f = 0.8 + 0.1*abs(dot(L,N)) + 0.1*(1+dot(V,R))/2
                f = min(f, 1.0)

                # f = 0.9+0.05*abs(dot(L,N)) + 0.05*(1+dot(V,R))/2
                # f = 1.0
                fc = (f*fc[1], f*fc[2], f*fc[3], opacity)
            elseif shape.family==LINE_SHAPE
                verts = plot_data_for_cell3d(points, shape, V, 0.0075*ll*rodlw)

                if has_line_field
                    v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                    ec = cmap(v)
                else
                    ec = rodcolor 
                end
                fc = (0.0, 0.0, 0.0, 0.0)
                fc = ec
                ew = rodlw  #! lineweight is not working in Poly3DCollection
                ew = lw
                # ew = 0.0
            end
            
            push!(all_verts, verts)
            push!(edgecolor, ec)
            push!(facecolor, fc)
            push!(lineweight, ew)
        end

        for cell in outlinecells
            ΔX = 0.01*ll*V
            points = [ node.coord+ΔX for node in cell.nodes ]
            # verts = plot_data_for_cell3d(points, cell.shape)
            verts = plot_data_for_cell3d(points, cell.shape, V, 0.001)

            # ec = "black"
            ec = "dimgray"
            fc = (0.0, 0.0, 0.0, 0.0)

            push!(all_verts, verts)
            push!(lineweight, lw*1.0)
            push!(edgecolor, ec)
            push!(facecolor, fc)
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
                push!(facecolor, fc)
            end
        end

        cltn = @eval art3D[:Poly3DCollection]($all_verts, facecolor=$facecolor, edgecolor=$edgecolor, lw=$lineweight, alpha=$opacity) # ! lineweight is not working in Poly3DCollection
        ax.add_collection3d(cltn)

        if has_field && colorbarscale>0
            normalized = matplotlib.colors.Normalize(fieldlims[1], fieldlims[2])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalized)
            sm.set_array([])

            cbar = plt.colorbar(sm, label=colorbarlabel, shrink=colorbarscale, aspect=10*colorbarscale*figsize[2], 
                                pad=colorbarpad, location=colorbarlocation)

            cbar.ax.tick_params(labelsize=6)
            cbar.outline.set_linewidth(0.0)
            cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
            cbar.update_ticks()
            cbar.solids.set_alpha(1)
        end


    elseif ndim==2

        all_patches = []
        edgecolor = []
        facecolor = []
        lineweight = []

        # for i=1:ncells
        for (i,cell) in enumerate(cells)
            shape = cell.shape
            points = [ node.coord[1:2] for node in cell.nodes ]
            verts, codes = plot_data_for_cell2d(points, shape)
            path  = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path)


            if shape.family==SOLID_SHAPE
                if !has_field || has_line_field
                    fc = (0.94, 0.97, 1.0, 1.0)
                    ec = (0.4, 0.4, 0.4, 1.0)

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
            elseif shape.family==LINE_SHAPE
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
            push!(facecolor, fc)
            push!(all_patches, patch)
        end


        cltn = matplotlib.collections.PatchCollection(all_patches, edgecolor=edgecolor, facecolor=facecolor, lw=lineweight)
        ax.add_collection(cltn)
        
        if has_field && colorbarscale>0
            normalized = matplotlib.colors.Normalize(fieldlims[1], fieldlims[2])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalized)
            sm.set_array([])

            h = colorbarorientation=="vertical" ? figsize[2] : figsize[1]
            h = norm(figsize)
            cbar = plt.colorbar(sm, label=colorbarlabel, shrink=colorbarscale, aspect=4*colorbarscale*h, 
                                # format="%.2f", 
                                pad=colorbarpad, orientation=colorbarorientation)
            cbar.ax.tick_params(labelsize=6)
            cbar.outline.set_linewidth(0.0)
            cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
            cbar.update_ticks()
            cbar.solids.set_alpha(1)
        end
    end

    # Draw nodes
    if markers
        if ndim==3
            X = [ node.coord.x for node in nodes ]
            Y = [ node.coord.y for node in nodes ]
            Z = [ node.coord.z for node in nodes ]
            ax.scatter(X, Y, Z, color="black", marker="o", s=ms)
        else
            X = [ node.coord.x for cell in areacells for node in cell.nodes ]
            Y = [ node.coord.y for cell in areacells for node in cell.nodes ]
            plt.plot(X, Y, color="black", marker="o", mfc="black", fillstyle="full", markersize=ms, mew=0, lw=0, clip_on=false)
        end
    end

   # Node markers on line cells
    if rodmarkers && ndim==2
        X = [ node.coord.x for line in linecells for node in line.nodes ]
        Y = [ node.coord.y for line in linecells for node in line.nodes ]
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

  #  # Draw node numbers
  #  if nodelabels
  #      nnodes = length(X)
  #      for i=1:nnodes
  #          x = X[i] + 0.01*L
  #          y = Y[i] - 0.01*L
  #          z = Z[i] - 0.01*L
  #          if ndim==3
  #              ax.text(x, y, z, i, va="center", ha="center", backgroundcolor="none")
  #          else
  #              ax.text(x, y, i, va="top", ha="left", backgroundcolor="none")
  #          end
  #      end
  #  end

  #  # Draw cell numbers
  #  if celllabels && ndim==2
  #      for i=1:ncells
  #          coo = getcoords(cells[i])
  #          x = mean(coo[:,1])
  #          y = mean(coo[:,2])
  #          ax.text(x, y, i, va="top", ha="left", color="blue", backgroundcolor="none", size=8)
  #      end
  #  end

    if ndim==3
        ax.view_init(elev=elev, azim=azim)
        ax.dist = dist
    end

    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.00, format="pdf")

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

        verbose && info("file $filename saved")

        if copypath!=""
            if isdir(copypath)
                copyfile = joinpath(copypath, basename(filename))
            else
                copyfile = copypath
            end
            cp(filename, copyfile, force=true)
            verbose && info("file $copyfile saved")
        end
    end

    # Do not close if in IJulia
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return
    end

    leaveopen || plt.close("all")

    return

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
    #colormaplims        = (0.0,1.0),
    #shrinkcolor         = false,
    #darkcolor           = false,
    #lightcolor          = false,
    #vividcolor          = false,
    #divergingcolor      = false,
)    

    # Lazy import of PyPlot
    @eval import PyPlot:plt, matplotlib, figure, art3D, Axes3D, ioff, ColorMap, gca, gcf
    @eval ioff()

    # fix PyPlot
    @eval import PyPlot:getproperty, LazyPyModule
    if ! @eval hasmethod(getproperty, (LazyPyModule, AbstractString))
        @eval Base.getproperty(lm::LazyPyModule, s::AbstractString) = getproperty(PyObject(lm), s)
    end

    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", lw=0.5)
    plt.rc("legend", fontsize=fontsize)
    plt.rc("figure", figsize=figsize) # suggested size (4.5,3)

    fig   = @eval plt.figure()
    scale = 0.5
    axes  = @eval gcf().add_axes([0, 0, $scale/$aspect, $scale])

    cmap = matplotlib.cm.get_cmap(colormap)
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
    mplot_linear(elems, filename="", kwargs...)

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
    axis             = true,
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
    
    @eval import PyPlot:plt, matplotlib, figure, gca, ioff
    @eval ioff()

    @assert barscale>0

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("figure", figsize=figsize)
    maxv = 0.0

    lines = [ elem.shape.family==LINE_SHAPE ? elem : elem.linked_elems[2] for elem in elems ]
    coords = get_coords(get_nodes(lines))
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
    
    #@show xidx
    #@show yidx


    # Get maximum values
    vmax = 0.0
    vmin = 0.0
    xmin = 0.0
    xmax = 0.0
    for elem in elems
        for ip in elem.ips
            v = ip_state_vals(elem.mat, ip.state)[Symbol(field)]
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
    dn = 0.01*xwidth

    for elem in elems
        line = elem.shape.family==LINE_SHAPE ? elem : elem.linked_elems[2]
        X = [ node.coord[xidx] for node in line.nodes ]
        Y = [ node.coord[yidx] for node in line.nodes ]
        plt.plot(X, Y, "tab:gray", lw=3, solid_capstyle="round", zorder=2) # plot line

        field===nothing && continue

        V = [ X[2]-X[1], Y[2]-Y[1] ]
        normalize!(V)
        N = [ -V[2], V[1] ]

        ips = elem.ips
        for ip in ips
            x = ip.coord[xidx]
            y = ip.coord[yidx]
            v = ip_state_vals(elem.mat, ip.state)[Symbol(field)]*sign(fieldmult)
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

    if showscale
        xmin, xmax = gca().get_xlim()
        xwidth = xmax-xmin # recalculate xwidth according to plot

        refval = round_for_scale(maxv)
        scalelen = refval*fieldtrans*barscale*1.0/xwidth

        xpos, ypos = scalepos

        scaletitle = "Scale"
        if fullscale
            scaletitle *= fieldunits!="" ? " [$fieldunits]" : "" 
        end
        plt.text(xpos, ypos, scaletitle, transform = gca().transAxes, ha="left", va="top")
        
        ypos -= 0.06
        ypos0 = ypos
        if vmin<0
            plt.plot([xpos, xpos+scalelen], [ypos, ypos], "tab:blue", ls="-", lw=1, solid_capstyle="butt", zorder=1, transform = gca().transAxes)
            ypos -= 0.02
        end
        if vmax>0
            plt.plot([xpos, xpos+scalelen], [ypos, ypos], "tab:red", ls="-", lw=1, solid_capstyle="butt", zorder=1, transform = gca().transAxes)
            ypos -= 0.02
        end
        

        if fullscale
            ticks = (0.0, 0.5, 1.0)
        else
            ticks = (0.0, 1.0)
        end
        for f in ticks
            plt.plot([xpos+f*scalelen, xpos+f*scalelen], [ypos+0.01, ypos], "black", ls="-", lw=0.5, zorder=1, transform = gca().transAxes)
        end

        ypos -= 0.01
        if fullscale
            for f in ticks
                plt.text(xpos+f*scalelen, ypos, "$(f*refval*abs(fieldmult))", transform=gca().transAxes, ha="center", va="top")
            end
        else
            plt.text(xpos+0.5*scalelen, ypos, "$(refval*abs(fieldmult)) $fieldunits", transform=gca().transAxes, ha="center", va="top")
        end
    end

    # gca().set_aspect("equal", "datalim")
    gca().set_aspect("equal", "box")

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

        if minmax
            vrmin = round(vmin*fieldmult, sigdigits=4)
            vrmax = round(vmax*fieldmult, sigdigits=4)
            if fieldmult<0
                vrmin, vrmax = vrmax, vrmin
            end

            tmin  = "\$$(vrmin>0 ? "+" : "" )$vrmin\$ $fieldunits"
            tmax  = "\$$(vrmax>0 ? "+" : "" )$vrmax\$ $fieldunits"
            plt.text(0.5, 0.98, "min: $tmin   max: $tmax", transform = gca().transAxes, ha="center", va="top")
        end
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

end

