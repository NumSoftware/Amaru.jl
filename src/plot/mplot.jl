export mplot

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
        error("plot_data_for_cell2d: Not implemented for ", shape)
    end

    return verts, codes
end

function plot_data_for_cell3d(points::Array{Array{Float64,1},1}, shape::ShapeType)
    if shape == LIN2
        verts = points
    elseif shape == LIN3
        verts = points[[1,3,2]]
    elseif shape in (TRI3, QUAD4)
        verts = points
    elseif shape == TRI6
        verts = points[[1,4,2,5,3,6]]
    elseif shape in (QUAD8, QUAD9)
        verts = points[[1,5,2,6,3,7,4,8]]
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
            if edge0==nothing
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


function get_facet_normal(face::AbstractCell)
    ndim = 1 + face.shape.ndim
    C = get_coords(face, ndim)

    if ndim==2
        C .+= [pi pi^1.1]
    else
        C .+= [pi pi^1.1 pi^1.2]
    end

    # calculate the normal
    I = ones(size(C,1))
    N = pinv(C)*I # best fit normal
    normalize!(N) # get unitary vector

    return N
end


function get_surface_based_on_displacements(mesh::Mesh)

    surf_dict = Dict{UInt64, Cell}()
    W = mesh.node_data["wn"]
    #U = mesh.node_data["U"]
    maxW = maximum(W)

    disp(face) = begin
        map = [p.id for p in face.nodes]
        #d=mean(W[map])
        d=max(maximum(W[map]), 0)

        return d
    end
    #disp(face) = mean( W[p.id] for p in face.nodes )
    #disp(face) = begin
        #map = [p.id for p in face.nodes]
        #mean(U[map,:])
    #end

    # Get only unique faces. If dup, original and dup are deleted
    for cell in mesh.elems
        for face in get_faces(cell)
            hs = hash(face)
            dface = disp(face)
            #dface < 0.0000000*maxW && continue
            if haskey(surf_dict, hs)
                d = norm(dface-disp(surf_dict[hs]))
                #d = disp(face)
                #d = abs(disp(face)-disp(surf_dict[hs]))
                #if d<0.1*maxW
                #if d<0.0000
                if false
                    delete!(surf_dict, hs)
                else
                    surf_dict[hs] = face
                end
            else
                surf_dict[hs] = face
            end
        end
    end

    return [ face for face in values(surf_dict) ]

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

`nodemarkers  = false` : If true, shows node markers

`nodelabels   = false` : If true, shows node labels

`celllabels    = false` : If true, shows cell labels

`opacity       = 1.0`   : Opacity,

`field         = nothing` : If provided, plots corresponding field

`fieldscale    = 1.0` : Factor multiplied to `field` values

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

`highlightcell = 0` : Cell number to be highlighted

`elev          = 30.0` : 3D plot elevation

`azim          = 45.0` : 3D plot azimute

`dist          = 10.0` : 3D plot distance from observer

`outline       = true` : Highlight main edges of 3D meshes in the pdf output

`outlineangle  = 100` : Limit angle to identify main edges

`figsize       = (3,2.5)` : Figure size

`leaveopen     = false` : If true, leaves the plot open so other drawings can be added
"""
function mplot(
               mesh::AbstractMesh,
               filename::String = "";
               axis             = true,
               lw               = 0.3,
               lineslw          = 1.0,
               nodemarkers     = false,
               ms              = 0.5,
               nodelabels      = false,
               celllabels       = false,
               field            = nothing,
               fieldscale       = 1.0,
               fieldlims        = (),
               vectorfield      = nothing,
               arrowscale       = 0.0,
               opacity          = 1.0,
               colormap         = nothing,
               colormaplims     = (0.0,1.0),
               shrinkcolors      = false,
               darkcolors        = false,
               lightcolors       = false,
               vividcolors       = false,
               divergingcolors   = false,
               colorbarscale    = 0.9,
               colorbarlabel    = "",
               colorbarlocation = "right",
               colorbarorientation = "vertical",
               colorbarpad      = 0.0,
               warpscale        = 0.0,
               highlightcell    = 0,
               elev             = 30.0,
               azim             = 45.0,
               dist             = 10.0,
               outline          = true,
               outlineangle     = 100,
               figsize          = (3,2.5),
               leaveopen        = false,
               crop             = false,
               verbose          = true,
              )


    # Get initial info from mesh
    ndim = mesh.ndim
    if ndim==2
        node_data = mesh.node_data
        elem_data  = mesh.elem_data
        nodes  = mesh.nodes
        cells   = mesh.elems
        connect = []
        connect = [ [ p.id for p in c.nodes ] for c in cells ] # Do not use type specifier inside comprehension to avoid problem with Revise
        id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(nodes) )
    else
        node_data = OrderedDict{String,Array}()
        elem_data  = OrderedDict{String,Array}()

        # get surface cells and update
        volume_cells = [ elem for elem in mesh.elems if elem.shape.ndim==3 ]
        area_cells = [ elem for elem in mesh.elems if elem.shape.ndim==2 ]
        scells = get_surface(volume_cells)
        #if haskey(mesh.node_data, "U") && haskey(mesh.node_data, "wn")
            # special case when using cohesive elements
            #scells = get_surface_based_on_displacements(mesh)
        #else
            #scells = get_surface(mesh.elems)
        #end
        linecells = [ cell for cell in mesh.elems if cell.shape.family==LINE_SHAPE]

        oc_ids = [ [c.oelem.id for c in scells]; [c.id for c in linecells]; [c.id for c in area_cells] ]

        newcells = [ scells; area_cells; linecells ]
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
        cells  = newcells

        # connectivities
        id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(nodes) )
        connect = [ [ id_dict[p.id] for p in c.nodes ] for c in cells ]
    end

    if verbose
        if filename==""
            printstyled("mplot: generating plot\n", color=:cyan )
        else
            printstyled("mplot: generating plot to file $filename\n", color=:cyan )
        end
        wrap(str::String) = (str=replace(str, r"(\s|\n)+" => " "); replace(str, r".{1,60}( |$)" => s"    \0\n");)
        options = "axis, lw, nodemarkers, nodelabels, celllabels, opacity, field,
                   fieldscale, fieldlims, vectorfield, arrowscale, colormap, colorbarscale,
                   colorbarlabel, colorbarlocation, colorbarorientation, colorbarpad, 
                   warpscale, highlightcell, elev, azim, dist, outline, outlineangle,
                   figsize, leaveopen, verbose"
        printstyled("  Options:\n", wrap(options), color=:light_black)
        printstyled("  Available node fields:\n", wrap(join(keys(node_data), ", ")), color=:light_black)
        printstyled("  Available element fields:\n", wrap(join(keys(elem_data), ", ")), color=:light_black)

    end

    ncells  = length(cells)
    nnodes = length(nodes)
    pts = [ [p.coord.x, p.coord.y, p.coord.z] for p in nodes ]
    XYZ = [ pts[i][j] for i=1:nnodes, j=1:3]

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


    # All nodes coordinates
    if warpscale>0 
        if haskey(node_data, "U")
            XYZ .+= warpscale.*node_data["U"]
        else
            @warn "mplot: Vector field U not found for warp"
        end
    end
    X = XYZ[:,1]
    Y = XYZ[:,2]
    Z = XYZ[:,3]

    limX = collect(extrema(X))
    limY = collect(extrema(Y))
    limZ = collect(extrema(Z))
    limX = limX + 0.05*[-1, 1]*norm(limX)
    limY = limY + 0.05*[-1, 1]*norm(limY)
    limZ = limZ + 0.05*[-1, 1]*norm(limZ)
    L = max(norm(limX), norm(limY), norm(limZ))

    # Configure plot
    if ndim==3
        ax = @eval Axes3D(figure())
        try
            ax.set_aspect("equal")
            #ax.set_proj_type("ortho") # do not ortho
            #ax.set_proj_type("3d") # not necessary
            #axis("scaled") ?
        catch err
            @warn "mplot: Could not set aspect ratio to equal"

            #dump(err)
        end

        # Set limits
        meanX = mean(limX)
        meanY = mean(limY)
        meanZ = mean(limZ)
        limX = [meanX-L/2, meanX+L/2]
        limY = [meanY-L/2, meanY+L/2]
        limZ = [meanZ-L/2, meanZ+L/2]
        ax.set_xlim( meanX-L/2, meanX+L/2)
        ax.set_ylim( meanY-L/2, meanY+L/2)
        ax.set_zlim( meanZ-L/2, meanZ+L/2)
        #ax.scatter](limX, limY, limZ, color="w", marker="o", alpha=0.0)

        # Labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        if axis == false
            ax.set_axis_off()
        end
    else
        #ax = plt.gca()
        #plt.axes().set_aspect("equal", "datalim")
        ax = plt.axes()
        ax.set_aspect("equal", "datalim")

        # Set limits
        ax.set_xlim(limX[1], limX[2])
        ax.set_ylim(limY[1], limY[2])

        # Labels
        ax.set_xlabel.("x")
        ax.set_ylabel.("y")
        if axis == false
            plt.axis("off")
        end
    end

    if colormap==nothing # colormap may be "bone", "plasma", "inferno", etc.
        #cm = colors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],256)
        #colors =  matplotlib.colors]
        #cm = matplotlib.colors].ListedColormap]([(1,0,0),(0,1,0),(0,0,1)],256)

        cdict = Dict("red"   => [(0.0,  0.8, 0.8), (0.5, 0.7, 0.7), (1.0, 0.0, 0.0)],
                     "green" => [(0.0,  0.2, 0.2), (0.5, 0.7, 0.7), (1.0, 0.2, 0.2)],
                     "blue"  => [(0.0,  0.0, 0.0), (0.5, 0.7, 0.7), (1.0, 0.6, 0.6)])

        cmap = matplotlib.colors.LinearSegmentedColormap("my_colormap",cdict,256)
    else
        if colormap isa ColorMap
            cmap = colormap
        else
            colormaps = matplotlib.pyplot.colormaps()
            if colormap in colormaps
                cmap = matplotlib.cm.get_cmap(colormap)
            else
                @error "mplot: Invalid colormap $colormap"
                error("  colormap should be one of:\n", colormaps)
            end
        end
    end

    has_field = field != nothing
    if has_field

        if colormaplims!=(0.0,1.0)
            allcolors = [ cmap(p) for p in range(colormaplims[1], colormaplims[2], length=21) ]
            Q = range(0,1,length=21)

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap",cdict,256)
        end

        if shrinkcolors
            allcolors = [ cmap(p) for p in range(0,1,length=21) ]
            Q = [ 0.5+0.5*(abs((p-0.5)/0.5))^2.5*sign(p-0.5) for p in range(0,1,length=21) ]

            cdict = Dict("red"   => [ (q, c[1], c[1]) for (q,c) in zip(Q,allcolors) ],
                         "green" => [ (q, c[2], c[2]) for (q,c) in zip(Q,allcolors) ],
                         "blue"  => [ (q, c[3], c[3]) for (q,c) in zip(Q,allcolors) ])

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap",cdict,256)

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

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap",cdict,256)
        end


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
        fvals *= fieldscale
        fieldlims==() && (fieldlims = extrema(fvals))

        if !(fieldlims[1]<0 && fieldlims[2]>0)
            divergingcolors = false
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

            cmap = matplotlib.colors.LinearSegmentedColormap("modified_colormap",cdict,256)
        end
    end

    # Plot cells
    if ndim==3
        # Plot line cells
        all_verts  = []
        lfvals = Float64[]
        for i=1:ncells
            shape = cells[i].shape
            shape.family == LINE_SHAPE || continue # only line cells
            con = connect[i]
            nodes = [ XYZ[i,1:3] for i in con ]
            verts = plot_data_for_cell3d(nodes, shape)
            push!(all_verts, verts)
            if has_field
                push!(lfvals, fvals[i])
            end
        end
        has_line_field = sum(abs, lfvals)!=0.0

        line_cltn = @eval art3D[:Line3DCollection]($all_verts, cmap=$cmap, lw=$lineslw, edgecolor="red")
        if has_line_field
            line_cltn.set_array(lfvals)
            line_cltn.set_clim(fieldlims)
        end
        ax.add_collection3d(line_cltn)

        # Plot main edges and surface cells
        all_verts  = []

        # Plot surface cells
        sfvals = Float64[]
        nsurf = 0
        for i=1:ncells
            shape = cells[i].shape
            shape.family == SOLID_SHAPE || continue # only surface cells
            con = connect[i]
            points = [ XYZ[i,1:3] for i in con ]
            verts = plot_data_for_cell3d(points, shape)
            push!(all_verts, verts)
            if has_field
                push!(sfvals, fvals[i])
            end
            nsurf+=1
        end

        has_surf_field = sum(abs, sfvals)!=0.0
        has_field && has_line_field==false && has_surf_field==false && (has_surf_field=true)

        edgecolor = []
        for i=1:ncells
            shape = cells[i].shape
            shape.family == SOLID_SHAPE || continue # only surface cells
            if has_surf_field
                v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                col = Tuple( (cmap(v) .+ [0.2, 0.2, 0.2, opacity])./2 )
            else
                col = (0.4, 0.4, 0.4, 1-0.75*(1-opacity))
            end
            push!(edgecolor, col)
        end

        # Plot main edges
        filename == "" && (outline=false)
        if outline
            edges = get_main_edges(cells, outlineangle)
            θ, γ = (azim+0)*pi/180, elev*pi/180
            ΔX = [ cos(θ)*cos(γ), sin(θ)*cos(γ), sin(γ) ]*0.015*L

            for edge in edges
                id1 = edge.nodes[1].id
                id2 = edge.nodes[2].id
                verts = [ XYZ[id_dict[id1],:], XYZ[id_dict[id2],:] ]
                for v in verts
                    v .+= ΔX
                end
                push!(all_verts, verts)
            end
            #edgecolor = [ fill((0.4, 0.4, 0.4, 1.0), ncells) ; fill((0.0, 0.0, 0.0, 1.0), length(edges)) ]
            edgecolor = [ edgecolor ; fill((0.0, 0.0, 0.0, 1.0), length(edges)) ]
        end

        facecolor = (0.94, 0.97, 1.0, opacity)
        surf_cltn = @eval art3D[:Poly3DCollection]($all_verts, cmap=$cmap, facecolor=$facecolor, edgecolor=$edgecolor, lw=$lw, alpha=$opacity)

        if has_field
            if outline
                sfvals = [ sfvals; fill(mean(sfvals), length(edges)) ]
            end

            if has_surf_field
            #if sum(abs, sfvals) != 0.0
                surf_cltn.set_array(sfvals)
                surf_cltn.set_clim(fieldlims)
                #has_line_field = false
            else
                #has_line_field = true
            end
        end
        ax.add_collection3d(surf_cltn)

        # plot colorbar
        if has_field && colorbarscale>0
            cltn = has_line_field ? line_cltn : surf_cltn
            cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=10*colorbarscale*figsize[2], 
                                #format="%.2f", 
                                pad=colorbarpad, location=colorbarlocation)
            cbar.ax.tick_params(labelsize=6)
            cbar.outline.set_linewidth(0.0)
            cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
            cbar.update_ticks()
            cbar.solids.set_alpha(1)
        end


    elseif ndim==2

        # Plot line cells
        for i=1:ncells
            cells[i].shape.family == LINE_SHAPE || continue # only line cells
            con = connect[i]
            Xs = XYZ[con, 1]
            Ys = XYZ[con, 2]
            color = "red"

            if has_field
                v = fvals[i]
                idx = (v-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                color = cmap(idx)
            end
            plt.plot(Xs, Ys, color=color, lw=lineslw)
        end


        all_patches = []
        for i=1:ncells
            shape = cells[i].shape
            shape.family == SOLID_SHAPE || continue # only surface cells

            con = connect[i]
            points = [ XYZ[i,1:2] for i in con ]
            verts, codes = plot_data_for_cell2d(points, shape)
            path  = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path)
            push!(all_patches, patch)

            if highlightcell==i
                patch = matplotlib.patches.PathPatch(path, facecolor="cadetblue", edgecolor="black", lw=0.5)
                ax.add_patch(patch)
            end
        end

        edgecolor = (0.4, 0.4 ,0.4, 1.0)

        # find edgecolors
        if has_field
            edgecolor = []
            for i=1:ncells
                v = (fvals[i]-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                col = Tuple( (cmap(v) .+ [0.3, 0.3, 0.3, 1.0])./2 )
                push!(edgecolor, col)
            end
        end

        #edgecolor = (0.3, 0.3 ,0.3, 0.6)
        cltn = matplotlib.collections.PatchCollection(all_patches, cmap=cmap, edgecolor=edgecolor, facecolor="aliceblue", lw=lw)
        if has_field && colorbarscale>0
            cltn.set_array(fvals)
            cltn.set_clim(fieldlims)
            h = colorbarorientation=="vertical" ? figsize[2] : figsize[1]
            h = norm(figsize)
            #cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=0.9*20*colorbarscale, format="%.1f", pad=colorbarpad, orientation=colorbarorientation)
            #cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=15, format="%.1f", pad=colorbarpad, orientation=colorbarorientation)
            cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=4*colorbarscale*h, 
                                format="%.2f", pad=colorbarpad, orientation=colorbarorientation)
            cbar.ax.tick_params(labelsize=6)
            cbar.outline.set_linewidth(0.0)
            cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
            cbar.update_ticks()
        end
        ax.add_collection(cltn)
    end

    # Draw nodes
    if nodemarkers
        if ndim==3
            ax.scatter(X, Y, Z, color="black", marker="o", s=ms)
        else
            plt.plot(X, Y, color="black", marker="o", markersize=ms, lw=0)
        end
    end

    # Node markers on line cells
    ids = unique!([ id for i=1:ncells for id in connect[i] if cells[i].shape.family==LINE_SHAPE ])
    if length(ids)>0
        if ndim==3
            ax.scatter(X[ids], Y[ids], Z[ids], marker="o", s=ms, facecolor="black", edgecolor="none", depthshade=false)
        else
            plt.plot(X[ids], Y[ids], color="black", marker="o", markersize=ms, lw=0)
        end
    end

    # Draw arrows
    if vectorfield!=nothing && ndim==2
        data = node_data[vectorfield]
        color = "blue"
        if arrowscale==0
            plt.quiver(X, Y, data[:,1], data[:,2], color=color)
        else
            plt.quiver(X, Y, data[:,1], data[:,2], color=color, scale=1.0/arrowscale)
        end
    end

    # Draw node numbers
    if nodelabels
        nnodes = length(X)
        for i=1:nnodes
            x = X[i] + 0.01*L
            y = Y[i] - 0.01*L
            z = Z[i] - 0.01*L
            if ndim==3
                ax.text(x, y, z, i, va="center", ha="center", backgroundcolor="none")
            else
                ax.text(x, y, i, va="top", ha="left", backgroundcolor="none")
            end
        end
    end

    # Draw cell numbers
    if celllabels && ndim==2
        for i=1:ncells
            coo = getcoords(cells[i])
            x = mean(coo[:,1])
            y = mean(coo[:,2])
            ax.text(x, y, i, va="top", ha="left", color="blue", backgroundcolor="none", size=8)
        end
    end

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
                @warn "crop option is not available on $(Sys.KERNEL)"
            end
        end

        verbose && printstyled("  file $filename saved\n", color=:cyan)
    end

    # Do not close if in IJulia
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return
    end

    leaveopen || plt.close("all")

    return

end


function mplot_colorbar(
                        filename::String = "";
                        colormap         = nothing,
                        fieldlims        = (0.0,1.0),
                        #colormaplims     = (0.0,1.0),
                        #shrinkcolor      = false,
                        #darkcolor        = false,
                        #lightcolor       = false,
                        #vividcolor       = false,
                        #divergingcolor   = false,
                        colorbarscale    = 0.9,
                        colorbarlabel    = "",
                        colorbarlocation = "right",
                        colorbarorientation = "vertical",
                        figsize          = (3,2.5),
                        crop             = false,
                       )

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


    fig = plt.figure()
    colorbarscale = 0.5
    aspect = 25
    axes = fig.add_axes([0, 0, colorbarscale/aspect, colorbarscale])

    cmap = matplotlib.cm.get_cmap(colormap)

    cbar = matplotlib.colorbar.ColorbarBase(axes, 
                                            label=colorbarlabel,
                                            norm=matplotlib.colors.Normalize(fieldlims[1], fieldlims[2]),
                                            #values=fieldlims,
                                            drawedges=false,
                                            #boundaries=fieldlims,
                                            orientation=colorbarorientation, 
                                            cmap=cmap)

    cbar.ax.tick_params(labelsize=7)
    cbar.outline.set_linewidth(0.0)
    cbar.locator = matplotlib.ticker.MaxNLocator(nbins=8)
    cbar.update_ticks()
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
                @warn "crop option is not available on $(Sys.KERNEL)"
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
            n=2.5
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

function mplot_linear(
                      elems,
                      filename   = "";
                      field      = nothing,
                      fieldscale = 1.0,
                      fieldunits = "",
                      barscale   = 1.0,
                      axis       = true,
                      xlabel     = "",
                      ylabel     = "",
                      xlim       = nothing,
                      ylim       = nothing,
                      legendlabels = [],
                      showscale = true,
                      scalepos  = (0.6, 0.05)
                     )
    @eval import PyPlot:plt, matplotlib, figure, gca, ioff
    @eval ioff()

    @assert barscale>0

    plt.rc("font", family="STIXGeneral", size=6)
    plt.rc("mathtext", fontset="cm")
    plt.rc("figure", figsize=(2.5, 2.5))
    maxv = 0.0

    lines = [ elem.shape.family==LINE_SHAPE ? elem : elem.linked_elems[2] for elem in elems ]
    coords = get_coords(get_nodes(lines))
    #sumx = maximum(abs, coords[:,1])
    #sumy = maximum(abs, coords[:,2])
    #sumz = maximum(abs, coords[:,3])

    n = size(coords, 1)
    sumx = sum(coords[:,1])
    sumy = sum(coords[:,2])
    #sumz = sum(coords[:,3])
    avgx = sumx/n
    avgy = sumy/n
    #avgz = sumz/n
    devx = sum((coords[:,1].-avgx).^2)/n
    devy = sum((coords[:,2].-avgy).^2)/n
    #devz = sum((coords[:,3].-avgz)^2)
    #@show devx
    #@show devy
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
    else
        xidx = 1
        yidx = 2
        xlabel=="" && (xlabel=raw"$x$")
        ylabel=="" && (ylabel=raw"$y$")
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

    if field == nothing 
        printstyled("mplot_linear:\n", color=:cyan )
    else
        printstyled("mplot_linear: plotting field $field\n", color=:cyan )
        printstyled("  (min,max): ($(vmin*abs(fieldscale)), $(vmax*abs(fieldscale)))\n", color=:light_black)
    end

    maxv = max(abs(vmin), abs(vmax))
    xwidth = 1.3*(xmax-xmin) # estimative of xwidth
    fieldmult = 0.2*xwidth/maxv

    for elem in elems
        line = elem.shape.family==LINE_SHAPE ? elem : elem.linked_elems[2]
        X = [ node.coord[xidx] for node in line.nodes ]
        Y = [ node.coord[yidx] for node in line.nodes ]
        plt.plot(X, Y, "tab:gray", lw=3, solid_capstyle="round", zorder=2) # plot line

        field==nothing && continue

        V = [ X[2]-X[1], Y[2]-Y[1] ]
        normalize!(V)
        N = [ -V[2], V[1] ]

        ips = elem.ips
        for ip in ips
            x = ip.coord[xidx]
            y = ip.coord[yidx]
            v = ip_state_vals(elem.mat, ip.state)[Symbol(field)]*sign(fieldscale)
            X = [ x, x+N[1]*v*fieldmult*barscale ]
            Y = [ y, y+N[2]*v*fieldmult*barscale ]
            color = v>0 ? "tab:red" : "tab:blue"
            plt.plot(X,Y, color, ls="-", lw=1, solid_capstyle="round", zorder=1)
            plt.plot(x,y,"k+", ms=3, mew=0.5, zorder=3)
        end
    end

    if showscale
        xmin, xmax = gca().get_xlim()
        xwidth = xmax-xmin # recalculate xwidth according to plot

        refval = round_for_scale(maxv)
        scalelen = refval*fieldmult*barscale*1.0/xwidth

        xpos, ypos = scalepos

        plt.text(xpos, ypos, "$(refval*abs(fieldscale)) $fieldunits", transform = gca().transAxes)
        dy = 0.05
        if vmin<0
            plt.plot([xpos, xpos+scalelen], [ypos+dy, ypos+dy], "tab:blue", ls="-", lw=1, solid_capstyle="round", zorder=1, transform = gca().transAxes)
            dy += 0.02
        end
        if vmax>0
            plt.plot([xpos, xpos+scalelen], [ypos+dy, ypos+dy], "tab:red", ls="-", lw=1, solid_capstyle="round", zorder=1, transform = gca().transAxes)
            dy += 0.02
        end
        plt.text(xpos, ypos+dy, "Scale", transform = gca().transAxes)
    end

    #gca().set_aspect("equal", "datalim")
    gca().set_aspect("equal", "box")

    # Set axis limits
    xlim!=nothing && plt.xlim(xlim) 
    ylim!=nothing && plt.ylim(ylim) 

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

        if length(legendlabels)==2
            legendlabels = legendlabels[idxs]
        end

        plt.legend(lines, legendlabels,
               loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), 
               edgecolor="k", ncol=2, 
              )
    end
    

    ax = plt.axes()
    ax.xaxis.set_tick_params(width=0.3)
    ax.yaxis.set_tick_params(width=0.3)
    ax.xaxis.set_tick_params(size=2.5)
    ax.yaxis.set_tick_params(size=2.5)
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
    end

end

