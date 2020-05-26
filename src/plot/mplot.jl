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
                #if α!=0
                    #@show α
                #end
                #@show angle
                angle<=α && push!(main_edges, edge)
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


import PyCall.PyObject # required

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

`alpha         = 1.0`   : Opacity,

`field         = nothing` : If provided, plots corresponding field

`fieldscale    = 1.0` : Factor multiplied to `field` values

`fieldlims     = ()` : Tuple `(min, max)` with field limits

`vectorfield   = nothing` : If provided, plots corresponding vector field

`arrowscale    = 0.0` : Factor multiplied to `vectorfield` values

`colormap      = nothing` : Colormap accrding to PyPlot

`colorbarscale = 0.9` : Scale of color bar

`colorbarlabel = ""` : Label at color bar

`colorbarlocation = ""` : Location of color bar (top, bottom, left and right)

`colorbarpad   = 0.0` : Separation of color bar from the plot

`warpscale     = 0.0` : Factor multiplied to "U" field when available

`highlightcell = 0` : Cell number to be highlighted

`elev          = 30.0` : 3D plot elevation

`azim          = 45.0` : 3D plot azimute

`dist          = 10.0` : 3D plot distance from observer

`outline       = true` : Highlight main edges of 3D meshes in the pdf output

`outlineangle  = 30` : Limit angle to identify main edges

`figsize       = (3,2.5)` : Figure size

`leaveopen     = false` : If true, leaves the plot open so other drawings can be added
"""
function mplot(
               mesh::AbstractMesh,
               filename::String = "";
               axis             = true,
               lw               = 0.3,
               nodemarkers     = false,
               nodelabels      = false,
               celllabels       = false,
               alpha            = 1.0,
               field            = nothing,
               fieldscale       = 1.0,
               fieldlims        = (),
               vectorfield      = nothing,
               arrowscale       = 0.0,
               opacity          = 1.0,
               colormap         = nothing,
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
               outlineangle     = 30,
               figsize          = (3,2.5),
               leaveopen        = false,
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
            node_data[field] = data[pt_ids]
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
        options = "axis, lw, nodemarkers, nodelabels, celllabels, alpha, field,
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
    @eval import PyPlot:plt, matplotlib, figure, art3D, Axes3D, ioff
    @eval ioff()

    # fix PyPlot
    @eval import PyPlot:getproperty, LazyPyModule
    if ! @eval hasmethod(getproperty, (LazyPyModule, AbstractString))
        @eval Base.getproperty(lm::LazyPyModule, s::AbstractString) = getproperty(PyObject(lm), s)
    end

    plt.close("all")

    plt.rc("font", family="serif", size=6)
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
        colormaps = matplotlib.pyplot.colormaps()
        if colormap in colormaps
            cmap = matplotlib.cm.get_cmap(colormap)
        else
            @error "mplot: Invalid colormap $colormap"
            error("  colormap should be one of:\n", colormaps)
        end
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
        fvals *= fieldscale
        fieldlims==() && (fieldlims = extrema(fvals))
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
            continue


            length(con)==3 && ( con = con[[1,3,2]] )
            X = XYZ[con, 1]
            Y = XYZ[con, 2]
            Z = XYZ[con, 3]
            color = "red"

            if has_field
                v = fvals[i]
                idx = (v-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                color = cmap(idx)
            end
            plt.plot(X, Y, Z, color=color, lw=1.0)
        end
        cltn = @eval art3D[:Line3DCollection]($all_verts, cmap=$cmap, lw=0.7, edgecolor="red")
        if has_field
            cltn.set_array(lfvals)
            cltn.set_clim(fieldlims)
        end
        ax.add_collection3d(cltn)

        # Plot main edges and surface cells
        all_verts  = []

        # Plot surface cells
        for i=1:ncells
            shape = cells[i].shape
            shape.family == SOLID_SHAPE || continue # only surface cells
            con = connect[i]
            points = [ XYZ[i,1:3] for i in con ]
            verts = plot_data_for_cell3d(points, shape)
            push!(all_verts, verts)
        end

        edgecolor = (0.4, 0.4, 0.4, 1.0)

        # Plot main edges
        filename == "" && (outline=false)
        if outline
            edges = get_main_edges(cells, outlineangle)
            θ, γ = (azim+0)*pi/180, elev*pi/180
            ΔX = [ cos(θ)*cos(γ), sin(θ)*cos(γ), sin(γ) ]*0.01*L

            for edge in edges
                #p1 = edge.nodes[1]
                #p2 = edge.nodes[2]
                #verts = [ p1.coord, p2.coord ]
                id1 = edge.nodes[1].id
                id2 = edge.nodes[2].id
                verts = [ XYZ[id_dict[id1],:], XYZ[id_dict[id2],:] ]
                for v in verts
                    v .+= ΔX
                end
                push!(all_verts, verts)
            end
            edgecolor = [ fill((0.4, 0.4, 0.4, 1.0), ncells) ; fill((0.20, 0.20, 0.20, 1.0), length(edges)) ]
        end

        facecolor = (0.94, 0.97, 1.0, opacity)
        cltn = @eval art3D[:Poly3DCollection]($all_verts, cmap=$cmap, facecolor=$facecolor, edgecolor=$edgecolor, lw=$lw, alpha=$alpha)

        if has_field
            if outline
                fvals = [ fvals; fill(mean(fvals), length(edges)) ]
            end

            cltn.set_array(fvals)
            cltn.set_clim(fieldlims)
        end
        ax.add_collection3d(cltn)

        # plot colorbar
        if has_field
            cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=10*colorbarscale*figsize[2], format="%.1f", pad=colorbarpad, location=colorbarlocation)
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
            X = XYZ[con, 1]
            Y = XYZ[con, 2]
            color = "red"

            if has_field
                v = fvals[i]
                idx = (v-fieldlims[1])/(fieldlims[2]-fieldlims[1])
                color = cmap(idx)
            end
            plt.plot(X, Y, color=color, lw=1.0)
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
        if has_field
            cltn.set_array(fvals)
            cltn.set_clim(fieldlims)
            h = colorbarorientation=="vertical" ? figsize[2] : figsize[1]
            h = norm(figsize)
            #cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=0.9*20*colorbarscale, format="%.1f", pad=colorbarpad, orientation=colorbarorientation)
            #cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=15, format="%.1f", pad=colorbarpad, orientation=colorbarorientation)
            cbar = plt.colorbar(cltn, label=colorbarlabel, shrink=colorbarscale, aspect=4*colorbarscale*h, format="%.1f", pad=colorbarpad, orientation=colorbarorientation)
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
            ax.scatter(X, Y, Z, color="black", marker="o", s=1)
        else
            plt.plot(X, Y, color="black", marker="o", markersize=3, lw=0)
        end
    end

    # Node markers on line cells
    ids = unique!([ id for i=1:ncells for id in connect[i] if cells[i].shape.family==LINE_SHAPE ])
    if length(ids)>0
        if ndim==3
            ax.scatter(X[ids], Y[ids], Z[ids], color="black", marker="o", s=1)
        else
            plt.plot(X[ids], Y[ids], color="black", marker="o", markersize=3, lw=0)
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
        verbose && printstyled("  file $filename saved\n", color=:cyan)
    end

    # Do not close if in IJulia
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return
    end

    leaveopen || plt.close("all")

    return

end
