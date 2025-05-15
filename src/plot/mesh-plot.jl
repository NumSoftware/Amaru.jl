# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


MeshPlot_params = [
    FunInfo(:MeshPlot, "Creates a customizable `MeshPlot` instance used to plot finite element meshes."),
    ArgInfo(:mesh, "Finite element mesh or model to plot"),
    KwArgInfo(:size, "Figure size in dpi", (220,150), length=2), # (400,250)
    KwArgInfo(:face_color, "Face color", :aliceblue),
    KwArgInfo(:warp, "Warping scale", 0.0 ),
    KwArgInfo((:lw, :line_weight), "Edge weight", 0.3,  cond=:(lw>=0) ),
    KwArgInfo(:field, "Scalar field", "" ),
    KwArgInfo(:limits, "Limits for the scalar field", [0.0,0.0], length=2 ),
    KwArgInfo(:field_factor, "Field multiplier", 1.0),
    KwArgInfo(:label, "Colorbar label", "", type=AbstractString ),
    KwArgInfo(:colormap, "Colormap for field display", :coolwarm),
    KwArgInfo(:diverging, "Modifies colormap to diverge from zero", false, type=Bool ),
    KwArgInfo(:colorbar, "Colorbar location", :right, values=(:none, :right, :bottom) ),
    KwArgInfo(:colorbar_scale, "Colorbar scale", 0.9, cond=:(colorbar_scale>0) ),
    KwArgInfo((:label, :colorbar_label), "Colorbar label", "" ),
    KwArgInfo(:bins, "Number of bins in the colorbar", 6 ),
    KwArgInfo(:font, "Font name", "NewComputerModern", type=AbstractString),
    KwArgInfo(:font_size, "Colorbar font size", 7.0, cond=:(font_size>0)),
    KwArgInfo(:gradient_mode, "Sets the gradient mode for surfaces", :linear, values=(:constant,:linear,:nonlinear)),
    KwArgInfo(:azimut, "Azimut angle for 3d in degrees", 30 ),
    KwArgInfo(:elevation, "Elevation angle for 3d in degrees", 30 ),
    KwArgInfo(:distance, "Distance from camera in 3d", 0.0, cond=:(distance>=0) ),
    KwArgInfo(:feature_edges, "Flag to enhance feature lines", true, type=Bool ),
    KwArgInfo(:view_mode, "Type of view", :surface_with_edges, values=(:surface_with_edges, :surface, :wireframe, :outline) ),
    KwArgInfo(:light_vector, "Light direction vector", [0.0,0.0,0.0], length=3 ),
    KwArgInfo(:node_labels, "Flag to show node labels", false, type=Bool ),
    KwArgInfo(:axes, "Flag to show the coordinate axes", :none, type=Symbol, values=_axes_widget_locations ),
    KwArgInfo(:axis_labels, "Custom axes labels", String[], type=Array ),
    KwArgInfo(:quiet, "Flag to set silent mode", false),
]
@doc docstring(MeshPlot_params) MeshPlot(mesh; kwargs...)


mutable struct MeshPlot<:Figure
    mesh::AbstractDomain
    canvas::Union{FigureComponent, Nothing}
    colorbar::Union{FigureComponent, Nothing}
    axes::Union{FigureComponent, Nothing}
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    feature_edges_d::Dict
    values::Vector{Float64}
    outerpad::Float64
    shades::Vector{Float64}
    azimut::Float64
    elevation::Float64
    distance::Float64
    
    gradient_mode::Symbol
    light_vector::Vector{Float64}

    width::Float64
    height::Float64

    lw::Float64
    outline_lw::Float64
    face_color::Tuple
    field::String
    limits::Vector{Float64}
    field_factor::Float64
    warp::Float64
    label::AbstractString
    colormap::Colormap
    diverging::Bool
    colorbar_loc::Symbol
    colorbar_scale::Float64
    bins::Int
    font::String
    font_size::Float64
    show_feature_edges::Bool
    view_mode::Symbol
    node_labels::Bool
    axes_loc::Symbol
    axis_labels::Vector{AbstractString}
    quiet::Bool

    # @params function MeshPlot(mesh =>, 
    #                             size::Tuple                     => "sdfsaldf sadf";
    #                             face_color::Symbol = :aliceblue => ("fsdf sdasdf sd fsadf", face_color>2),
    #                             face_colo::Symbol = :aliceblue  => ("asdf asdf sa fasdf sdf", (:a, :b))
    # ) "Creates a customizable `MeshPlot` instance used to plot finite element meshes."

    # end
    function MeshPlot(mesh; args...)
        args = checkargs(args, MeshPlot_params, aliens=false)

        this = new()
        this.mesh = mesh
        this.canvas = nothing
        this.colorbar = nothing
        this.axes = nothing
        this.nodes = []
        this.elems = []
        this.feature_edges_d = Dict()
        this.values = []
        this.shades = []
        this.outerpad = 0.0

        this.width        = args.size[1]
        this.height       = args.size[2]
        this.lw           = args.lw
        this.outline_lw   = 1.4*args.lw
        this.face_color   = _colors_dict[args.face_color]
        this.field        = string(args.field)
        this.limits       = vec(args.limits)
        this.field_factor = args.field_factor
        this.warp         = args.warp
        this.label        = args.label

        colormap            = args.colormap isa Symbol ? Colormap(args.colormap) : args.colormap
        this.colormap       = colormap
        this.diverging      = args.diverging
        this.colorbar_loc   = args.colorbar
        this.colorbar_scale = args.colorbar_scale
        this.bins           = args.bins
        this.gradient_mode  = args.gradient_mode

        this.font               = args.font
        this.font_size          = args.font_size
        this.azimut             = args.azimut
        this.elevation          = args.elevation
        this.distance           = args.distance
        this.show_feature_edges = (args.feature_edges || args.view_mode==:outline)
        this.view_mode          = args.view_mode
        this.light_vector       = args.light_vector

        this.node_labels = args.node_labels
        this.axes_loc    = args.axes
        this.axis_labels = length(args.axis_labels)==0 ? [L"x", L"y", L"z"] : args.axis_labels

        this.quiet = args.quiet

        if !this.quiet 
            printstyled("Mesh/Model plot\n", bold=true, color=:cyan)
            println("  size: $(this.width) x $(this.height) dpi")
            println("  view mode: $(this.view_mode)")
            this.field != "" && println("  field: $(this.field)")
        end

        return this
    end
end


function bezier_points(edge::AbstractCell)
    p1 = edge.nodes[1].coord[1:2]
    p4 = edge.nodes[2].coord[1:2]
    ξ1 = -1/3
    ξ2 = +1/3
    C = getcoords(edge.nodes, 2)
    p2 = C'*edge.shape.func([ξ1])
    p3 = C'*edge.shape.func([ξ2])

    cp2 = 1/6*(-5*p1+18*p2- 9*p3+2*p4)
    cp3 = 1/6*( 2*p1- 9*p2+18*p3-5*p4)

    return [p1, cp2, cp3, p4]
end


function project_to_2d!(nodes, azimut, elevation, distance)
    # Find bounding box
    xmin, xmax = extrema( node.coord[1] for node in nodes)
    ymin, ymax = extrema( node.coord[2] for node in nodes)
    zmin, zmax = extrema( node.coord[3] for node in nodes)
    reflength = max(xmax-xmin, ymax-ymin, zmax-zmin)

    # Centralize 
    center = 0.5*Vec3(xmin+xmax, ymin+ymax, zmin+zmax)
    for node in nodes
        node.coord = node.coord - center
    end

    # Rotation around z axis
    θ = -azimut*pi/180
    R = Quaternion(cos(θ/2), 0, 0, sin(θ/2))
    for node in nodes
        node.coord = (R*node.coord*conj(R))[2:4]
    end

    # Rotation around y axis
    θ = elevation*pi/180
    R = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
    for node in nodes
        node.coord = (R*node.coord*conj(R))[2:4]
    end

    # Set projection values
    distance==0 && (distance=reflength*3)
    distance = max(distance, reflength)
    focal_length = 0.1*distance

    # Make projection
    for node in nodes
        x´ = node.coord[1]
        y′ = node.coord[2]*focal_length/(distance-x´)
        z′ = node.coord[3]*focal_length/(distance-x´)
        node.coord = Vec3(y′, z′, distance-x´)
    end

    xmin, xmax = extrema( node.coord[1] for node in nodes)
    ymin, ymax = extrema( node.coord[2] for node in nodes)

    # normalize
    l = max(xmax-xmin, ymax-ymin)
    for node in nodes
        x = (node.coord[1]-xmin)/l
        y = (node.coord[2]-ymin)/l
        node.coord = Vec3(x, y, node.coord[3])
    end
end


function configure!(mplot::MeshPlot)
    mesh = copy(mplot.mesh)
    ndim = mesh.ctx.ndim

    # Change coords if warping
    if mplot.warp>0.0
        U = get(mesh.node_data, "U", nothing)
        if U === nothing
            alert("MeshPlot: Vector field U not found for warping.")
        else
            for node in mesh.nodes
                node.coord = node.coord + mplot.warp*U[node.id,:]  
            end
        end
    end

    if ndim==2
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]
        feature_edges = get_outer_facets(areacells)

        elems    = [ areacells; linecells ]
        nodes    = getnodes(elems)
    else
        # get surface cells and update
        volcells  = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==3 ]
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==2 ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]
        surfcells = get_outer_facets(volcells)
        feature_edges = get_feature_edges(surfcells)
        
        for c in surfcells
            c.id = c.owner.id # useful to get values but not shades
        end

        elems    = [ surfcells; areacells; linecells ]

        # unique by id
        nodes = unique( n -> n.id, [ node for elem in elems for node in elem.nodes ] )
    end

    node_data = mesh.node_data
    elem_data = mesh.elem_data

    # populate feature_edges_d
    for edge in feature_edges
        node_idxs = sort([ node.id for node in edge.nodes ])
        mplot.feature_edges_d[node_idxs] = edge
    end

    # 3D -> 2D projection
    if mesh.ctx.ndim==3
        # compute shades (before 2d projection)
        V = Vec3( cosd(mplot.elevation)*cosd(mplot.azimut), cosd(mplot.elevation)*sind(mplot.azimut), sind(mplot.elevation) ) # observer vector
        norm(mplot.light_vector)==0 && (mplot.light_vector = V)
        L = mplot.light_vector
        # shades = zeros(length(mesh.elems))
        shades = zeros(length(elems))
        for (i,elem) in enumerate(elems)
            elem.shape.family==BULKCELL || continue
            N = get_facet_normal(elem)
            dot(N,L)<0 && (N = -N)
            R = normalize(2*N*dot(L,N) - L)  # reflection vector
            dot(V,R)<0 && (R = -R)

            Ia = 0.68 # ambient
            Id = 0.34 # diffuse
            Is = 0.20 # specular
            sh = 2.00 # shininess

            shades[i] = Ia + Id*max(0,dot(L,N)) + Is*max(0,dot(V,R))^sh
        end
        mplot.shades = shades

        # compute projection
        project_to_2d!(nodes, mplot.azimut, mplot.elevation, mplot.distance)
        zmin, zmax = extrema(node.coord[3] for node in nodes)

        # distances = [ sum(node.coord[3] for node in elem.nodes)/length(elem.nodes)  for elem in mesh.elems ]
        # distances = [ minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        distances = [ 0.9*sum(node.coord[3] for node in elem.nodes)/length(elem.nodes) + 0.1*minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        perm = sortperm(distances, rev=true)
        elems = elems[perm]
        mplot.shades = mplot.shades[perm]
    else
        # sort elems by area (show largest first)
        areas = [ cell_extent(elem) for elem in elems ]
        perm = sortperm(areas, rev=true)
        elems = elems[perm]
    end

    # Field 
    has_field = mplot.field != ""

    if has_field
        field = string(mplot.field)
        mplot.label == ""  && (mplot.label = field)

        # available fields
        fields = [ string(field) for field in Iterators.flatten([keys(node_data), keys(elem_data)]) ]
        field in fields || error("mplot: field $field not found. Available fields are: $(join(fields, ", ", " and "))")
        
        if haskey(elem_data, field)
            fvals = elem_data[field].*mplot.field_factor
            fmax = maximum(fvals)
            fmin = minimum(fvals)
        end
        if haskey(node_data, field)
            fvals = node_data[field].*mplot.field_factor
            fmax = maximum(fvals[node.id] for node in nodes)
            fmin = minimum(fvals[node.id] for node in nodes)
        end
        if fmin==fmax
            fmin -= 1
            fmax += 1
        end

        # Colormap
        mplot.values = fvals
        if mplot.limits==[0.0, 0.0]
            mplot.limits = [fmin, fmax]
        else
            fmin, fmax = mplot.limits
        end
        mplot.colormap = resize(mplot.colormap, fmin, fmax, diverging=mplot.diverging)

    end

    mplot.nodes = nodes
    mplot.elems = elems

    # Canvas 
    mplot.canvas = Canvas()
    width, height = mplot.width, mplot.height
    mplot.outerpad = 0.01*min(width, height)

    rpane = 0.0
    bpane = 0.0
    
    # Colorbar
    if has_field
        mplot.colorbar = Colorbar(;
            location  = mplot.colorbar_loc,
            label     = mplot.label,
            scale     = mplot.colorbar_scale,
            colormap  = mplot.colormap,
            limits    = mplot.limits,
            font_size = mplot.font_size,
            font      = mplot.font,
            bins      = mplot.bins,
        )
        configure!(mplot, mplot.colorbar)

        if mplot.colorbar.location==:right
            rpane = mplot.colorbar.width
        elseif mplot.colorbar.location==:bottom
            bpane = mplot.colorbar.height
        end
    end

    # Canvas box
    canvas = mplot.canvas
    
    canvas.width = width - rpane - 2*mplot.outerpad
    canvas.height = height - bpane - 2*mplot.outerpad
    Xmin = mplot.outerpad
    Ymin = mplot.outerpad
    Xmax = width - rpane - mplot.outerpad
    Ymax = height - bpane - mplot.outerpad
    canvas.box = [ Xmin, Ymin, Xmax, Ymax ] # in user coordinates

    xmin, xmax = extrema( node.coord[1] for node in mplot.nodes)
    ymin, ymax = extrema( node.coord[2] for node in mplot.nodes)

    ratio = min((Xmax-Xmin)/(xmax-xmin), (Ymax-Ymin)/(ymax-ymin) )
    dX = 0.5*((Xmax-Xmin)-ratio*(xmax-xmin))
    dY = 0.5*((Ymax-Ymin)-ratio*(ymax-ymin))

    xmin = xmin - dX/ratio
    xmax = xmax + dX/ratio
    ymin = ymin - dY/ratio
    ymax = ymax + dY/ratio

    mplot.canvas.limits = [xmin, ymin, xmax, ymax] # in data coordinates

    # Axes widget
    if mplot.axes_loc != :none
        mplot.axes = AxisWidget(
            location  = mplot.axes_loc,
            font_size = mplot.font_size,
            font      = mplot.font,
            azimut    = mplot.azimut,
            elevation = mplot.elevation,
            labels    = mplot.axis_labels[1:ndim],
        )
        configure!(mplot.axes)
    end
end


function iscounterclockwise(points::Array{Node,1})
    val = 0.0
    n = length(points)
    for i in 1:n
        x1, y1 = points[i].coord
        x2, y2 = points[mod1(i+1, n)].coord
        val += (x1*y2) - (x2*y1)
    end
    return val > 0
end


function draw!(mplot::MeshPlot, ctx::CairoContext)
    # Cairo.push_group(ctx)
    # set_operator(ctx, Cairo.OPERATOR_SOURCE)

    set_line_join(ctx, Cairo.CAIRO_LINE_JOIN_ROUND)
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)
    font = get_font(mplot.font)
    select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(ctx, mplot.font_size)
 
    has_field = mplot.field != ""
    is_nodal_field = has_field && haskey(mplot.mesh.node_data, mplot.field)

    Xmin, Ymin, Xmax, Ymax = mplot.canvas.box
    xmin, ymin, xmax, ymax = mplot.canvas.limits

    ratio = (Xmax-Xmin)/(xmax-xmin) # == (Ymax-Ymin)/(ymax-ymin)
    set_matrix(ctx, CairoMatrix([ratio, 0, 0, -ratio, Xmin-xmin*ratio, Ymax+ymin*ratio]...))

    # Draw elements
    for (i,elem) in enumerate(mplot.elems)
        
        elem.shape.family in ( LINECELL, BULKCELL ) || continue

        # culling back faces
        if elem.owner !== nothing && mplot.view_mode != :outline
            if elem.owner.shape.ndim==3
                !iscounterclockwise( elem.nodes ) && continue
            end
        end

        if elem.shape.family==LINECELL
            x, y = elem.nodes[1].coord
            move_to(ctx, x, y)
            color = Vec3(0.8, 0.2, 0.1)
            if has_field && !is_nodal_field
                color = mplot.colormap(mplot.values[elem.id])
            end

            pts = bezier_points(elem)
            curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(ctx, 2*mplot.lw)
            set_source_rgb(ctx, color...) # gray
            stroke(ctx)
        else    
            shade = mplot.mesh.ctx.ndim==3 ? mplot.shades[i] : 1.0
            draw_surface_cell!(ctx, mplot, elem, has_field, shade)
        end

        # point labels
        Cairo.save(ctx)
        set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
        if mplot.node_labels 
            for node in elem.nodes
                x, y = data2user(mplot.canvas, node.coord[1], node.coord[2])
                set_source_rgb(ctx, _colors_dict[:blue]...)
                draw_text(ctx, x, y, string(node.id), halign="center", valign="center", angle=0)
            end
        end
        Cairo.restore(ctx)
    end

    # grp = Cairo.pop_group(ctx)
    # Cairo.set_source(ctx, grp)

    # Composite all at once
    # set_operator(ctx, Cairo.OPERATOR_OVER)
    # paint(ctx)

    # draw colorbar
    has_field && draw!(mplot, ctx, mplot.colorbar)

    # draw axes
    if mplot.axes_loc != :none
        if mplot.axes_loc == :bottomleft
            x = mplot.outerpad
            y = mplot.canvas.box[4] - mplot.axes.height
        elseif mplot.axes_loc == :topleft
            x = mplot.outerpad
            y = mplot.canvas.box[2]
        else
            error("MeshPlot: axes location $(mplot.axes_loc) not implemented")
        end
        move_to(ctx, x, y)
        draw!(ctx, mplot.axes)
    end

end




function draw_surface_cell!(ctx::CairoContext, mplot::MeshPlot, elem::AbstractCell, has_field::Bool, shade::Float64)
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)

    is_nodal_field = has_field && haskey(mplot.mesh.node_data, mplot.field)

    constant_field = !is_nodal_field || mplot.gradient_mode==:constant

    # compute linear gradients
    if constant_field
        if has_field
            if is_nodal_field
                val = sum( mplot.values[node.id] for node in elem.nodes)/length(elem.nodes)
                color = mplot.colormap(val)
            else
                color = mplot.colormap(mplot.values[elem.id])
            end
        else
            color = mplot.face_color
        end

        # edges
        factor = 0.55
        if mplot.view_mode == :surface # disguise edges in surface view
            factor = 1.0
        end
        edge_color    = color.*factor.*shade
        outline_color = color.*0.55.*shade
    else
        cm = mplot.colormap
        npoints = elem.shape.basic_shape.npoints # use only corner nodes
        nodes = elem.nodes[1:npoints]

        # regression plane
        values = [ mplot.values[node.id] for node in nodes ]
        X = [ elem.nodes[i].coord[j] for i in 1:npoints, j=1:3]
        X[:,3] .= 1.0 # for regression
        A = pinv(X)*values # regression plane coefficients
        V = normalize(A[1:2]) # gradient direction
        
        # distances along V and sorting
        d = [ dot(V, node.coord[1:2]) for node in nodes ] 
        idxs = sortperm(d)
        nodes = nodes[idxs]
        values = values[idxs]

        # compute gradient xtremes
        (xa, ya) = nodes[1].coord
        (xb, yb) = nodes[end].coord
        (xb, yb) = [xa, ya] + dot([xb - xa, yb - ya], V )*V
        
        vmin = dot(A, (xa, ya, 1.0))
        vmax = dot(A, (xb, yb, 1.0))

        # compute stops
        stops = Float64[]
        for i in 1:npoints
            val  = values[i]
            stop = clamp(round((val-vmin)/(vmax-vmin), digits=8), 0.0, 1.0) # clamp is required
            push!(stops, stop)
        end

        pat = pattern_create_linear(xa, ya, xb, yb)
        for (val, stop) in zip(values, stops)
            color = cm(val).*shade
            pattern_add_color_stop_rgb(pat, stop, color...)
        end
        
        # edges
        factor = 0.55
        if mplot.view_mode == :surface # disguise edges in surface view
            factor = 1.0
            mplot.lw = 0.
        end
        
        edge_pat = pattern_create_linear(xa, ya, xb, yb)
        for (val, stop) in zip(values, stops)
            color = cm(val).*shade*factor
            pattern_add_color_stop_rgb(edge_pat, stop, color...)
        end

        # outline
        outline_pat = pattern_create_linear(xa, ya, xb, yb)
        for (val, stop) in zip(values, stops)
            color = cm(val).*shade*0.55
            pattern_add_color_stop_rgb(outline_pat, stop, color...)
        end
    end


    # draw cells face
    center       = sum(getcoords(elem), dims=1)[1:2]/length(elem.nodes)
    edges        = getedges(elem)
    show_surface = mplot.view_mode in (:surface, :surface_with_edges)
    show_edges   = mplot.view_mode in (:surface, :surface_with_edges, :wireframe) # edges are disguised in surface view
    if mplot.view_mode != :outline

        if constant_field
            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(ctx)

            if show_surface
                set_source_rgb(ctx, color...)
                fill_preserve(ctx)
            end

            if show_edges
                set_source_rgb(ctx, edge_color...)
                set_line_width(ctx, mplot.lw)
                stroke(ctx)
            end
        elseif mplot.gradient_mode==:linear

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                # mplot.view_mode == :surface && expand_points!(center, pts)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(ctx)

            if show_surface
                set_source(ctx, pat)
                fill_preserve(ctx)
            end

            if show_edges
                set_source(ctx, edge_pat)
                set_line_width(ctx, mplot.lw)
                stroke(ctx)
            end
        else # nonlinear
            # set pattern mesh for nonlinear gradient
            pattern = CairoPatternMesh()
            mesh_pattern_begin_patch(pattern)

            x, y = edges[1].nodes[1].coord
            mesh_pattern_move_to(pattern, x, y)
            nedges = length(edges)
            
            if length(edges)==3
                for edge in edges[1:2]
                    x, y = edge.nodes[2].coord
                    mesh_pattern_line_to(pattern, x, y)
                end
            else
                for edge in edges
                    pts = bezier_points(edge)
                    mesh_pattern_curve_to(pattern, pts[2]..., pts[3]..., pts[4]...)
                end
            end
            
            # elem colors
            color = mplot.face_color
            if has_field && !is_nodal_field
                color = mplot.colormap(mplot.values[elem.id])
            end

            # set nodal colors
            for (i,node) in enumerate(elem.nodes[1:nedges])
                if has_field && is_nodal_field
                    color = mplot.colormap(mplot.values[node.id])
                end
                scolor = color.*shade # apply shade
                mesh_pattern_set_corner_color_rgb(pattern, i-1, scolor...)
            end
            mesh_pattern_end_patch(pattern)

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            if show_surface
                set_source(ctx, pattern)
                fill_preserve(ctx)
                # paint(ctx) # only if using pattern instead of a path
            end

            if show_edges
                set_source(ctx, edge_pat)
                set_line_width(ctx, mplot.lw)
                stroke(ctx)
            end
        end
    end

    # draw feature edges
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)

    if mplot.show_feature_edges
        new_path(ctx) # clear path, e.g. when last command used preserve
        set_line_width(ctx, mplot.outline_lw)
        if constant_field
            set_source_rgb(ctx, outline_color...)
        else
            set_source(ctx, outline_pat)
        end

        for edge in edges
            node_idxs = sort([ node.id for node in edge.nodes ])
            haskey(mplot.feature_edges_d, node_idxs) || continue

            x, y = edge.nodes[1].coord
            move_to(ctx, x, y)
            pts = bezier_points(edge)
            curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            stroke(ctx)
        end
    end

end

