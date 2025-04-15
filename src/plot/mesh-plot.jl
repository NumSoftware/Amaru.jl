# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


MeshPlot_params = [
    FunInfo(:MeshPlot, "Creates a customizable `MeshPlot` instance used to plot finite element meshes."),
    ArgInfo(:mesh, "Finite element mesh or model to plot"),
    KwArgInfo((:size, :figsize), "Mesh drawing size in dpi", (220,150), length=2),
    KwArgInfo(:facecolor, "Face color", :aliceblue),
    KwArgInfo(:warp, "Warping scale", 0.0 ),
    KwArgInfo((:lw, :lineweight), "Edge weight", 0.4,  cond=:(lw>0) ),
    KwArgInfo(:field, "Scalar field", "" ),
    KwArgInfo(:limits, "Limits for the scalar field", [0.0,0.0], length=2 ),
    KwArgInfo(:mult, "Field multiplier", 1.0),
    KwArgInfo(:label, "Colorbar label", "", type=AbstractString ),
    KwArgInfo(:colormap, "Colormap for field display", :coolwarm),
    KwArgInfo((:colorbarloc,:colorbar), "Colorbar location", :right, values=(:none, :right, :bottom) ),
    KwArgInfo((:colorbarscale, :cbscale), "Colorbar scale", 0.9, cond=:(colorbarscale>0) ),
    KwArgInfo((:label, :colorbarlabel, :cblabel, :colorbartitle), "Colorbar label", "" ),
    KwArgInfo((:fontsize, :colorbarfontsize, :cbfontsize), "Colorbar font size", 7.0, cond=:(fontsize>0)),
    KwArgInfo(:gradientmode, "Sets the gradient mode for surfaces", :linear, values=(:constant,:linear,:nonlinear)),
    KwArgInfo(:font, "Font name", "NewComputerModern", type=AbstractString),
    KwArgInfo(:azimut, "Azimut angle for 3d in degrees", 30 ),
    KwArgInfo(:elevation, "Elevation angle for 3d in degrees", 30 ),
    KwArgInfo(:distance, "Distance from camera in 3d", 0.0, cond=:(distance>=0) ),
    KwArgInfo(:featurelines, "Flag to enhance feature lines", true, type=Bool ),
    # KwArgInfo(:wireframe, "Flag to show a wireframe", false, type=Bool ),
    KwArgInfo(:viewmode, "Type of view", :surface, values=(:surface, :wireframe, :outline) ),
    KwArgInfo((:lightvector, :lv), "Light direction vector", [0.0,0.0,0.0], length=3 )
]
@doc docstring(MeshPlot_params) MeshPlot(mesh; kwargs...)


mutable struct MeshPlot<:Figure
    mesh::AbstractDomain
    canvas::Union{FigureComponent, Nothing}
    colorbar::Union{FigureComponent, Nothing}
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    values::Vector{Float64}
    outerpad::Float64
    shades::Vector{Float64}
    args::NamedTuple
    azimut::Float64
    elevation::Float64
    distance::Float64
    
    gradientmode::Symbol
    lightvector::Vector{Float64}

    width::Float64
    height::Float64

    lw::Float64
    facecolor::Tuple
    field::String
    limits::Vector{Float64}
    mult::Float64
    warp::Float64
    label::AbstractString
    colormap::Colormap
    colorbarloc::Symbol
    colorbarscale::Float64
    colorbarfontsize::Float64
    featurelines::Bool
    # wireframe::Bool
    viewmode::Symbol

    function MeshPlot(mesh; args...)
        args = checkargs(args, MeshPlot_params, aliens=false)

        this = new()
        this.mesh = mesh
        this.canvas = nothing
        this.colorbar = nothing
        this.nodes = []
        this.elems = []
        this.values = []
        this.shades = []
        this.outerpad = 0.0

        this.width     = args.size[1]
        this.height    = args.size[2]
        this.lw        = args.lw
        this.facecolor = _colors_dict[args.facecolor]
        this.field     = string(args.field)
        this.limits    = args.limits
        this.mult      = args.mult
        this.warp      = args.warp
        this.label     = args.label

        colormap          = args.colormap isa Symbol ? Colormap(args.colormap) : args.colormap
        this.colormap     = colormap
        this.gradientmode = args.gradientmode

        this.azimut      = args.azimut
        this.elevation   = args.elevation
        this.distance    = args.distance
        this.featurelines = (args.featurelines || args.viewmode==:outline)
        this.viewmode    = args.viewmode
        this.lightvector = args.lightvector

        this.args = args

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

    return p1, cp2, cp3, p4
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
        x = node.coord[1]
        y′ = node.coord[2]*focal_length/(distance-x)
        z′ = node.coord[3]*focal_length/(distance-x)
        node.coord = Vec3(y′, z′, distance-x)
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

    if ndim==2
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]
        feature_edges = mplot.featurelines ? get_feature_edges(areacells) : Cell[]

        elems    = [ areacells; linecells; feature_edges ]
        nodes    = getnodes(elems)
    else
        # get surface cells and update
        volcells  = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==3 ]
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==2 ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]
        surfcells = get_outer_facets(volcells)
        feature_edges = mplot.featurelines ? [ get_feature_edges(surfcells); get_outer_facets(areacells) ] : Cell[]

        for c in surfcells
            c.id = c.owner.id # useful to get values but not shades
        end

        # detach nodes
        for edge in feature_edges
            edge.nodes = copy.(edge.nodes) 
        end

        tag!(feature_edges, "_feature_edge")
        elems    = [ surfcells; areacells; linecells; feature_edges ]
        nodes    = getnodes(elems) # todo: check for joined nodes in cracks
    end

    node_data = mesh.node_data
    elem_data = mesh.elem_data

    # Change coords if warping
    if mplot.warp>0.0
        if haskey(node_data, "U")
            U = node_data["U"]
            for node in nodes
                node.coord = node.coord + mplot.warp*U[node.id,:]  
            end
        else
            alert("MeshPlot: Vector field U not found for warping.")
        end
    end

    # 3D -> 2D projection
    if mesh.ctx.ndim==3
        # compute shades (before 2d projection)
        V = Vec3( cosd(mplot.elevation)*cosd(mplot.azimut), cosd(mplot.elevation)*sind(mplot.azimut), sind(mplot.elevation) ) # observer vector
        norm(mplot.lightvector)==0 && (mplot.lightvector = V)
        L = mplot.lightvector
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

        # raise featurelines cells
        out_nodes = getnodes(feature_edges)
        for node in out_nodes
            node.coord = node.coord - Vec3(0,0,0.01*(zmax-zmin))
        end

        # distances = [ sum(node.coord[3] for node in elem.nodes)/length(elem.nodes)  for elem in mesh.elems ]
        # distances = [ minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        distances = [ 0.9*sum(node.coord[3] for node in elem.nodes)/length(elem.nodes) + 0.1*minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        perm = sortperm(distances, rev=true)
        elems = elems[perm]
        mplot.shades = mplot.shades[perm]
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
            fvals = elem_data[field].*mplot.mult
            fmax = maximum(fvals)
            fmin = minimum(fvals)
        end
        if haskey(node_data, field)
            fvals = node_data[field].*mplot.mult
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
        end
        # mplot.colormap = resize(mplot.colormap, fmin, fmax, divergefromzero=mplot.args.divergefromzero)
        mplot.colormap = resize(mplot.colormap, fmin, fmax)

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
            location   = mplot.args.colorbarloc,
            label      = mplot.label,
            scale      = mplot.args.colorbarscale,
            colormap   = mplot.colormap,
            limits     = mplot.limits,
            fontsize   = mplot.args.fontsize,
            font       = mplot.args.font,
        )
        configure!(mplot, mplot.colorbar)

        if has_field
            if mplot.colorbar.location==:right
                rpane = mplot.colorbar.width
            elseif mplot.colorbar.location==:bottom
                bpane = mplot.colorbar.height
            end
        end
    end

    # Canvas box
    canvas = mplot.canvas

    canvas.width = width - rpane - 2*mplot.outerpad
    canvas.height = height - bpane - 2*mplot.outerpad
    canvas.box = [ mplot.outerpad, mplot.outerpad, width-rpane-mplot.outerpad, height-bpane-mplot.outerpad ]

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


function draw!(mplot::MeshPlot, cc::CairoContext)
    set_line_join(cc, Cairo.CAIRO_LINE_JOIN_ROUND)
 
    xmin, xmax = extrema( node.coord[1] for node in mplot.nodes)
    ymin, ymax = extrema( node.coord[2] for node in mplot.nodes)

    Xmin, Ymin, Xmax, Ymax = mplot.canvas.box
    has_field = mplot.field != ""
    is_nodal_field = has_field && haskey(mplot.mesh.node_data, mplot.field)

    ratio = min((Xmax-Xmin)/(xmax-xmin), (Ymax-Ymin)/(ymax-ymin) )
    dx = 0.5*((Xmax-Xmin)-ratio*(xmax-xmin))
    dy = 0.5*((Ymax-Ymin)-ratio*(ymax-ymin))
    set_matrix(cc, CairoMatrix([ratio, 0, 0, -ratio, Xmin+dx-xmin*ratio, Ymax-dy+ymin*ratio]...))

    # gray = sum(mplot.facecolor)/3*Vec3(0.5,0.5,0.5)

    # Draw elements
    for (i,elem) in enumerate(mplot.elems)

        if elem.tag=="_feature_edge"
            # gray = sum(mplot.facecolor)/3*Vec3(0.5,0.5,0.5)
            if has_field 
                if is_nodal_field
                    id = elem.nodes[1].id
                    color = mplot.colormap(mplot.values[id]).*0.45
                else
                    color = mplot.colormap(mplot.values[elem.owner.id]).*0.45
                end
                # color = Vec3(0.3,0.3,0.3) .+ mplot.colormap(mplot.values[id]).*0.3
                # color = sum(mplot.colormap(mplot.values[id]))/3*Vec3(0.5,0.5,0.5)
            else
                color = mplot.facecolor.*0.45
            end
            # color = _colors_dict[:red]

            x, y = elem.nodes[1].coord
            move_to(cc, x, y)
            # color = Vec3(0.4, 0.4, 0.4)
            pts = bezier_points(elem)
            curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(cc, 1.4*mplot.lw)
            set_source_rgb(cc, color...) 
            stroke(cc)
            continue
        end

        if elem.shape.family==LINECELL
            x, y = elem.nodes[1].coord
            move_to(cc, x, y)
            color = Vec3(0.8, 0.2, 0.1)
            if has_field && !is_nodal_field
                color = mplot.colormap(mplot.values[elem.id])
            end

            pts = bezier_points(elem)
            curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(cc, 2*mplot.lw)
            set_source_rgb(cc, color...) # gray
            stroke(cc)
            continue
        end

        # culling back faces
        if elem.owner !== nothing
            if elem.owner.shape.ndim==3
                !iscounterclockwise( elem.nodes ) && continue
            end
        end

        shade = mplot.mesh.ctx.ndim==3 ? mplot.shades[i] : 1.0
        draw_surface_cell!(cc, mplot, elem, has_field, shade)
    end

    # draw colorbar
    has_field && draw!(mplot, cc, mplot.colorbar)
end


function draw_surface_cell!(cc::CairoContext, mplot::MeshPlot, elem::AbstractCell, has_field::Bool, shade::Float64)
    is_nodal_field = has_field && haskey(mplot.mesh.node_data, mplot.field)

    # draw cells face
    edges = getedges(elem)
    if mplot.viewmode==:surface

        if !has_field || !is_nodal_field || mplot.gradientmode==:constant
            if has_field
                if is_nodal_field
                    val = sum( mplot.values[node.id] for node in elem.nodes)/length(elem.nodes)
                    color = mplot.colormap(val)
                else
                    color = mplot.colormap(mplot.values[elem.id])
                end
            else
                color = mplot.facecolor
            end

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(cc)
            move_to(cc, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(cc)
            set_source_rgb(cc, (color.*shade)...)
            fill(cc)

        elseif mplot.gradientmode==:linear
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

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(cc)
            move_to(cc, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(cc)
            set_source(cc, pat)
            fill(cc)
        else
            pattern = CairoPatternMesh()
            mesh_pattern_begin_patch(pattern)

            x, y = edges[1].nodes[1].coord
            mesh_pattern_move_to(pattern, x, y)
            nedges = length(edges)

            # draw elements
            color = mplot.facecolor
            if has_field && !is_nodal_field
                color = mplot.colormap(mplot.values[elem.id])
            end

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

            # set nodal colors
            for (i,node) in enumerate(elem.nodes[1:nedges])
                if has_field && is_nodal_field
                    color = mplot.colormap(mplot.values[node.id])
                end
                scolor = color.*shade # apply shade
                mesh_pattern_set_corner_color_rgb(pattern, i-1, scolor...)
            end

            mesh_pattern_end_patch(pattern)
            set_source(cc, pattern)
            paint(cc)
        end
    end

    # draw cell edges for viewmode=:surface and viewmode=:wireframe
    if mplot.viewmode != :outline
        x, y = edges[1].nodes[1].coord
        move_to(cc, x, y)
        color = mplot.facecolor.*0.45
        for edge in edges
            pts = bezier_points(edge)
            x, y = edge.nodes[1].coord
            move_to(cc, x, y)
            curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            id = edge.nodes[1].id
            if has_field 
                if is_nodal_field
                    color = mplot.colormap(mplot.values[id]).*0.55
                else
                    # color = mplot.facecolor
                    color = mplot.colormap(mplot.values[elem.id]).*0.55
                end
            end
            set_line_width(cc, mplot.lw)
            set_source_rgb(cc, color...)
            stroke(cc)
        end
    end
end

