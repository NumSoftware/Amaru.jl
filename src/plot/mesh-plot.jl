# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MeshPlot<:AbstractChart
    mesh::AbstractDomain
    canvas::Union{ChartComponent, Nothing}
    colorbar::Union{ChartComponent, Nothing}
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    values::Vector{Float64}
    outerpad::Float64
    args::NamedTuple

    figsize::Union{Array, Tuple}
    lw::Float64
    facecolor::Tuple
    field::String
    limits::Vector{Float64}
    warp::Float64
    label::AbstractString
    colormap::Colormap
    colorbarloc::Symbol
    colorbarscale::Float64
    colorbarfontsize::Float64
    azimut::Float64
    elevation::Float64
    distance::Float64
    outline::Bool

    function MeshPlot(mesh; args...)
        args = checkargs( args, 
            ArgInfo( :figsize, "Mesh drawing size in dpi", default=(300,200), length=2),
            ArgInfo( :facecolor, "Surface color", default=:aliceblue),
            ArgInfo( :warp, "Warping scale", default=1.0 ),
            ArgInfo( (:lw, :lineweight), "Line weight", default=0.5,  condition=:(lw>0) ),
            ArgInfo( :field, "Scalar field", default="" ),
            ArgInfo( :limits, "Field limit values", default=[0.0,0.0], length=2 ),
            ArgInfo( :label, "Colorbar label", default="", type=AbstractString ),
            ArgInfo( :colormap, "Colormap for field display", default=:coolwarm, values=_colormaps_list),
            ArgInfo( :divergefromzero, "Sets if colormap will diverge from zero", default=false, type=Bool),
            ArgInfo( (:colorbarloc,:colorbar), "Colorbar location", default=:right, values=(:right, :bottom) ),
            ArgInfo( (:colorbarscale, :cbscale), "Colorbar scale", default=0.9, condition=:(colorbarscale>0) ),
            ArgInfo( (:label, :colorbarlabel, :cblabel, :colorbartitle), "Colorbar label", default="" ),
            ArgInfo( (:fontsize, :colorbarfontsize, :cbfontsize), "Colorbar font size", default=9.0, condition=:(fontsize>0)),
            ArgInfo( :font, "Font name", default="NewComputerModern", type=AbstractString),
            ArgInfo( :azimut, "Azimut angle for 3d in degrees", default=30 ),
            ArgInfo( :elevation, "Elevation angle for 3d in degrees", default=30 ),
            ArgInfo( :distance, "Distance from camera in 3d", default=1.0, condition=:(distance>0) ),
            ArgInfo( :outline, "Flag to show outline", default=true, type=Bool ),
        )

        this = new()
        this.mesh = mesh
        this.canvas = nothing
        this.colorbar = nothing
        this.nodes = []
        this.elems = []
        this.values = []
        this.outerpad = 0.0

        this.figsize = args.figsize
        this.lw = args.lw
        this.facecolor = _colors_dict[args.facecolor]
        this.field = args.field
        this.limits = args.limits
        this.warp = args.warp
        this.label = args.label
        this.colormap = _colormaps_dict[args.colormap]
        # this.colorbarloc = args.colorbarloc
        # this.colorbarscale = args.colorbarscale
        # this.colorbarfontsize = args.colorbarfontsize
        this.azimut = args.azimut
        this.elevation = args.elevation
        this.distance = args.distance
        this.outline = args.outline
        this.args = args

        return this
    end

end


function bezier_points(edge)
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
    isnan(distance) && (distance=reflength*6)
    distance = max(distance, reflength)
    focal_length = 0.1*distance

    # Make projection
    for node in nodes
        x = node.coord[1]
        y′ = node.coord[2]*focal_length/(distance-x)
        z′ = node.coord[3]*focal_length/(distance-x)
        node.coord =  Vec3(y′, z′, distance-x)
    end
end


function configure!(mplot::MeshPlot)
    orig_mesh = mplot.mesh
    mesh = copy(mplot.mesh)
    ndim = mesh.env.ndim

    if ndim==2
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]

        elems = [ areacells; linecells ]
        nodes = getnodes(elems)
    else
        # get surface cells and update

        volcells  = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==3 ]
        areacells = [ elem for elem in mesh.elems.active if elem.shape.family==BULKCELL && elem.shape.ndim==2 ]
        linecells = [ cell for cell in mesh.elems.active if cell.shape.family==LINECELL]
        surfcells = get_surface(volcells)
        # outline_edges = outline ? copy.(get_outline_edges(surfcells)) : Cell[] # copies nodes
        outline_edges = mplot.outline ? get_outline_edges(surfcells) : Cell[] # copies nodes
        # detach nodes
        for edge in outline_edges
            edge.nodes = copy.(edge.nodes)
        end

        tag!(outline_edges, "_outline")
        elems = [ surfcells; areacells; linecells; outline_edges ]
        nodes = getnodes(elems)

        # observer and light vectors
        # V = Vec3( cosd(elev)*cosd(azim), cosd(elev)*sind(azim), sind(elev) )

        # lightvector===nothing && (lightvector=V) 
        # if lightvector isa AbstractArray
        #     L = lightvector
        # else
        #     error("mplot: lightvector must be a vector.")
        # end
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
    if mesh.env.ndim==3
        project_to_2d!(nodes, mplot.azimut, mplot.elevation, mplot.distance)
        zmin, zmax = extrema(node.coord[3] for node in nodes)

        # raise outline cells
        out_nodes = getnodes(outline_edges)
        for node in out_nodes
            node.coord = node.coord - Vec3(0,0,0.01*(zmax-zmin))
        end

        # distances = [ sum(node.coord[3] for node in elem.nodes)/length(elem.nodes)  for elem in mesh.elems ]
        distances = [ minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        perm = sortperm(distances, rev=true)
        elems = elems[perm]
    end

    # Field 
    has_field = mplot.field != ""

    if has_field
        mplot.label == ""  && (mplot.label = mplot.field)
        
        field = string(mplot.field)
        found = false
        if haskey(elem_data, field)
            fvals = elem_data[field]
            fmax = maximum(fvals[elem.id] for elem in elems)
            fmin = minimum(fvals[elem.id] for elem in elems)
            found = true
        end
        if haskey(node_data, field)
            fvals = node_data[field]
            found = true
            fmax = maximum(fvals[node.id] for node in nodes)
            fmin = minimum(fvals[node.id] for node in nodes)
        end
        found || error("mplot: field $field not found")
        if fmin==fmax
            fmin -= 1
            fmax += 1
        end

        # Colormap
        mplot.values = fvals
        mplot.limits = [fmin, fmax]
        mplot.colormap = resize(mplot.colormap, fmin, fmax, divergefromzero=mplot.args.divergefromzero)

    else
        # Solid colormap
        # mplot.values = 
    end

    mplot.nodes = nodes
    mplot.elems = elems

    # Canvas 
    mplot.canvas = GeometryCanvas()
    mplot.outerpad = 0.01*minimum(mplot.figsize)


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
            else
                bpane = mplot.colorbar.height
            end
        end
    end

    # Canvas box
    canvas = mplot.canvas
    width, height = mplot.figsize
    # pad = mplot.outerpad

    canvas.width = width - rpane - 2*mplot.outerpad
    canvas.height = height - bpane - 2*mplot.outerpad
    canvas.box = [ mplot.outerpad, mplot.outerpad, width-rpane-mplot.outerpad, height-bpane-mplot.outerpad ]

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

    # Draw elements
    for elem in mplot.elems
        # @show elem.shape.name
        if elem.tag=="_outline"
            x, y = elem.nodes[1].coord
            move_to(cc, x, y)
            color = Vec3(0.4, 0.4, 0.4)
            pts = bezier_points(elem)
            curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(cc, 1.4*mplot.lw)
            set_source_rgb(cc, color...) # gray
            stroke(cc)
            continue
        end

        pattern = CairoPatternMesh()
        mesh_pattern_begin_patch(pattern)

        edges = getfaces(elem)
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

        for (i,node) in enumerate(elem.nodes[1:nedges])
            id = node.id
            if has_field && is_nodal_field
                color = mplot.colormap(mplot.values[id])
            end
            mesh_pattern_set_corner_color_rgb(pattern, i-1, color...)
        end

        mesh_pattern_end_patch(pattern)
        set_source(cc, pattern)
        paint(cc)
        
        # draw edges
        gray = sum(mplot.facecolor)/3*Vec3(0.5,0.5,0.5)
        if has_field && !is_nodal_field
            gray = sum(mplot.colormap(mplot.values[elem.id]))/3*Vec3(0.5,0.5,0.5)
        end
        
        x, y = edges[1].nodes[1].coord
        move_to(cc, x, y)
        for edge in edges
            pts = bezier_points(edge)
            x, y = edge.nodes[1].coord
            move_to(cc, x, y)
            curve_to(cc, pts[2]..., pts[3]..., pts[4]...)
            id = edge.nodes[1].id
            if has_field && is_nodal_field
                gray = sum(mplot.colormap(mplot.values[id]))/3*Vec3(0.5,0.5,0.5)
            end
            set_line_width(cc, mplot.lw)
            set_source_rgb(cc, gray...) # gray
            stroke(cc)
        end
    end
    
    # draw colorbar
    has_field && draw!(mplot, cc, mplot.colorbar)
end


function save(mplot::MeshPlot, filename::String)
    width, height = mplot.figsize
    surf = CairoPDFSurface(filename, width, height)
    cc = CairoContext(surf)
    
    configure!(mplot)
    draw!(mplot, cc)
    
    finish(surf)
end