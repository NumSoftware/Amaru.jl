# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


GeoPlot_params = [
    FunInfo(:GeometryPlot, "Creates a customizable `GeometryPlot` instance used to plot a model geometry."),
    ArgInfo(:geo, "2D/3D model geometry"),
    KwArgInfo((:size, :figsize), "Chart size in dpi", (220,150), length=2),
    KwArgInfo(:facecolor, "Face color", :aliceblue),
    KwArgInfo((:lw, :lineweight), "Edge weight", 0.4,  cond=:(lw>0) ),
    KwArgInfo(:font, "Font name", "NewComputerModern", type=AbstractString),
    KwArgInfo((:fontsize, :colorbarfontsize, :cbfontsize), "Colorbar font size", 7.0, cond=:(fontsize>0)),
    KwArgInfo(:azimut, "Azimut angle for 3d in degrees", 30 ),
    KwArgInfo(:elevation, "Elevation angle for 3d in degrees", 30 ),
    KwArgInfo(:distance, "Distance from camera in 3d", 0.0, cond=:(distance>=0) ),
    KwArgInfo((:lightvector, :lv), "Light direction vector", [0.0,0.0,0.0], length=3 ),
    KwArgInfo((:pointlabels, :nodelabels), "Flag to show point labels", true, type=Bool ),
    KwArgInfo(:linelabels, "Flag to show line labels", true, type=Bool ),
    KwArgInfo(:surfacelabels, "Flag to show surface labels", true, type=Bool ),
    KwArgInfo(:volumelabels, "Flag to show volume labels", true, type=Bool ),
]
@doc docstring(GeoPlot_params) GeometryPlot(geo; kwargs...)


mutable struct GeometryPlot<:Figure
    geo::GeoModel
    points::Vector{Point}
    faces::Vector{Face}
    pointlabels::Bool
    linelabels::Bool
    surfacelabels::Bool
    volumelabels::Bool
    canvas::Union{FigureComponent, Nothing}
    outerpad::Float64
    shades::Vector{Float64}
    args::NamedTuple
    font::String
    fontsize::Float64
    azimut::Float64
    elevation::Float64
    distance::Float64
    
    lightvector::Vector{Float64}

    width::Float64
    height::Float64

    lw::Float64
    facecolor::Tuple

    function GeometryPlot(geo; args...)
        args = checkargs(args, GeoPlot_params, aliens=false)

        this = new()
        this.geo = copy(geo)
        # this.geo = geo

        this.canvas = nothing
        this.pointlabels = args.pointlabels
        this.linelabels = args.linelabels
        this.surfacelabels = args.surfacelabels
        this.volumelabels = args.volumelabels
        this.outerpad = 0.0

        this.width     = args.size[1]
        this.height    = args.size[2]
        this.lw        = args.lw
        this.facecolor = _colors_dict[args.facecolor]

        this.font        = args.font
        this.fontsize    = args.fontsize
        this.azimut      = args.azimut
        this.elevation   = args.elevation
        this.distance    = args.distance
        this.lightvector = args.lightvector

        this.args = args

        return this
    end
end


# function bezier_points(edge)
#     p1 = edge.points[1].coord[1:2]
#     p4 = edge.points[2].coord[1:2]
#     ξ1 = -1/3
#     ξ2 = +1/3
#     C = getcoords(edge.points, 2)
#     p2 = C'*edge.shape.func([ξ1])
#     p3 = C'*edge.shape.func([ξ2])

#     cp2 = 1/6*(-5*p1+18*p2- 9*p3+2*p4)
#     cp3 = 1/6*( 2*p1- 9*p2+18*p3-5*p4)

#     return p1, cp2, cp3, p4
# end


function project_to_2d!(points::Vector{Point}, azimut, elevation, distance)
    # Find bounding box
    xmin, xmax = extrema( point.coord[1] for point in points)
    ymin, ymax = extrema( point.coord[2] for point in points)
    zmin, zmax = extrema( point.coord[3] for point in points)
    reflength = max(xmax-xmin, ymax-ymin, zmax-zmin)

    # Centralize 
    center = 0.5*Vec3(xmin+xmax, ymin+ymax, zmin+zmax)
    for point in points
        point.coord = point.coord - center
    end

    # Rotation around z axis
    θ = -azimut*pi/180
    R = Quaternion(cos(θ/2), 0, 0, sin(θ/2))
    for point in points
        point.coord = (R*point.coord*conj(R))[2:4]
    end

    # Rotation around y axis
    θ = elevation*pi/180
    R = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
    for point in points
        point.coord = (R*point.coord*conj(R))[2:4]
    end

    # Set projection values
    distance==0 && (distance=reflength*3)
    distance = max(distance, reflength)
    focal_length = 0.1*distance

    # Make projection
    for point in points
        x = point.coord[1]
        y′ = point.coord[2]*focal_length/(distance-x)
        z′ = point.coord[3]*focal_length/(distance-x)
        point.coord =  Vec3(y′, z′, distance-x)
    end

    xmin, xmax = extrema( point.coord[1] for point in points)
    ymin, ymax = extrema( point.coord[2] for point in points)

    # normalize
    l = max(xmax-xmin, ymax-ymin)
    for point in points
        x = (point.coord[1]-xmin)/l
        y = (point.coord[2]-ymin)/l
        point.coord =  Vec3(x, y, point.coord[3])
    end
end


function configure!(gplot::GeometryPlot)
    geo = gplot.geo
    ndim = geo.ndim
    points = geo.points
    edges = geo.edges
    faces = geo.faces

    # add extra points from arcs
    for line in edges
        if line isa Arc
            for point in line.extrapoints
                push!(points, point)
            end
        end
    end

    # 3D -> 2D projection
    if geo.ndim==3

        # compute projection
        project_to_2d!(points, gplot.azimut, gplot.elevation, gplot.distance)
        zmin, zmax = extrema(point.coord[3] for point in points)

        # reorder faces
        distances = Float64[]
        for surf in faces
            spoints = getpoints(surf)
            d = 0.9*sum(point.coord[3] for point in spoints)/length(spoints) + 0.1*minimum(point.coord[3] for point in spoints) 
            push!(distances, d)
        end
        perm = sortperm(distances, rev=true)
        faces = faces[perm]
    end

    gplot.points = points
    gplot.faces = faces

    # Canvas 
    gplot.canvas = Canvas()
    width, height = gplot.width, gplot.height
    gplot.outerpad = 0.025*min(width, height)

    rpane = 0.0
    bpane = 0.0

    # Canvas box
    canvas = gplot.canvas

    canvas.width = width - rpane - 2*gplot.outerpad
    canvas.height = height - bpane - 2*gplot.outerpad
    Xmin = gplot.outerpad
    Ymin = gplot.outerpad
    Xmax = width - rpane - gplot.outerpad
    Ymax = height - bpane - gplot.outerpad
    canvas.box = [ gplot.outerpad, gplot.outerpad, width-rpane-gplot.outerpad, height-bpane-gplot.outerpad ]

    # Canvas limits
    xmin, xmax = extrema( point.coord[1] for point in gplot.points)
    ymin, ymax = extrema( point.coord[2] for point in gplot.points)

    ratio = min((Xmax-Xmin)/(xmax-xmin), (Ymax-Ymin)/(ymax-ymin) )
    dX = 0.5*((Xmax-Xmin)-ratio*(xmax-xmin))
    dY = 0.5*((Ymax-Ymin)-ratio*(ymax-ymin))

    xmin = xmin - dX/ratio
    xmax = xmax + dX/ratio
    ymin = ymin - dY/ratio
    ymax = ymax + dY/ratio
    gplot.canvas.limits = [xmin, ymin, xmax, ymax]

end


function bezier_points_arc(P1, P2, P3, P4)
    p1 = P1.coord[1:2]
    p2 = P2.coord[1:2]
    p3 = P3.coord[1:2]
    p4 = P4.coord[1:2]
    
    cp2 = 1/6*(-5*p1+18*p2- 9*p3+2*p4)
    cp3 = 1/6*( 2*p1- 9*p2+18*p3-5*p4)

    return p1, cp2, cp3, p4
end


function draw!(gplot::GeometryPlot, ctx::CairoContext)
    set_line_join(ctx, Cairo.CAIRO_LINE_JOIN_ROUND)
 
    Xmin, Ymin, Xmax, Ymax = gplot.canvas.box
    xmin, ymin, xmax, ymax = gplot.canvas.limits

    ratio = (Xmax-Xmin)/(xmax-xmin) # == (Ymax-Ymin)/(ymax-ymin)
    set_matrix(ctx, CairoMatrix([ratio, 0, 0, -ratio, Xmin-xmin*ratio, Ymax+ymin*ratio]...))

    font = get_font(gplot.font)
    select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(ctx, gplot.fontsize)
 
    function getpoints(loop::Loop)
        points = [ loop.darts[1].points[1] ]
        for dir_edge in loop.darts
            # @show dir_edge.forward
            # @show dir_edge.points[1].id
            # @show dir_edge.points[2].id
            if dir_edge isa Arc
                if dir_edge.forward
                    push!(points, dir_edge.edge.extrapoints...)
                else
                    push!(points, reverse(dir_edge.edge.extrapoints)...)
                end
            end
            push!(points, dir_edge.points[end])
        end
        return points
    end

    # function getpoints(surf::Face)
    #     points = Set{Point}()
    #     for dir_edge in surf.loops[1].darts
    #         push!(points, dir_edge.points[1])
    #         push!(points, dir_edge.points[end])
    #     end
    #     return collect(points)
    # end

    function getpoints(vol::Volume)
        points = Set{Point}()
        for s in vol.spins
            push!(points, getpoints(s.face.loops[1])...)
        end
        return collect(points)
    end


    # Draw orphan edges
    for edge in gplot.geo.edges
        length(edge.faces)==0 || continue

        x, y = edge.points[1].coord
        move_to(ctx, x, y)

        if edge isa Line
            x, y = edge.points[end].coord
            line_to(ctx, x, y)
        elseif edge isa Arc
            p1 = edge.points[1]
            p4 = edge.points[3]
            p2 = edge.extrapoints[1]
            p3 = edge.extrapoints[3]
    
            pts = bezier_points_arc(p1, p2, p3, p4)
    
            curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
        end
        edgecolor = gplot.facecolor.*0.45
        set_source_rgb(ctx, edgecolor...) # edges
        set_line_width(ctx, gplot.lw)
        stroke(ctx)

        
        # draw labels
        Cairo.save(ctx)
        set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
        
        # point labels
        for point in edge.points
            x, y = data2user(gplot.canvas, point.coord[1], point.coord[2])
            set_source_rgb(ctx, _colors_dict[:blue]...)
            gplot.pointlabels && draw_text(ctx, x, y, string(point.id), halign="center", valign="center", angle=0)
        end

        # edge labels
        if edge isa Line
            x1, y1 = data2user(gplot.canvas, edge.points[1].coord[1], edge.points[1].coord[2])
            x2, y2 = data2user(gplot.canvas, edge.points[end].coord[1], edge.points[end].coord[2])
            x = 0.5*(x1+x2)
            y = 0.5*(y1+y2)
            set_source_rgb(ctx, _colors_dict[:green]...)
            draw_text(ctx, x, y, string(edge.id), halign="center", valign="center", angle=0)
        else
            x, y = data2user(gplot.canvas, edge.extrapoints[2].coord[1], edge.extrapoints[2].coord[2])
            set_source_rgb(ctx, _colors_dict[:green]...)
            draw_text(ctx, x, y, string(edge.id), halign="center", valign="center", angle=0)
        end

        Cairo.restore(ctx)

    end


    # Draw elements
    for face in gplot.faces
        draw_face!(ctx, gplot, face)
        
        Cairo.save(ctx)

        # Draw surface label
        set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
        points = getpoints(face.loops[1])
        # points = [ face.loops[1].darts[1].points[1] ]
        # for dir_edge in face.loops[1].darts
        #     # push!(points, dir_edge.points[1])
        #     push!(points, dir_edge.points[end])
        #     if dir_edge isa Arc
        #         push!(points, dir_edge.edge.extrapoints...)
        #     end
        # end
        C = sum(point.coord for point in points)/length(points)
        x, y = data2user(gplot.canvas, C[1], C[2])
        set_source_rgb(ctx, _colors_dict[:red]...)
        draw_text(ctx, x, y, string(face.id), halign="center", valign="center", angle=0)

        # Draw points labels
        for loop in face.loops
            points = getpoints(loop)
            for point in points
                x, y = data2user(gplot.canvas, point.coord[1], point.coord[2])
                set_source_rgb(ctx, _colors_dict[:blue]...)
                gplot.pointlabels && draw_text(ctx, x, y, string(point.id), halign="center", valign="center", angle=0)
            end
        end

        # Draw edges label
        for loop in face.loops
            for dir_edge in loop.darts
                if dir_edge.edge isa Line
                    x1, y1 = data2user(gplot.canvas, dir_edge.points[1].coord[1], dir_edge.points[1].coord[2])
                    x2, y2 = data2user(gplot.canvas, dir_edge.points[end].coord[1], dir_edge.points[end].coord[2])
                    x = 0.5*(x1+x2)
                    y = 0.5*(y1+y2)
                    set_source_rgb(ctx, _colors_dict[:green]...)
                    draw_text(ctx, x, y, string(dir_edge.edge.id), halign="center", valign="center", angle=0)
                else
                    x, y = data2user(gplot.canvas, dir_edge.edge.extrapoints[2].coord[1], dir_edge.edge.extrapoints[2].coord[2])
                    set_source_rgb(ctx, _colors_dict[:green]...)
                    draw_text(ctx, x, y, string(dir_edge.edge.id), halign="center", valign="center", angle=0)
                end
            end
        end

        Cairo.restore(ctx)
    end

    # Draw volumes label
    Cairo.save(ctx)
    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    for vol in gplot.geo.volumes
        points = getpoints(vol)
        C = sum(point.coord for point in points)/length(points)
        x, y = data2user(gplot.canvas, C[1], C[2])
        set_source_rgb(ctx, _colors_dict[:orange]...)
        draw_text(ctx, x, y, string(vol.id), halign="center", valign="center", angle=0)
    end
    Cairo.restore(ctx)

end


function draw_face!(ctx::CairoContext, gplot::GeometryPlot, face::Face)

    ndim = gplot.geo.ndim
    if length(face.volumes)==0
        color = (0.911, 0.973, 1.0)
    else
        color = (0.9, 0.9, 1)
    end

    on_holes = false
    loop = face.loops[1]

    shade = 1.0

    on_main_loop = true

    for loop in face.loops
        darts = loop.darts

        x, y = darts[1].points[1].coord
        move_to(ctx, x, y)
        for dir_edge in darts
            if dir_edge.edge isa Line
                x, y = dir_edge.points[end].coord
                line_to(ctx, x, y)
                continue
            end
            if dir_edge.edge isa Arc
                p1 = dir_edge.points[1]
                p4 = dir_edge.points[3]
                if dir_edge.forward
                    p2 = dir_edge.edge.extrapoints[1]
                    p3 = dir_edge.edge.extrapoints[3]
                else
                    p2 = dir_edge.edge.extrapoints[3]
                    p3 = dir_edge.edge.extrapoints[1]
                end

                pts = bezier_points_arc(p1, p2, p3, p4)

                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
                continue
            end
        end
        
    end
    
    Cairo.set_source_rgba(ctx, color..., 0.5) 
    set_fill_type(ctx, Cairo.CAIRO_FILL_RULE_EVEN_ODD)
    fill_preserve(ctx)

    edgecolor = gplot.facecolor.*0.45
    set_source_rgb(ctx, edgecolor...) # edges
    set_line_width(ctx, gplot.lw)
    stroke(ctx)
end

