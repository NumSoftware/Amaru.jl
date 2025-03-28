
export array!, polar!, revolve!, pull!

function move!(subpath::SubPath; dx::Real=0.0, dy::Real=0.0, dz::Real=0.0)
    for p in subpath.path.points
        p.coord = p.coord + Vec3(dx, dy, dz)
    end
end



function array!(geo::GeoModel, subpath::SubPath; nx=1, ny=1, nz=1, dx=0.0, dy=0.0, dz=0.0)
    for k in 0:nz-1
        for j in 0:ny-1
            for i in 0:nx-1
                i==j==k==0 && continue
                cp = copy(subpath)
                move!(cp, dx=i*dx, dy=j*dy, dz=k*dz)
                addsubpath!(geo, cp)
            end
        end
    end
end


function revolve!(
    geo::GeoModel, 
    line::Line; 
    base::AbstractArray{<:Real,1} = Float64[],
    axis::AbstractArray{<:Real,1} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -180 < angle < 180

    axis = Vec3(normalize(axis))
    base = Vec3(base)

    θ = angle*pi/180
    R = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
    
    C1 = line.points[1].coord
    C2 = line.points[2].coord
    C3 = base + R*(C1 - base)*conj(R)
    C4 = base + R*(C2 - base)*conj(R)

    sz1 = line.points[1].size
    sz2 = line.points[2].size

    p1 = Point(C1, size=sz1)
    p2 = Point(C2, size=sz2)
    p3 = Point(C3, size=sz1)
    p4 = Point(C4, size=sz2)

    if p1==p3 && p2==p4 # line is on the axis
        return nothing
    end
    
    C5 = base + ( dot(C1-base, axis) )*axis
    C6 = base + ( dot(C2-base, axis) )*axis
    
    p5 = addsinglepoint!(geo, Point(C5)) 
    p6 = addsinglepoint!(geo, Point(C6))

    l2 = addline!(geo, p3, p4)
    # l2 = addsingleline!(geo, p3, p4)

    if p1==p3 # p1 is on the axis
        a2 = addarc!(geo, p2, p6, p4)
        if p5==p6 # same center
            loop = PlaneLoop(line, a2, l2)
            s = addplanesurface!(geo, loop)
        else
            loop = Loop(line, a2, l2)
            s = addsurface!(geo, loop)
        end    
    elseif p2==p4 # p2 is on the axis
        a1 = addarc!(geo, p1, p5, p3)
        if p5==p6 # same center
            loop = PlaneLoop(line, a1, l2)
            s = addplanesurface!(geo, loop)
        else
            loop = Loop(line, a1, l2)
            s = addsurface!(geo, loop)
        end
    else # p1 and p2 are not on the axis
        a1 = addarc!(geo, p1, p5, p3)
        a2 = addarc!(geo, p2, p6, p4)
        if p5==p6 # same center
            loop = PlaneLoop(line, a1, l2, a2)
            s = addplanesurface!(geo, loop)
        else
            loop = Loop(line, a1, l2, a2)
            s = addsurface!(geo, loop)
        end
    end

    return s
end


function revolve!(
    geo::GeoModel, 
    surf::PlaneFace; 
    base::AbstractArray{<:Real,1} = Float64[],
    axis::AbstractArray{<:Real,1} = Float64[],
    angle   ::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -180 < angle < 180

    # revolve all lines
    for lo in surf.loops
        for line in lo.lines
            revolve!(geo, line, base=base, axis=axis, angle=angle)
        end
    end

    # add lid
    
end


function revolve!(
    geo::GeoModel, 
    surfs::Vector{<:AbstractFace};
    base::AbstractArray{<:Real,1} = Float64[],
    axis::AbstractArray{<:Real,1} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -180 < angle < 180

    # revolve all lines
    surfs = copy(surfs) # detach from geo.surfaces since the latter gets updated
    for s in surfs
        revolve!(geo, s, base=base, axis=axis, angle=angle)
    end
end


function LinearAlgebra.rotate!(path::Path; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )
    axis = normalize(Vec3(axis))
    base = Vec3(base)
    θ    = angle*pi/180
    R    = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
    digs = 8

    local X
    for node in path.points
        X          = base + R*(node.coord-base)*conj(R)
        node.coord = round.(X, digits=digs)
    end
end


function LinearAlgebra.rotate!(subpath::SubPath; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0 )
    rotate!(subpath.path, base=base, axis=axis, angle=angle)
end


function polar!(geo::GeoModel, subpath::SubPath; base=[0.0,0,0], axis=[0.0,0,1], angle=360.0, n=2)
    Δθ = angle/n
    for i in 1:n-1
        cp = copy(subpath)
        rotate!(cp, base=base, axis=axis, angle=Δθ*i)
        addsubpath!(geo, cp)
    end
end



function pull!(geo::GeoModel, line::Line; axis=[0.,0,1], length=1.0)
    p1, p2 = line.points
    dx = length*axis[1]
    dy = length*axis[2]
    dz = length*axis[3]
    p3 = copy!(geo, p2, dx=dx, dy=dy, dz=dz)
    p4 = copy!(geo, p1, dx=dx, dy=dy, dz=dz)
    l1 = addline!(geo, p1, p4, tag=line.tag)
    l2 = addline!(geo, p2, p3, tag=line.tag)
    l3 = addline!(geo, p3, p4, tag=line.tag)

    # lo = PlaneLoop(l1, l2, l3, l4)
    # s = addplanesurface!(geo, lo)

    # @assert s!==nothing
    # return s
end


# function pull_old!(geo::GeoModel, line::Line; axis=[0.,0,1], length=1.0)
#     p1, p2 = line.points
#     dx = length*axis[1]
#     dy = length*axis[2]
#     dz = length*axis[3]
#     p3 = copy!(geo, p2, dx=dx, dy=dy, dz=dz)
#     p4 = copy!(geo, p1, dx=dx, dy=dy, dz=dz)
#     l1 = addsingleline!(geo, p1, p2, tag=line.tag)
#     l2 = addsingleline!(geo, p2, p3, tag=line.tag)
#     l3 = addsingleline!(geo, p3, p4, tag=line.tag)
#     l4 = addsingleline!(geo, p4, p1, tag=line.tag)

#     lo = PlaneLoop(l1, l2, l3, l4)
#     s = addplanesurface!(geo, lo)

#     @assert s!==nothing
#     return s
# end


function pull!(geo::GeoModel, arc::Arc; axis=[0,0,1], length=1.0)
    p1, p2, p3 = arc.points
    dx = length*axis[1]
    dy = length*axis[2]
    dz = length*axis[3]
    p4 = copy!(geo, p3, dx=dx, dy=dy, dz=dz)
    p5 = copy!(geo, p2, dx=dx, dy=dy, dz=dz)
    p6 = copy!(geo, p1, dx=dx, dy=dy, dz=dz)

    c1 = addsinglearc!(geo, p1, p2, p3, tag=arc.tag)
    c2 = addsingleline!(geo, p3, p4, tag=arc.tag)
    c3 = addsinglearc!(geo, p4, p5, p6, tag=arc.tag)
    c4 = addsingleline!(geo, p6, p1, tag=arc.tag)

    lo = Loop(c1, c2, c3, c4)

    geo._id +=1
    lo.id = geo._id
    push!(geo.loops, lo)

    s = addsurface!(geo, lo)
    return s
end



function pull!(geo::GeoModel, surf::PlaneFace; axis=[0.,0,1], length=1.0)

    length==0 && return

    ntotalvols = Base.length(geo.volumes)
    nsurfvols = Base.length(surf.volumes)

    # pull lateral lines
    for lo in surf.loops
        for line in lo.lines
            s = pull!(geo, line, axis=axis, length=length)
        end
    end

    nnewvols = Base.length(geo.volumes)-ntotalvols
    if nnewvols>1
        error("more than one volume found in pull operation")
    elseif nnewvols==0
        error("no volume found in pull operation")
    end

    # join volumes if surf was part of a volume
    if nsurfvols>0
        # if Base.length(geo.volumes)-ntotalvols>1
            # error("more than one volume found in pull operation")
        # end
        
        delete!(geo, surf)

        for l in getlines(surf)
            # check for lines in coplanar faces
            if Base.length(l.surfaces)==2 && l.surfaces[1].plane==l.surfaces[2].plane
                delete!(geo, l)
            end
            if Base.length(l.surfaces)==1 # for push
                delete!(geo, l)
            end
        end
        
    end

end

# function pull_old!(geo::GeoModel, surf::PlaneFace; axis=[0.,0,1], length=1.0)

#     # check if surf is an inner surface
#     # innersurf = nothing
#     # for s in geo.surfaces
#     #     for lo in s.loops[2:end]
#     #         if lo==surf.loop
#     #             innersurf=s
#     #             break
#     #         end
#     #     end
#     # end

#     surfs = AbstractFace[]

#     # extrude lateral lines
#     for lo in surf.loops
#         for line in lo.lines
#             s = pull!(geo, line, axis=axis, length=length)
#             s.tag = surf.tag
#             push!(surfs, s)
#         end
#     end

#     # find lid loops (outer and inner loops if existent)
#     loops = PlaneLoop[]
#     for lo in surf.loops
#         lines = AbstractLine[]
#         for line in lo.lines

#             points = Point[]
#             for p in line.points
#                 pp = Point(p.coord .+ length.*axis)
#                 pp = getpoint(geo, pp) # point should exists
#                 push!(points, pp)
#             end

#             if Base.length(points)==2
#                 l = getline(geo, Line(points...))
#             else
#                 l = getline(geo, Arc(points...))
#             end
#             push!(lines, l)
#         end
#         lo = PlaneLoop(lines...)
#         push!(loops, lo)
#     end

#     # add single loops
#     for lo in loops
#         geo._id +=1
#         lo.id = geo._id
#         push!(geo.loops, lo)
#     end

#     # add closing lid
#     for lo in loops
#         s = addplanesurface!(geo, lo, tag=surf.tag)
#         push!(surfs, s)
#     end

#     # check if surf is part of a volume
#     volume = nothing
#     for v in geo.volumes
#         if surf in v.surfaces
#             volume = v
#             break
#         end
#     end

#     if volume===nothing
#         push!(surfs, surf)
#         volume = addvolume!(geo, surfs, tag=surf.tag)
#     else
#         # remove surf from geo
#         filter!(!=(surf), geo.surfaces)
        
#         # remove surf from volume
#         filter!(!=(surf), volume.surfaces)
        
#         # remove unused surface loops from geo
#         for lo in surf.loops
#             # for each curve remove links to this surface
#             for l in lo.lines
#                 filter!(!=(surf), l.surfaces)
#             end 
#         end

#         # add new surfaces to volume
#         for s in surfs
#             push!(volume.surfaces, s)
#         end
#     end

#     for s in surfs
#         push!(s.volumes, volume)
#     end

#     return volume
# end


function pull!(geo::GeoModel, surfs::Vector{<:AbstractFace}; axis=[0.,0,1], length=1.0)
    surfs = copy(surfs) # make a copy
    for s in surfs
        pull!(geo, s; axis=axis, length=length)
    end
end


function pull!(m::GeoModel; nargs...)
    pull!(m, m.surfaces; nargs...)
end
