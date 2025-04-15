
export array!, polar!, revolve!, pull!



function reverse_loop(loop::Loop)
    LoopType = typeof(loop)
    darts = [ flip(dart) for dart in loop.darts ]
    return  LoopType(reverse(darts), id=loop.id)
end





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


function LinearAlgebra.rotate!(path::Path; base=[0.0,0,0], axis=[0.0,0,1], angle=90.0)
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



export picksurface

function picksurface(geo::GeoModel, p::Point)
    return picksurface(geo, p.coord...)
end


function picksurface(geo::GeoModel, x::Real, y::Real, z::Real=0.0)
    p = Point(x, y, z)
    for s in geo.faces
        isin = inside(p, s.loops[1])
        if isin && length(s.loops)>=2
            for lo in s.loops[2:end]
                if inside(p, lo)
                    isin = false
                    break
                end
            end
        end
        isin && return s
    end
    return nothing
end


function tag!(s::Face, tag::String)
    s.tag = tag
end





function join_volumes!(geo::GeoModel, v1::Volume, v2::Volume)

    # select find inner and outer faces
    innerfaces = intersect(v1.faces, v2.faces)
    outerfaces = setdiff(union(v1.faces, v2.faces), innerfaces)

    # remove innerfaces
    for s in innerfaces
        delete!(geo, s)
    end

    # add new volume
    v = addvolume!(geo, collect(outerfaces))

    return v
end


function join_planar_surfaces!(geo::GeoModel, s1::Face, s2::Face)
    @assert s1.flat && s2.flat

    # check if faces are coplanar
    s1.plane != s2.plane && return

    holes = [ s1.loops[2:end]; s2.loops[2:end] ]

    # get inner and outer lines
    innerlines = intersect(s1.loops[1].edges, s2.loops[1].edges)
    outerlines = setdiff(union(s1.loops[1].edges, s2.loops[1].edges), innerlines)

    # remove inner lines
    for l in innerlines
        delete!(geo, l)
    end

    # find loop
    loops = find_flat_loops(outerlines[1], lines=outerlines; exclude_inner_edges=false)
    loop = addloop!(geo, loops[1])
    s = addface!(geo, loops[1])
    # geo.quiet || println("  Added face $s.id")

    return s
end




function split_flat_face!(geo::GeoModel, face::Face, loop1::Loop, loop2::Loop)
    @assert face.flat && loop1.flat && loop2.flat

    # add new faces
    newfaces = Face[]
    for loop in ( loop1, loop2 )
        loop = addloop!(geo, loop)
        f = Face(loop, tag=face.tag)
        
        # add face
        geo._id +=1
        f.id = geo._id
        push!(geo.faces, f)
        push!(newfaces, f)
    end

    face1, face2 = newfaces

    # remove face from geo
    filter!(!=(face), geo.faces)
    geo.quiet || println("  Deleted $(face.flat ? "flat" : "curved") face $(face.id)")

    # remove loop from geo if not used in other faces
    delete!(geo, face.loops[1])

    # remove references to face in face edges, including holes
    for lo in face.loops
        for dart in lo.darts
            filter!(!=(face), dart.edge.faces)
        end
    end
    
    # update holes and references in hole edges
    for f in [face1, face2]
        for loop in face.loops[2:end]
            if inside(loop, f.loops[1])
                push!(f.loops, loop)
                for dart in loop.darts
                    push!(dart.edge.faces, f)
                end
            end
        end
    end

    # update references to f in edges
    for f in (face1, face2)
        for loop in f.loops
            for dart in loop.darts
                push!(dart.edge.faces, f)
            end
        end
    end

    # update volumes
    if length(face.volumes)>0
        for volume in face.volumes
            # update spins
            idx = findfirst(s->s.face==face, volume.spins)
            face_spin = volume.spins[idx]
            filter!(!=(face_spin), volume.spins) # remove spin from volume
            normal = face_spin.normal
            spin1 = FaceSpin(face1, normal)
            spin2 = FaceSpin(face2, normal)
            push!(volume.spins, spin1)
            push!(volume.spins, spin2)

            # update face1 and face2 volumes
            push!(face1.volumes, volume)
            push!(face2.volumes, volume)
        end
    end

    if !geo.quiet 
        for face in newfaces
            println("  Added $(face.flat ? "flat" : "curved") face $(face.id)")
        end
    end

    return newfaces
end



# function picksurface(geo::GeoModel, p::Point)
#     return picksurface(geo, p.coord...)
# end


# function picksurface(geo::GeoModel, x::Real, y::Real, z::Real=0.0)
#     p = Point(x, y, z)
#     for s in geo.faces
#         isin = inside(p, s.loops[1])
#         if isin && length(s.loops)>=2
#             for lo in s.loops[2:end]
#                 if inside(p, lo)
#                     isin = false
#                     break
#                 end
#             end
#         end
#         isin && return s
#     end
#     return nothing
# end


# function tag!(s::Face, tag::String)
#     s.tag = tag
# end


# function is_neighbor(surf1::Face, surf2::Face)
#     for dart in surf1.loops[1].darts
#         for s in dart.edge.faces
#             s==surf2 && return true
#         end
#     end
#     return false
# end
