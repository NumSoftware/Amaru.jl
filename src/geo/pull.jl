

function pull!(geo::GeoModel, line::Line; axis=[0.,0,1], length=1.0)

    p1, p2 = line.points
    dx, dy, dz = length*axis
    p3 = copy!(geo, p2, dx=dx, dy=dy, dz=dz)
    p4 = copy!(geo, p1, dx=dx, dy=dy, dz=dz)
    
    addline!(geo, p1, p4, tag=line.tag)
    addline!(geo, p2, p3, tag=line.tag)
    addline!(geo, p3, p4, tag=line.tag)

end



function pull!(geo::GeoModel, arc::Arc; axis=[0.,0,1], length=1.0)

    p1, p2, p3 = arc.points
    dx, dy, dz = length*axis

    p4 = copy!(geo, p3, dx=dx, dy=dy, dz=dz)
    p5 = copy!(geo, p2, dx=dx, dy=dy, dz=dz)
    p6 = copy!(geo, p1, dx=dx, dy=dy, dz=dz)

    c1 = addsinglearc!(geo, p1, p2, p3, tag=arc.tag)
    c2 = addsingleline!(geo, p3, p4, tag=arc.tag)
    c3 = addsinglearc!(geo, p4, p5, p6, tag=arc.tag)
    c4 = addsingleline!(geo, p1, p6, tag=arc.tag)

    loop = Loop([c1, c2, c3, c4], flat=false)

    geo._id +=1
    loop.id = geo._id
    push!(geo.loops, loop)
    addface!(geo, loop)

end


function pull!(geo::GeoModel, face::Face; axis=[0.,0,1], length=1.0)
    # disable automatic volume detection
    geo._volume_detection = false 

    face.flat || error("pull! only works for flat faces")

    length==0.0 && return

    # pull lateral lines
    for loop in face.loops
        for dart in loop.darts
            pull!(geo, dart.edge, axis=axis, length=length)
        end
    end

    # get a seed point to find the lid loop
    p = Point(face.loops[1].points[1].coord + length.*axis)
    p = getpoint(geo, p) # point should exists

    # find normal
    normal = getnormal(face.loops[1])
    dot(normal, axis) > 0.0 && (normal = -normal)

    loops = find_flat_loops(p, Edge[], normal, exclude_inner_edges=true)
    loop = addloop!(geo, loops[1])
    lid_face = getface(geo, loop)

    # add lid face
    if lid_face===nothing
        addface!(geo, loop)
    end

    # remove extra faces if face has holes
    for loop in face.loops[2:end]
        # find a face with the same loop
        found = false
        for f in geo.faces
            if loop==f.loops[1]
                found = true
                break
            end
        end

        found || continue

        # get hash of moved points
        _points = [ Point(p.coord + length.*axis) for p in loop.points ]
        _hs = sum( hash(p) for p in _points )
        
        # find a face with a loop with same hash
        for f in geo.faces
            hs = sum( hash(p) for p in f.loops[1].points )
            if hs==_hs
                delete!(geo, f)
                break
            end
        end
    end


    # make a new volume

    # find_face_loops(FaceSpin(face, normal))

    loop = find_face_loops(FaceSpin(face, normal))[1]
    addvolume!(geo, loop.spins)

    # restore automatic volume detection
    geo._volume_detection = false 

end


function pull!(geo::GeoModel, surfs::Vector{<:Face}; axis=[0.,0,1], length=1.0)
    surfs = copy(surfs) # make a copy
    for s in surfs
        pull!(geo, s; axis=axis, length=length)
    end
end


function pull!(m::GeoModel; nargs...)
    pull!(m, m.surfaces; nargs...)
end
