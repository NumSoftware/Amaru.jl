

function Base.delete!(geo::GeoModel, line::Edge)
    filter!(!=(line), geo.edges)

    p1 = line.points[1]
    p2 = line.points[end]

    filter!(!=(line), p1.edges)
    filter!(!=(line), p2.edges)

    # check if line is between two coplanar faces
    if length(line.faces)==2 && line.faces[1].plane==line.faces[2].plane
        s1 = line.faces[1] # to be deleted
        s2 = line.faces[2] # to be extended

        holes = [ s1.loops[2:end]; s2.loops[2:end] ]

        # get inner and outer lines
        innerlines = intersect(s1.loops[1].edges, s2.loops[1].edges)
        outerlines = setdiff(union(s1.loops[1].edges, s2.loops[1].edges), innerlines)

        # find loop and update s2 outer loop
        loops = find_flat_loops(outerlines[1], lines=outerlines; exclude_inner_edges=false)
        loop = addloop!(geo, loops[1])
        s2.loops[1] = loop

        # remove s1 from geo
        filter!(!=(s1), geo.faces)
        filter!(!=(s1.loops[1]), geo.loops)

        # remove references to s1 in s1 lines and update with s2
        for lo in s1.loops
            for l in lo.darts
                l in innerlines && continue
                filter!(!=(s1), l.faces)
                push!(l.faces, s2)
            end
        end

        # remove references to s1 and s2 in inner lines
        for l in innerlines
            filter!(!=(s1), l.faces)
            filter!(!=(s2), l.faces)
        end

        # remove s1 from volumes associated to s2
        for v in s2.volumes
            filter!(!=(s1), v.faces)
        end

        # fix holes
    
    else
        # delete faces associated with that line
        for s in line.faces
            for lo in s.loops
                line in lo.darts && delete!(geo, s)
            end
        end
    end
end


function Base.delete!(geo::GeoModel, face::Face)

    # remove face from geo
    filter!(!=(face), geo.faces)

    if !face.flat
        # delete loop if not flat face
        delete!(geo.loops, face.loops[1]) # todo: check
    else
        # Remove references and check for loop usage
        for loop in face.loops
            # remove face for each edge
            for dart in loop.darts
                filter!(!=(face), dart.edge.faces)
            end 
    
            # delete loop if not used elsewhere
            delete!(geo, loop)
        end
    end

    # remove volumes that uses this face
    if length(face.volumes)==1
        delete!(geo, face.volumes[1])
    elseif length(face.volumes)==2 # join volumes
        v1 = face.volumes[1] # volume to be removed
        v2 = face.volumes[2] # volume to be extended

        innerfaces = intersect(v1.faces, v2.faces)
        outerfaces = setdiff(union(v1.faces, v2.faces), innerfaces)

        # update faces in v2
        v2.faces = outerfaces

        # remove v1 from geo
        filter!(!=(v1), geo.volumes)

        # remove references to v1 in v1 faces and update with v2
        for s in v1.faces
            s in innerfaces && continue
            filter!(!=(v1), s.volumes)
            push!(s.volumes, v2)
        end

        # remove v1 and v2 from innerfaces
        for s in innerfaces
            filter!(!=(v1), s.volumes)
            filter!(!=(v2), s.volumes)
        end
    end

    geo.quiet || println("  Deleted face $(face.id)")

end


function Base.delete!(geo::GeoModel, loop::Loop)
    # deletes a loop if not used by any face

    found = false
    for f in geo.faces
        f.flat || continue
        if loop in f.loops
            found = true
            break
        end
    end

    !found && filter!(!=(loop), geo.loops)
    
end

function Base.delete!(geo::GeoModel, v::Volume)
    filter!(!=(v), geo.volumes)

    # filter v from v faces
    for s in v.faces
        filter!(!=(v), s.volumes)
    end

    geo.quiet || println("  Deleted volume $(v.id)")

end