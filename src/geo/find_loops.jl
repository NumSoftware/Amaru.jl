
function get_angle(dart1::Dart, dart2::Dart, normal::AbstractArray{Float64})
    V1 = dir(dart1)
    V2 = dir(dart2)

    θ = acos(dot(V1,V2))
    if dot(cross(V1,V2), normal) > 0
        return pi - θ
    else
        return pi + θ
    end
    
end


# Gets the numnber of flat faces adjacent to a dart perpendicular to normal
function get_n_flat_faces(dart::Dart, normal::Maybe{AbstractArray{Float64}})
    tol = 1e-6
    n = 0
    for face in dart.edge.faces
        face.flat || continue
        if normal!==nothing 
            # check if face normal is parallel to normal
            face_normal = getnormal(face.loops[1])
            norm(cross(face_normal, normal)) < tol || continue
            # norm(cross(face.plane.normal, normal)) < tol || continue
        end
        n = n+1
    end
    return n
end


function get_next_darts(dart::Dart, subset::Vector{<:Edge}, normal::Maybe{AbstractArray{Float64}}, has_side::Bool, exclude_inner_edges::Bool)
    tol = 1e-6
    cands = Dart[]
    last_point = dart.points[end] 
    for edge in last_point.edges
        edge == dart.edge && continue
        
        # check if edge is in the subset
        length(subset) > 0 && (edge in subset || continue)

        forward = last_point==edge.points[1]
        new_dart = Dart(edge, forward=forward)

        # check if new_dart is not a dead end
        length(new_dart.points[end].edges) > 1 || continue

        if normal!==nothing 
            # check if edge is coplanar
            coplanar(last_point, normal, edge.points) || continue

            # check if new_dart is perpendicular to the normal
            # abs(dot(dir(new_dart), normal)) < tol || continue

            # check if new_dart has less than 2 flat faces
            exclude_inner_edges && ( get_n_flat_faces(new_dart, normal) < 2 || continue )
        end
        push!(cands, new_dart)
    end


    # filter according to angle and normal
    if normal!==nothing && length(cands)>1
        # get angles for dart with each candidate
        angles = Float64[]
        for cand in cands
            angle = get_angle(dart, cand, normal)
            push!(angles, angle)
        end

        # sort
        idxs = sortperm(angles)
        cands = cands[idxs]

        if has_side
            cands = [ cands[1] ] # take the first candidate
        else
            cands = [ cands[1], cands[end] ] # take the first and last candidates
        end
    end

    return cands
end



# Algorithm for finding loops on flat faces
# 1. Start with a edge -> dart. If the optional list of edges is provided, the loops will only be searched inside those edges
# 2. Starting from the last point of visited darts, generate a list of darts with a point equal to the last point of the dart
# 2. Set a list of visited darts (chord)
# 3. Check if the new dart is coplanar with the previous darts
# 4. If the new dart is coplanar, add it to the chord and go to step 2
# 6. If the chord is closed, add it to the list of loops
# 7. If the first dart repeats with the in the chord, discard the chord
# 8. Continue the process after updating the list of visited darts

function find_flat_loops(point::Point, subset::Vector{<:Edge}=Edge[], normal::Maybe{AbstractArray{Float64}}=nothing; exclude_inner_edges::Bool=false)
    tol = 1e-6
    darts = Dart[]

    for edge in point.edges
        
        # check if edge is in the subset
        length(subset) > 0 && (edge in subset || continue)

        forward = point==edge.points[1]
        dart = Dart(edge, forward=forward)

        # check if dart is not a dead end
        length(dart.points[end].edges) > 1 || continue

        # check if dart is perpendicular to the normal
        abs(dot(dir(dart), normal)) < tol || continue

        # check if dart has less than 2 flat faces
        exclude_inner_edges && ( get_n_flat_faces(dart, normal) < 2 || continue )

        # check if dart is coplanar
        normal==nothing || (coplanar(point, normal, edge.points) || continue)

        push!(darts, dart)
    end

    return find_flat_loops(darts[1], subset, normal; exclude_inner_edges=exclude_inner_edges)

end



function find_flat_loops(dart::Dart, subset::Vector{<:Edge}=Edge[], normal::Maybe{AbstractArray{Float64}}=nothing; exclude_inner_edges::Bool=false)

    function find_loops(visited::Vector{Dart}, dart::Dart, normal::Maybe{AbstractArray{Float64}}, has_side::Bool, exclude_inner_edges::Bool)

        # length(visited)>0 &&coplanar(visited, dart)
        # normal==nothing || (coplanar(visited[end].points[end], normal, dart.edge.points) || return Loop[])

        if length(visited)>=2
            dart==visited[1] && return [ Loop(visited, flat=true) ]
            dart in visited && return Loop[]

            dart.points[end] == visited[1].points[end] && return Loop[] # bad closing
        end


        visited = [visited; dart] # make a new list

        if normal===nothing && length(visited)>1
            # find a non-parallel dart and find a normal vector
            for d in visited[end-1:-1:1]
                coplanar(d, dart) || return Loop[]
                if !parallel(d, dart) 
                    normal = normalize(cross(dir(d), dir(dart)))
                    break
                end
            end
            exclude_inner_edges && ( get_n_flat_faces(dart, normal) > 2 && return Loop[] )
        end

        cands = get_next_darts(dart, subset, normal, has_side, exclude_inner_edges)

        length(cands)==0 && return Loop[]
        
        loops = Loop[]
        if !has_side && normal!=nothing && length(cands)==2
            append!(loops, find_loops(visited, cands[1], normal, true, exclude_inner_edges))
            append!(loops, find_loops(visited, cands[2], -normal, true, exclude_inner_edges))
        else
            for cand in cands
                append!(loops, find_loops(visited, cand, normal, has_side, exclude_inner_edges))
            end
        end

        return loops
    end

    return find_loops(Dart[], dart, normal, false, exclude_inner_edges)

end





function get_next_spins(spins::Vector{FaceSpin}, border::Vector{Dart}, subset::Vector{Face})
    tol = 1e-6
    volfaces = [ spin.face for spin in spins ]

    next_spins = FaceSpin[]


    for dart in border

        faces = Face[]
        for face in dart.edge.faces
            # check if faces are not in current volume faces
            face in volfaces && continue

            # target && @show any(face==spin.face for spin in next_spins) 
            # any(face==spin.face for spin in next_spins) && continue #! wrong check location
            
            # check if faces are in the subset
            length(subset)>0 && (face in subset || continue)

            # check if faces has less than 2 volumes
            length(face.volumes)<2 || continue
            
            # check if faces do not have dangling edges
            any(length(dart.edge.faces)==1 for dart in face.loops[1].darts) && continue

            push!(faces, face)
        end

        length(faces)==0 && continue

        # get current spin
        current_spin = nothing
        for spin in spins
            for face in dart.edge.faces
                if spin.face == face
                    current_spin = spin
                    break
                end
            end
        end

        # normals
        n = dir(dart)
        n1 = current_spin.normal # points outwards the volume
        
        # compute angles for each candidate face with the current spin
        angles = Float64[]
        normals = Vec3[]
        for face in faces
            # find loop where dart edge is
            loop = nothing
            loop_idx = 0
            for (i,lo) in enumerate(face.loops)
                if dart.edge in lo.edges
                    loop_idx = i
                    break
                end
            end

            loop = face.loops[loop_idx]

            n2 = getnormal(loop)

            if dart in loop.darts # if dart has the same direction 
                n2 = -n2 # n2 should point outwards the volume
            end
            if loop_idx>1 # hole
                n2 = -n2
            end

            θ = acos(clamp(dot(n1, n2), -1, 1))
            if dot(cross(n1, n2), n) > 0
                angle = π - θ
            else
                angle = π + θ
            end
            push!(angles, angle)
            push!(normals, n2)
        end

        idxs = sortperm(angles)
        faces = faces[idxs]
        normals = normals[idxs]
        spin = FaceSpin(faces[1], normals[1])

        !(spin in next_spins) && push!(next_spins, spin)

    end

    return next_spins
end



function find_face_loops(spin::FaceSpin, subset::Vector{Face} = Face[])
    #! Note: This provides limited volume detection

    # spin: seed face
    # returns a list of faces that define volumes (currently only one volume)

    # check if face does not add a volume
    any(length(dart.edge.faces)==1 for dart in spin.face.loops[1].darts) && return FaceLoop[]

    border = copy(spin.loops[1].darts)

    # add holes
    for lo in spin.face.loops[2:end]
        if norm(getnormal(lo) - spin.normal) < 1e-6
            lo = reverse_loop(lo)
        end
        push!(border, lo.darts...)
    end

    spins = [ spin ]
    visited = [ spin.face ]

    # grow faces
    while length(border)>0

        next_spins = get_next_spins(spins, border, subset)

        length(next_spins)==0 && return FaceLoop[]

        # add faces
        for spin in next_spins
            push!(spins, spin)

            # update border
            spin_darts = [ dart for loop in spin.loops for dart in loop.darts ] # include directed holes
            nborder = length(border)
            nspin = length(spin_darts)
            idxs_border = trues(nborder)
            idxs_spin = trues(nspin)

            for i in 1:length(spin_darts)
                dart = spin_darts[i]
                for j in 1:length(border)
                    if spin_darts[i].edge == border[j].edge
                        idxs_border[j] = false
                        idxs_spin[i] = false
                    end
                end
            end

            border = [ border[idxs_border]; spin_darts[idxs_spin] ]
            length(border)==0 && break

        end
    end

    res = [ FaceLoop(spins) ]
    return res

end
