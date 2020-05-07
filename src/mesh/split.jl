# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# This file includes the code for adding joints between cells

export split!
export generate_joints!
export generate_joints_by_tag!

function joint_shape(shape::ShapeType)
    if shape == LIN2  ; return JLIN2  end
    if shape == LIN3  ; return JLIN3  end
    if shape == LIN4  ; return JLIN4  end
    if shape == TRI3  ; return JTRI3  end
    if shape == TRI6  ; return JTRI6  end
    if shape == QUAD4 ; return JQUAD4 end
    if shape == QUAD8 ; return JQUAD8 end
    if shape == QUAD9 ; return JQUAD9 end
    error("No joint for shape $shape")
end

#=
# Adds joint cells over all shared faces
function generate_joints_old(mesh::Mesh; layers::Int64=2, verbose::Bool=true, tag="", midpointstag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)
    cells  = mesh.elems

    any(c.shape.family==JOINT_SHAPE for c in cells) && error("generate_joints!: mesh already contains joint elements.")
    solids = [ c for c in cells if c.shape.family==SOLID_SHAPE ]

    newpoints = Node[]

    # Splitting: generating new points
    for c in solids
        for (i,p) in enumerate(c.nodes)
            newp = Node([p.coord.x, p.coord.y, p.coord.z])
            push!(newpoints, newp)
            c.nodes[i] = newp
        end
    end

    # List all repeated faces
    face_pairs = Tuple{Cell, Cell}[]

    # Joints generation
    facedict = Dict{UInt64, Cell}()

    # Get paired faces
    for cell in solids
        for face in get_faces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f==nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Generate joint elements
    jcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, tag=tag)
        cell.linked_elems = [f1.oelem, f2.oelem]
        push!(jcells, cell)
    end

    # Generate inner points at joints (used in hydromechanical analyses)
    if layers==3

        auxpointdict = Dict{UInt64,Node}()

        for jcell in jcells
            npts = jcell.shape.basic_shape.npoints
            #npts = div(length(jcell.nodes),2) # number of points in one side
            sample_pts = jcell.nodes[1:npts]
            for p in sample_pts
                hs = hash(p)
                if haskey(auxpointdict, hs)
                    newp = auxpointdict[hs]
                    push!(jcell.nodes, newp)
                else
                    newp = Node(p.coord.x, p.coord.y, p.coord.z, tag=midpointstag)
                    auxpointdict[hs] = newp
                    push!(newpoints, newp)
                    push!(jcell.nodes, newp)
                end
            end

        end

    end

    # Fix JOINT1D_SHAPE cells connectivities
    for c in cells
        c.shape.family != JOINT1D_SHAPE && continue
        scell = c.linked_elems[1]
        nspts = length(scell.nodes)
        c.nodes[1:nspts] .= scell.nodes
    end

    if haskey(mesh.elem_data, "inset-data")
        idata = mesh.elem_data["inset-data"]
        mesh.elem_data["inset-data"] = [ idata; zeros(Int, length(jcells), 3) ]
    end

    # Get points from non-separated cells, e.g. lines, beams, etc.
    points_dict = Dict{UInt64, Node}()
    for c in mesh.elems
        c.shape.family == SOLID_SHAPE && continue # skip because they have points with same coordinates
        c.shape.family == JOINT1D_SHAPE && continue # skip because their points were already considered
        for p in c.nodes
            points_dict[hash(p)] = p
        end
    end

    # Add new cells
    mesh.elems  = vcat(cells, jcells)

    # Get points from solid and joint cells
    mesh.nodes = [ collect(values(points_dict)); newpoints ]

    # update and reorder mesh
    fixup!(mesh, reorder=true)

    # Add field for joints (1 or 0 values)
    ncells = length(mesh.elems)
    joint_data = zeros(Int, ncells, 3) # nlayers, first link, second link
    for i=1:ncells
        cell = mesh.elems[i]
        if cell.shape.family==JOINT_SHAPE
            joint_data[i,1] = layers
            joint_data[i,2] = cell.linked_elems[1].id
            joint_data[i,3] = cell.linked_elems[2].id
        end
    end
    mesh.elem_data["joint-data"] = joint_data

    if verbose
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" length(mesh.nodes)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d new joint cells\n" length(jcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
    end

    return mesh

end
=#


"""
    generate_joints!(mesh, filter=:(); layers=2, tag="", midnodestag="", verbose=true)

Adds joint elements between bulk elements in `mesh`. 
If `filter` is supplied (element tag or expression), the joints are generated over a specific region only.
Joints can be generated with 2 or 3 `layers`.  All generated joints get the supplied `tag`.
Also, all middle points in three-layered joints receive a `midnodestag`.
"""
function generate_joints!(mesh::Mesh, filter::Union{Expr,String}=:(); layers::Int64=2, verbose::Bool=true, tag="", midnodestag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)

    # Target and locked cells
    # locked cells include solids, lines, beams, etc.
    if filter==:()
        targetcells = mesh.elems.solids
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, lockedcells[:joints])
    else
        targetcells = mesh.elems[filter].solids
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, lockedcells[filter][:joints])
    end

    # Splitting: generating new nodes
    newnodes = Node[]
    for c in targetcells
        for (i,p) in enumerate(c.nodes)
            newp = Node([p.coord.x, p.coord.y, p.coord.z])
            push!(newnodes, newp)
            c.nodes[i] = newp
        end
    end

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for cell in targetcells
        for face in get_faces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f==nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Add pairs using surface faces from locked cells
    for face in get_surface(lockedcells)
        hs = hash(face)
        f  = get(facedict, hs, nothing)
        f==nothing && continue
        push!(face_pairs, (face, f))
        delete!(facedict, hs)
    end

    # Generate joint elements
    jointcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, tag=tag)
        cell.linked_elems = [f1.oelem, f2.oelem]
        push!(jointcells, cell)
    end

    # Generate inner nodes at joints (used in hydromechanical analyses)
    if layers==3

        auxnodedict = Dict{UInt64,Node}()

        for jcell in jointcells
            npts = jcell.shape.basic_shape.npoints
            sample_pts = jcell.nodes[1:npts]
            for p in sample_pts
                hs = hash(p)
                if haskey(auxnodedict, hs)
                    newp = auxnodedict[hs]
                    push!(jcell.nodes, newp)
                else
                    newp = Node(p.coord.x, p.coord.y, p.coord.z, tag=midnodestag)
                    auxnodedict[hs] = newp
                    push!(newnodes, newp)
                    push!(jcell.nodes, newp)
                end
            end
        end
    end

    # Fix JOINT1D_SHAPE cells connectivities
    for c in lockedcells
        c.shape.family == JOINT1D_SHAPE || continue
        scell = c.linked_elems[1]
        nspts = length(scell.nodes)
        c.nodes[1:nspts] .= scell.nodes
    end

    if haskey(mesh.elem_data, "inset-data")
        idata = mesh.elem_data["inset-data"]
        mesh.elem_data["inset-data"] = [ idata; zeros(Int, length(jointcells), 3) ]
    end

    # Locked nodes
    lockednodedict = Dict{Int,Node}()
    for cell in lockedcells
        for node in cell.nodes
            node.id == -1 && continue # skip in case of new nodes
            lockednodedict[node.id] = node
        end
    end
    lockednodes = collect(values(lockednodedict))

    # All nodes
    mesh.nodes = vcat(lockednodes, newnodes)

    # All cells
    mesh.elems  = vcat(lockedcells, targetcells, jointcells)

    # Update and reorder mesh
    fixup!(mesh, reorder=true)

    # Add field for joints (1 or 0 values)
    ncells = length(mesh.elems)
    joint_data = zeros(Int, ncells, 3) # nlayers, first link, second link
    for i=1:ncells
        cell = mesh.elems[i]
        if cell.shape.family==JOINT_SHAPE
            joint_data[i,1] = layers
            joint_data[i,2] = cell.linked_elems[1].id
            joint_data[i,3] = cell.linked_elems[2].id
        end
    end
    mesh.elem_data["joint-data"] = joint_data

    if verbose
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" length(mesh.nodes)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d new joint cells\n" length(jointcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
    end

    return mesh

end

# Deprecated function
function split!(mesh::Mesh)
    info("split function was deprecated. Use generate_joints! function instead.")
    generate_joints!(mesh)
end

#=
mutable struct FacePair
    face1::Face
    face2::Face
    idxs1::Array{Int,1}
    idxs2::Array{Int,1}
    FacePair() = new()
end

function generate_joints_candidate!(mesh::Mesh, expr::Expr, tag::String="") # TODO: needs checking
    solids = mesh.elems[:solids][expr]
    @assert length(solids)>0

    # List for all paired faces
    face_pairs = FacePair[]

    # Dict for face pairs: hash(face) => facepair
    fp_dict = Dict{UInt64, FacePair}()

    # Get paired faces
    for cell in solids
        faces_idxs = cell.shape.facet_idxs # vertex connectivities of all faces from cell
        for (i,face) in enumerate(get_faces(cell))
            hs = hash(face)
            fp = get(fp_dict, hs, nothing)
            if fp==nothing # fill first face in fp
                fp = FacePair()
                fp.face1 = face
                fp.idxs1 = faces_idxs[i]
                fp_dict[hs] = fp
            else # fill second face in fp
                fp.face2 = face
                fp.idxs2 = faces_idxs[i]
                push!(face_pairs, fp)
                delete!(fp_dict, hs)
            end
        end
    end

    # Filtering faces
    faces = [ fp.face1 for fp in face_pairs]
    for (i,face) in enumerate(faces); face.id = i end
    faces = faces[expr]

    in_idxs  = [ face.id for face in faces ]
    out_idxs = setdiff(1:length(face_pairs), in_idxs)

    # Generating new points
    for fp in face_pairs[in_idxs]
        for i in fp.idxs1
            p = fp.face1.oelem.nodes[i]
            newp = Node(p.coord.x, p.coord.y, p.coord.z)
            fp.face1.oelem.nodes[i] = newp
        end
        for i in fp.idxs2
            p = fp.face2.oelem.nodes[i]
            newp = Node(p.coord.x, p.coord.y, p.coord.z)
            fp.face2.oelem.nodes[i] = newp
        end
    end

    # Joining extra points
    points_dict = Dict{UInt64,Node}()
    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.oelem.nodes[i]
            points_dict[hash(p)] = p
        end
        for i in fp.idxs2
            p = fp.face2.oelem.nodes[i]
            points_dict[hash(p)] = p
        end
    end

    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.oelem.nodes[i]
            fp.face1.oelem.nodes[i] = points_dict[hash(p)]
        end
        for i in fp.idxs2
            p = fp.face2.oelem.nodes[i]
            fp.face2.oelem.nodes[i] = points_dict[hash(p)]
        end
    end

    # Generating Joints
    jcells = Cell[]
    for fp in face_pairs[in_idxs]
        n   = length(fp.face1.nodes)
        con = Array{Node}(2*n)
        k   = 0
        for i in fp.idxs1
            k += 1
            p1 = fp.face1.oelem.nodes[i]
            h1 = hash(p1)
            for j in fp.idxs2
                p2 = fp.face2.oelem.nodes[j]
                if h1==hash(p2)
                    con[k]   = p1
                    con[n+k] = p2
                    break
                end
            end
        end

        jshape = joint_shape(fp.face1.shape)
        cell = Cell(jshape, con, tag="")
        cell.linked_elems = [fp.face1.oelem, fp.face2.oelem]
        push!(jcells, cell)
    end

    mesh.nodes = collect(Set(p for c in mesh.elems for p in c.nodes))
    mesh.elems  = [mesh.elems; jcells]

    # update and reorder mesh
    fixup!(mesh, reorder=true)

    tag!(jcells, tag)
end
=#


#Generate joints taking account tags

function generate_joints_by_tag!(mesh::Mesh; layers::Int64=2, verbose::Bool=true, tag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)
    cells  = mesh.elems

    any(c.shape.family==JOINT_SHAPE for c in cells) && error("generate_joints!: mesh already contains joint elements.")
    solids = [ c for c in cells if c.shape.family==SOLID_SHAPE ]


    # List all repeated faces
    face_pairs = Tuple{Cell, Cell}[]

    # Joints generation
    facedict = Dict{UInt64, Cell}()

    # Get paired faces
    for cell in solids
        for face in get_faces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f==nothing
                facedict[hs] = face
            else
                if f.tag != face.tag #add to face_pairs only if are different materials. Delete of dictionary of faces them.
                    push!(face_pairs, (face, f))
                    delete!(facedict, hs)
                end
            end
        end
    end

    newpoints = Node[]
    # Splitting: generating new points

    cfp = [fp[1] for fp in face_pairs] #cells of face (pairs) shared between diffenrent materials
    cellsfp = [cfp;cfp]#cells of faces shared between different materials
    pointsdict_fp = Dict{UInt64, Node}()
    for c in cellsfp
        for (i,p) in enumerate(c.nodes)
            hs = hash(p)
            pointsdict_fp[hs] = p
            newp =Node([p.coord.x, p.coord.y, p.coord.z])
            push!(newpoints, newp)
        end
    end

    c_other = collect(values(facedict))
    pointsdict_other = Dict{UInt64, Node}()
    for c in c_other
        for (i,p) in enumerate(c.nodes)
            hs = hash(p)
            f = get(pointsdict_fp,hs,nothing)
            if f==nothing
                pointsdict_fp[hs] = p
                newp =Node([p.coord.x, p.coord.y, p.coord.z])
                push!(newpoints, newp)
            end
        end
    end

    # Generate joint elements
    jcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, tag=tag)
        cell.linked_elems = [f1.oelem, f2.oelem]
        push!(jcells, cell)
    end


    # Get points from non-separated cells
    points_dict = Dict{UInt64, Node}()
    for c in mesh.elems
        c.shape.family == SOLID_SHAPE && continue # skip because they have points with same coordinates
        c.shape.family == JOINT1D_SHAPE && continue # skip because their points were already considered
        for p in c.nodes
            points_dict[hash(p)] = p
        end
    end

    # Add new cells
    mesh.elems  = vcat(cells, jcells)


    # Get points from solid and joint cells
    mesh.nodes = [ collect(values(points_dict)); newpoints ]


    #Creating dictionaries to after assign mesh.nodes to mesh.elems.nodes

   pointsdict_1 =  Dict{UInt64, Node}()
   pointsdict_2 =  Dict{UInt64, Node}()
   pointsdict_3 =  Dict{UInt64, Node}()
   pointsdict_4 =  Dict{UInt64, Node}()

       for p in mesh.nodes
            hs = hash(p)
            f = get(pointsdict_1,hs,nothing)
            if f == nothing
                pointsdict_1[hs] = p
            else
                f = get(pointsdict_2,hs,nothing)
                if f==nothing
                    pointsdict_2[hs] = p
                else
                    f = get(pointsdict_3,hs,nothing)
                    if f==nothing
                        pointsdict_3[hs] = p
                    else
                        f = get(pointsdict_4,hs,nothing)
                        if f==nothing
                            pointsdict_4[hs] = p
                        end
                    end
                end
            end
        end

   pointsdict_1bk =  Dict{UInt64, Node}()
   pointsdict_2bk =  Dict{UInt64, Node}()
   pointsdict_3bk =  Dict{UInt64, Node}()
   pointsdict_4bk =  Dict{UInt64, Node}()

       for p in mesh.nodes
            hs = hash(p)
            f = get(pointsdict_1bk,hs,nothing)
            if f == nothing
                pointsdict_1bk[hs] = p
            else
                f = get(pointsdict_2bk,hs,nothing)
                if f==nothing
                    pointsdict_2bk[hs] = p
                else
                    f = get(pointsdict_3bk,hs,nothing)
                    if f==nothing
                        pointsdict_3bk[hs] = p
                    else
                        f = get(pointsdict_4bk,hs,nothing)
                        if f==nothing
                            pointsdict_4bk[hs] = p
                        end
                    end
                end
            end
        end

    #Copying points from mesh.nodes to mesh.elems[i].nodes if joints

    for c in mesh.elems
        c.shape.family == SOLID_SHAPE && continue
        for (i,p) in enumerate(c.nodes)
            hs = hash(p)
            f = get(pointsdict_4,hs,nothing)
            if f!=nothing
                c.nodes[i] = f
                delete!(pointsdict_4,hs)
            else
                f = get(pointsdict_3,hs,nothing)
                if f!=nothing
                    c.nodes[i] = f
                    delete!(pointsdict_3,hs)
                else
                     f = get(pointsdict_2,hs,nothing)
                    if f!=nothing
                        c.nodes[i] = f
                        delete!(pointsdict_2,hs)
                    else
                        f = get(pointsdict_1,hs,nothing)
                        if f!=nothing
                        c.nodes[i] = f
                        delete!(pointsdict_1,hs)
                        end
                    end
                end
            end
        end
    end


    #Copying points from mesh.nodes to mesh.elems[i].nodes if solids

    for c in mesh.elems
        c.shape.family == JOINT_SHAPE && continue
        for (i,p) in enumerate(c.nodes)
            hs = hash(p)
            f = get(pointsdict_1,hs,nothing)
            if f!=nothing
                c.nodes[i] = f
                delete!(pointsdict_1,hs)
            else
                f = get(pointsdict_4bk,hs,nothing)
                if f!=nothing
                    c.nodes[i] = f
                    delete!(pointsdict_4bk,hs)
                else
                    f = get(pointsdict_3bk,hs,nothing)
                    if f!=nothing
                        c.nodes[i] = f
                        delete!(pointsdict_3bk,hs)
                    else
                        f = get(pointsdict_2bk,hs,nothing)
                        if f!=nothing
                            c.nodes[i] = f
                            delete!(pointsdict_2bk,hs)
                        else
                            f = get(pointsdict_1bk,hs,nothing)
                            if f!=nothing
                                c.nodes[i] = f
                                delete!(pointsdict_1bk,hs)
                            end
                        end
                    end
                end
            end
        end
    end


    # update and reorder mesh
    fixup!(mesh, reorder=true)

    if verbose
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" length(mesh.nodes)
        @printf "  %5d total cells\n" length(mesh.elems)
        @printf "  %5d new joint cells\n" length(jcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
    end

    return mesh

end




