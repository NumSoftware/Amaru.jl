# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru
#

# This file includes the code for adding joints between cells

export split!
export generate_joints!
export generate_joints_by_tag!


"""
    generate_joints!(mesh, filter=nothing; layers=2, tag="", midnodestag="", onlybetweenregions=false, quiet=true)

Adds joint elements between bulk elements in `mesh`.  If `filter` is supplied
(element tag or expression), the joints are generated over a specific region
only.  All generated joints get the supplied `tag`.  If `onlybetweenregions` is
`true`, joints are generated only between regions specified by bulk element
tags.  In this case, the joints get specific tags according to the regions
they are linking.  Joints can be generated with 2 or 3 `layers`. 
For three layers, all middle points in joints receive a `midnodestag`, if provided.
"""
function generate_joints!(
    mesh              ::Mesh,
    filter            ::Union{Expr,String,Nothing}=nothing;
    onlybetweenregions::Bool=false,
    autotag           ::Bool=false,
    layers            ::Int64=2,
    tag               ::String="",
    quiet             ::Bool=false,
    midnodestag       ::String=""
)

    quiet || printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)

    layers in (2,3) || error("generate_joints!: wrong number of layers ($layers).")

    joint_shape_dict = Dict(
                            (LIN2 ,2) => JLIN2,
                            (LIN3 ,2) => JLIN3,
                            (LIN4 ,2) => JLIN4,
                            (TRI3 ,2) => JTRI3,
                            (TRI6 ,2) => JTRI6,
                            (QUAD4,2) => JQUAD4,
                            (QUAD8,2) => JQUAD8,
                            (LIN2 ,3) => J3LIN2,
                            (LIN3 ,3) => J3LIN3,
                            (LIN4 ,3) => J3LIN4,
                            (TRI3 ,3) => J3TRI3,
                            (TRI6 ,3) => J3TRI6,
                            (QUAD4,3) => J3QUAD4,
                            (QUAD8,3) => J3QUAD8,
                           )
    onlybetweenregions && (autotag=true)

    # Target and locked cells
    # locked cells include solids, lines, beams, etc.
    if filter===nothing
        targetcells = mesh.elems.solids
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, lockedcells.joints)
    else
        targetcells = mesh.elems[filter].solids
        if length(targetcells)==0
            error("generate_joints!: no targetcells found for filter $filter")
        end
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, lockedcells[filter].joints)
    end

    if !onlybetweenregions
        # Splitting: generating new nodes
        for c in targetcells
            for (i,p) in enumerate(c.nodes)
                newp = Node(p.coord, tag=p.tag)
                newp.id = -1
                c.nodes[i] = newp
            end
        end

        # Get paired faces
        facedict = Dict{UInt64, Cell}()
        face_pairs = Tuple{Cell, Cell}[]
        for cell in targetcells
            for face in getfaces(cell)
                hs = hash(face)
                f  = get(facedict, hs, nothing)
                if f===nothing
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
            f===nothing && continue
            push!(face_pairs, (face, f))
            delete!(facedict, hs)
        end
    else # generate joints between tagged regions only TODO:check
        # Get tags
        tag_set = Set{String}()
        for cell in targetcells
            # cell.tag != "" && push!(tag_set, cell.tag)
            push!(tag_set, cell.tag)
        end

        # @show tag_set

        # Get joint faces
        trial_faces = Face[]
        for tag in tag_set
            for face in get_surface(targetcells[tag])
                push!(trial_faces, face)
            end
        end

        # @show length(trial_faces)

        # Get nodes to duplicate
        nodes_to_dup = Set{Node}()
        facedict = Dict{UInt64, Cell}()
        ocells   = Cell[] # face owner cells
        for face in trial_faces
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
                facedict[hs] = face
            else # a matching face was found
                
                push!(ocells, face.owner)
                push!(ocells, f.owner)
                for node in face.nodes
                    push!(nodes_to_dup, node)
                end
                delete!(facedict, hs)
            end
        end

        # @show length(facedict)
        # @show length(nodes_to_dup)

        # Duplicate nodes
        for cell in targetcells
            for (i,node) in enumerate(cell.nodes)
                if node in nodes_to_dup
                    newnode = Node(node.coord)
                    newnode.id = -1
                    cell.nodes[i] = newnode
                end
            end
        end

        # Join nodes per tag
        for tag in tag_set
            nodedict = Dict{UInt64, Node}()
            for cell in targetcells[tag]
                for (i,node) in enumerate(cell.nodes)
                    hs = hash(node)
                    n  = get(nodedict, hs, nothing)
                    if n===nothing
                        nodedict[hs] = node
                    else
                        cell.nodes[i] = n
                    end
                end
            end
        end

        # Get joint faces (now with new nodes)
        trial_faces = Face[]
        for tag in tag_set
            for face in get_surface(ocells[tag])
                push!(trial_faces, face)
            end
        end

        # Get paired faces
        facedict = Dict{UInt64, Cell}()
        face_pairs = Tuple{Cell, Cell}[]
        for face in trial_faces
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Generate joint elements
    jointcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.nodes)
        con = Array{Node}(undef, 2*n)
        k = 0
        for (i,p1) in enumerate(f1.nodes)
            for p2 in f2.nodes
                if hash(p1)==hash(p2)
                    k += 1
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end
        k==n || error("generate_joints!: faces f1 and f2 are not coincident.")

        #jshape = joint_shape(f1.shape)
        jshape = joint_shape_dict[(f1.shape,layers)]

        # if tag=="" && autotag
        if autotag
            tagA, tagB = sort([f1.owner.tag, f2.owner.tag]) 
            if tagA==tagB
                tag = "joint-"*tagA
            else
                tag = "joint-"*tagA*"-"*tagB
            end
            # tag = join( sort([f1.owner.tag, f2.owner.tag]), "-" )
        end
        cell = Cell(jshape, con, tag=tag)
        cell.linked_elems = [f1.owner, f2.owner]
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
                    push!(jcell.nodes, newp)
                end
            end
        end
    end

    # Fix cells connectivities for special interface elements
    for c in lockedcells
        c.shape.family in (TIPJOINTCELL, LINEJOINTCELL) || continue
        scell = c.linked_elems[1]
        nspts = length(scell.nodes)
        c.nodes[1:nspts] .= scell.nodes
    end

    if haskey(mesh.elem_data, "inset-data")
        idata = mesh.elem_data["inset-data"]
        mesh.elem_data["inset-data"] = [ idata; zeros(Int, length(jointcells), 3) ]
    end

    # All cells
    mesh.elems  = vcat(lockedcells, targetcells, jointcells)

    # Nodes dict
    nodesdict = Dict{Int,Node}()
    idx = length(mesh.nodes)
    for cell in mesh.elems
        for node in cell.nodes
            if node.id==-1
                idx += 1
                node.id = idx # new id
            end
            nodesdict[node.id] = node
        end
    end

    # All nodes
    mesh.nodes = collect(values(nodesdict))

    # Update and reorder mesh
    fixup!(mesh, reorder=true)


    if !quiet
        @printf "  %4dd mesh                             \n" mesh.env.ndim
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

function generate_joints_by_tag!(mesh::Mesh; layers::Int64=2, verbose::Bool=true)
    # Get tags
    tag_set = Set{String}()
    for cell in mesh.cells
        push!(tag_set, cell.tag)
    end

    # Get joint faces
    joint_faces = Face[]
    for tag in tag_set
        for face in get_surface(mesh.cells[tag])
            push!(joint_faces, face)
        end
    end

    # Get nodes to duplicate
    nodes_to_dup = Set{Node}()
    for face in joint_faces
        for node in face
            push!(nodes_to_dup, node)
        end
    end

    # List of owner cells
    ocells = Cell[ f.owner for f in joint_faces ]

    # Duplicate nodes
    for cell in ocells
        for (i,node) in enumerate(cell.nodes)
            if node in nodes_to_dup
                cell.nodes[i] = copy(node)
            end
        end
    end

    # Get joint faces (now with new nodes)
    joint_faces = Face[]
    for tag in tag_set
        for face in get_surface(ocells[tag])
            push!(joint_faces, face)
        end
    end

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for face in joint_faces
        hs = hash(face)
        f  = get(facedict, hs, nothing)
        if f===nothing
            facedict[hs] = face
        else
            push!(face_pairs, (face, f))
            delete!(facedict, hs)
        end
    end

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
    solids = mesh.elems.solids[expr]
    @assert length(solids)>0

    # List for all paired faces
    face_pairs = FacePair[]

    # Dict for face pairs: hash(face) => facepair
    fp_dict = Dict{UInt64, FacePair}()

    # Get paired faces
    for cell in solids
        faces_idxs = cell.shape.facet_idxs # vertex connectivities of all faces from cell
        for (i,face) in enumerate(getfaces(cell))
            hs = hash(face)
            fp = get(fp_dict, hs, nothing)
            if fp===nothing # fill first face in fp
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
            p = fp.face1.owner.nodes[i]
            newp = Node(p.coord.x, p.coord.y, p.coord.z)
            fp.face1.owner.nodes[i] = newp
        end
        for i in fp.idxs2
            p = fp.face2.owner.nodes[i]
            newp = Node(p.coord.x, p.coord.y, p.coord.z)
            fp.face2.owner.nodes[i] = newp
        end
    end

    # Joining extra points
    points_dict = Dict{UInt64,Node}()
    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.owner.nodes[i]
            points_dict[hash(p)] = p
        end
        for i in fp.idxs2
            p = fp.face2.owner.nodes[i]
            points_dict[hash(p)] = p
        end
    end

    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.owner.nodes[i]
            fp.face1.owner.nodes[i] = points_dict[hash(p)]
        end
        for i in fp.idxs2
            p = fp.face2.owner.nodes[i]
            fp.face2.owner.nodes[i] = points_dict[hash(p)]
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
            p1 = fp.face1.owner.nodes[i]
            h1 = hash(p1)
            for j in fp.idxs2
                p2 = fp.face2.owner.nodes[j]
                if h1==hash(p2)
                    con[k]   = p1
                    con[n+k] = p2
                    break
                end
            end
        end

        jshape = joint_shape(fp.face1.shape)
        cell = Cell(jshape, con, tag="")
        cell.linked_elems = [fp.face1.owner, fp.face2.owner]
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

function generate_joints_by_tag_2!(mesh::Mesh; layers::Int64=2, verbose::Bool=true, tag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)
    cells  = mesh.elems

    any(c.shape.family==JOINTCELL for c in cells) && error("generate_joints!: mesh already contains joint elements.")
    solids = [ c for c in cells if c.shape.family==BULKCELL ]


    # List all repeated faces
    face_pairs = Tuple{Cell, Cell}[]

    # Joints generation
    facedict = Dict{UInt64, Cell}()

    # Get paired faces
    for cell in solids
        for face in getfaces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
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
            if f===nothing
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
        cell.linked_elems = [f1.owner, f2.owner]
        push!(jcells, cell)
    end


    # Get points from non-separated cells
    points_dict = Dict{UInt64, Node}()
    for c in mesh.elems
        c.shape.family == BULKCELL && continue # skip because they have points with same coordinates
        c.shape.family == LINEJOINTCELL && continue # skip because their points were already considered
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
                if f===nothing
                    pointsdict_2[hs] = p
                else
                    f = get(pointsdict_3,hs,nothing)
                    if f===nothing
                        pointsdict_3[hs] = p
                    else
                        f = get(pointsdict_4,hs,nothing)
                        if f===nothing
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
                if f===nothing
                    pointsdict_2bk[hs] = p
                else
                    f = get(pointsdict_3bk,hs,nothing)
                    if f===nothing
                        pointsdict_3bk[hs] = p
                    else
                        f = get(pointsdict_4bk,hs,nothing)
                        if f===nothing
                            pointsdict_4bk[hs] = p
                        end
                    end
                end
            end
        end

    #Copying points from mesh.nodes to mesh.elems[i].nodes if joints

    for c in mesh.elems
        c.shape.family == BULKCELL && continue
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
        c.shape.family == JOINTCELL && continue
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
        @printf "  %4dd mesh                             \n" mesh.env.ndim
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


function cracksmesh(mesh::Mesh, opening::Real)

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for cell in mesh.elems
        for face in getfaces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f===nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Get normals and distances
    U = mesh.node_data["U"]
    crack_faces = Cell[]
    for pair in face_pairs
        face1, face2 = pair

        #error()
        X1 = face1.nodes[1].coord
        X2 = face1.nodes[2].coord
        X3 = face1.nodes[3].coord
        n = cross(X2-X1, X3-X1)
        normalize!(n)
        nnodes = length(face1.nodes)
        node_map = [node.id for node in face1.nodes]
        U1 = mean(U[node_map,:], dims=1)
        node_map = [node.id for node in face2.nodes]
        U2 = mean(U[node_map,:], dims=1)

        #dn = maximum((U2-U1)*n) # normal distance
        dn = dot(U2-U1,n) # normal distance
        #d  = norm(mean(U2-U1, dims=1)) # total distance

        #U2 = mean(U2, dims=1)
        #U1 = mean(U1, dims=1)
        d  = norm(U2-U1) # total distance
        #d = maximum(norm.(eachrow(U2-U1)))
        if dn>0
            #display(U2-U1)
            #error()
        end

        if dn>0 && d>=opening
            push!(crack_faces, face1)
        end
    end

    nodes = getnodes(crack_faces)
    ids = [ node.id for node in nodes ]

    newsmesh = Mesh(crack_faces)

    for (k,v) in mesh.node_data
        k=="id" && continue
        newsmesh.node_data[k] = v[ids]
    end

    return newsmesh
end


