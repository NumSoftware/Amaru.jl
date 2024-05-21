# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export insert_cohesive_elements!


"""
    insert_cohesive_elements!(mesh, filter=nothing; layers=2, tag="", midnodestag="", onlybetweenregions=false, quiet=true)

Adds joint elements between bulk elements in `mesh`.  If `filter` is supplied
(element tag or expression), the joints are generated over a specific region
only.  All generated joints get the supplied `tag`.  If `onlybetweenregions` is
`true`, joints are generated only between regions specified by bulk element
tags.  In this case, the joints get specific tags according to the regions
they are linking.  Joints can be generated with 2 or 3 `layers`. 
For three layers, all middle points in joints receive a `midnodestag`, if provided.
"""
function insert_cohesive_elements!(
    mesh              ::Mesh,
    filter            ::Union{Expr,String,Nothing}=nothing;
    onlybetweenregions::Bool=false,
    autotag           ::Bool=false,
    layers            ::Int64=2,
    tag               ::String="",
    quiet             ::Bool=false,
    midnodestag       ::String=""
)

    

    quiet || printstyled("Insertion of cohesive elements:\n", bold=true, color=:cyan)

    layers in (2,3) || error("insert_cohesive_elements!: wrong number of layers ($layers).")

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
        # remove previous joints at locked region
        lockedcells = setdiff(lockedcells, lockedcells.joints)
    else
        targetcells = mesh.elems[filter].solids
        if length(targetcells)==0
            error("insert_cohesive_elements!: no targetcells found for filter $filter")
        end
        lockedcells = setdiff(mesh.elems, targetcells)
        # remove previous joints at filtered region
        lockedcells = setdiff(lockedcells, lockedcells[filter].joints)
    end

    if !onlybetweenregions
            # Splitting
            # generating new nodes at target cells
            for c in targetcells
                for (i,p) in enumerate(c.nodes)
                    newp = Node(p.coord, tag=p.tag)
                    newp.id = -1
                    c.nodes[i] = newp
                end
            end

            lockedouterfaces = get_outer_facets(lockedcells)

            # Get faces pairs
            facedict = Dict{UInt64, Cell}()
            face_pairs = Tuple{Cell, Cell}[]
            for cell in targetcells
                for face in getfacets(cell)
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
            for face in lockedouterfaces
                hs = hash(face)
                f  = get(facedict, hs, nothing)
                f===nothing && continue
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end

    else # generate joints between tagged regions only 
        # Get tags
        tag_set = Set{String}()
        for cell in targetcells
            # cell.tag != "" && push!(tag_set, cell.tag)
            push!(tag_set, cell.tag)
        end

        # Get joint faces
        trial_faces = Face[]
        for tag in tag_set
            for face in get_outer_facets(targetcells[tag])
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

        # Update joint faces (now with new nodes)
        trial_faces = Face[]
        for tag in tag_set
            for face in get_outer_facets(ocells[tag])
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
        k==n || error("insert_cohesive_elements!: faces f1 and f2 are not coincident.")

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
    # idx = 10000
    # idx = maximum( n.id for n in mesh.nodes )
    for cell in mesh.elems
        for node in cell.nodes
            if node.id<0
                idx += 1
                node.id = idx # new id
            end
            nodesdict[node.id] = node
        end
    end


    # All nodes
    mesh.nodes = collect(values(nodesdict))

    # Update and reorder mesh
    syncronize!(mesh, reorder=true, cleandata=true)


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


function generate_joints_by_tag!(mesh::Mesh; layers::Int64=2, verbose::Bool=true)
    # Get tags
    tag_set = Set{String}()
    for cell in mesh.cells
        push!(tag_set, cell.tag)
    end

    # Get joint faces
    joint_faces = Face[]
    for tag in tag_set
        for face in get_outer_facets(mesh.cells[tag])
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
        for face in get_outer_facets(ocells[tag])
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


function cracksmesh(mesh::Mesh, opening::Real)

    # Get paired faces
    facedict = Dict{UInt64, Cell}()
    face_pairs = Tuple{Cell, Cell}[]
    for cell in mesh.elems
        for face in getfacets(cell)
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


