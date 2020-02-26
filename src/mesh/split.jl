# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

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

# Adds joint cells over all shared faces
function generate_joints!(mesh::Mesh; layers::Int64=2, verbose::Bool=true, tag="", midpointstag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)
    cells  = mesh.cells

    any(c.shape.family==JOINT_SHAPE for c in cells) && error("generate_joints!: mesh already contains joint elements.")
    solids = [ c for c in cells if c.shape.family==SOLID_SHAPE ]

    newpoints = Point[]

    # Splitting: generating new points
    for c in solids
        for (i,p) in enumerate(c.points)
            newp = Point([p.x, p.y, p.z])
            push!(newpoints, newp)
            c.points[i] = newp
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
        n   = length(f1.points)
        con = Array{Point}(undef, 2*n)
        for (i,p1) in enumerate(f1.points)
            for p2 in f2.points
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, tag=tag)
        cell.linked_cells = [f1.ocell, f2.ocell]
        push!(jcells, cell)
    end

    # Generate inner points at joints (used in hydromechanical analyses)
    if layers==3

        auxpointdict = Dict{UInt64,Point}()

        for jcell in jcells
            npts = jcell.shape.basic_shape.npoints
            #npts = div(length(jcell.points),2) # number of points in one side
            sample_pts = jcell.points[1:npts]
            for p in sample_pts
                hs = hash(p)
                if haskey(auxpointdict, hs)
                    newp = auxpointdict[hs]
                    push!(jcell.points, newp)
                else
                    newp = Point(p.x, p.y, p.z, tag=midpointstag)
                    auxpointdict[hs] = newp
                    push!(newpoints, newp)
                    push!(jcell.points, newp)
                end
            end

        end

    end

    # Fix JOINT1D_SHAPE cells connectivities
    for c in cells
        c.shape.family != JOINT1D_SHAPE && continue
        scell = c.linked_cells[1]
        nspts = length(scell.points)
        c.points[1:nspts] .= scell.points
    end


    # Get points from non-separated cells
    points_dict = Dict{UInt64, Point}()
    for c in mesh.cells
        c.shape.family == SOLID_SHAPE && continue # skip because they have points with same coordinates
        c.shape.family == JOINT1D_SHAPE && continue # skip because their points were already considered
        for p in c.points
            points_dict[hash(p)] = p
        end
    end

    # Add new cells
    mesh.cells  = vcat(cells, jcells)

    # Get points from solid and joint cells
    mesh.points = [ collect(values(points_dict)); newpoints ]

    # update and reorder mesh
    fixup!(mesh, reorder=true)

    if verbose
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" length(mesh.points)
        @printf "  %5d total cells\n" length(mesh.cells)
        @printf "  %5d new joint cells\n" length(jcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
        println("  done.")
    end

    return mesh

end

# Deprecated function
function split!(mesh::Mesh)
    info("split function was deprecated. Use generate_joints! function instead.")
    generate_joints!(mesh)
end

mutable struct FacePair
    face1::Face
    face2::Face
    idxs1::Array{Int,1}
    idxs2::Array{Int,1}
    FacePair() = new()
end

function generate_joints_candidate!(mesh::Mesh, expr::Expr, tag::String="") # TODO: needs checking
    solids = mesh.cells[:solids][expr]
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
            p = fp.face1.ocell.points[i]
            newp = Point(p.x, p.y, p.z)
            fp.face1.ocell.points[i] = newp
        end
        for i in fp.idxs2
            p = fp.face2.ocell.points[i]
            newp = Point(p.x, p.y, p.z)
            fp.face2.ocell.points[i] = newp
        end
    end

    # Joining extra points
    points_dict = Dict{UInt64,Point}()
    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.ocell.points[i]
            points_dict[hash(p)] = p
        end
        for i in fp.idxs2
            p = fp.face2.ocell.points[i]
            points_dict[hash(p)] = p
        end
    end

    for fp in face_pairs[out_idxs]
        for i in fp.idxs1
            p = fp.face1.ocell.points[i]
            fp.face1.ocell.points[i] = points_dict[hash(p)]
        end
        for i in fp.idxs2
            p = fp.face2.ocell.points[i]
            fp.face2.ocell.points[i] = points_dict[hash(p)]
        end
    end

    # Generating Joints
    jcells = Cell[]
    for fp in face_pairs[in_idxs]
        n   = length(fp.face1.points)
        con = Array{Point}(2*n)
        k   = 0
        for i in fp.idxs1
            k += 1
            p1 = fp.face1.ocell.points[i]
            h1 = hash(p1)
            for j in fp.idxs2
                p2 = fp.face2.ocell.points[j]
                if h1==hash(p2)
                    con[k]   = p1
                    con[n+k] = p2
                    break
                end
            end
        end

        jshape = joint_shape(fp.face1.shape)
        cell = Cell(jshape, con, tag="")
        cell.linked_cells = [fp.face1.ocell, fp.face2.ocell]
        push!(jcells, cell)
    end

    mesh.points = collect(Set(p for c in mesh.cells for p in c.points))
    mesh.cells  = [mesh.cells; jcells]

    # update and reorder mesh
    fixup!(mesh, reorder=true)

    tag!(jcells, tag)
end


 #Generate joints taking account tags

function generate_joints_by_tag!(mesh::Mesh; layers::Int64=2, verbose::Bool=true, tag="")

    verbose && printstyled("Mesh generation of joint elements:\n", bold=true, color=:cyan)
    cells  = mesh.cells

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

    newpoints = Point[]
    # Splitting: generating new points

    cfp = [fp[1] for fp in face_pairs] #cells of face (pairs) shared between diffenrent materials
    cellsfp = [cfp;cfp]#cells of faces shared between different materials
    pointsdict_fp = Dict{UInt64, Point}()
    for c in cellsfp
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            pointsdict_fp[hs] = p
            newp =Point([p.x, p.y, p.z])
            push!(newpoints, newp)
        end
    end

    c_other = collect(values(facedict))
    pointsdict_other = Dict{UInt64, Point}()
    for c in c_other
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            f = get(pointsdict_fp,hs,nothing)
            if f==nothing
                pointsdict_fp[hs] = p
                newp =Point([p.x, p.y, p.z])
                push!(newpoints, newp)
            end
        end
    end

    # Generate joint elements
    jcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.points)
        con = Array{Point}(undef, 2*n)
        for (i,p1) in enumerate(f1.points)
            for p2 in f2.points
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, tag=tag)
        cell.linked_cells = [f1.ocell, f2.ocell]
        push!(jcells, cell)
    end


    # Get points from non-separated cells
    points_dict = Dict{UInt64, Point}()
    for c in mesh.cells
        c.shape.family == SOLID_SHAPE && continue # skip because they have points with same coordinates
        c.shape.family == JOINT1D_SHAPE && continue # skip because their points were already considered
        for p in c.points
            points_dict[hash(p)] = p
        end
    end

    # Add new cells
    mesh.cells  = vcat(cells, jcells)


    # Get points from solid and joint cells
    mesh.points = [ collect(values(points_dict)); newpoints ]


    #Creating dictionaries to after assign mesh.points to mesh.cells.points

   pointsdict_1 =  Dict{UInt64, Point}()
   pointsdict_2 =  Dict{UInt64, Point}()
   pointsdict_3 =  Dict{UInt64, Point}()
   pointsdict_4 =  Dict{UInt64, Point}()

       for p in mesh.points
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

   pointsdict_1bk =  Dict{UInt64, Point}()
   pointsdict_2bk =  Dict{UInt64, Point}()
   pointsdict_3bk =  Dict{UInt64, Point}()
   pointsdict_4bk =  Dict{UInt64, Point}()

       for p in mesh.points
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

    #Copying points from mesh.points to mesh.cells[i].points if joints

    for c in mesh.cells
        c.shape.family == SOLID_SHAPE && continue
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            f = get(pointsdict_4,hs,nothing)
            if f!=nothing
                c.points[i] = f
                delete!(pointsdict_4,hs)
            else
                f = get(pointsdict_3,hs,nothing)
                if f!=nothing
                    c.points[i] = f
                    delete!(pointsdict_3,hs)
                else
                     f = get(pointsdict_2,hs,nothing)
                    if f!=nothing
                        c.points[i] = f
                        delete!(pointsdict_2,hs)
                    else
                        f = get(pointsdict_1,hs,nothing)
                        if f!=nothing
                        c.points[i] = f
                        delete!(pointsdict_1,hs)
                        end
                    end
                end
            end
        end
    end


    #Copying points from mesh.points to mesh.cells[i].points if solids

    for c in mesh.cells
        c.shape.family == JOINT_SHAPE && continue
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            f = get(pointsdict_1,hs,nothing)
            if f!=nothing
                c.points[i] = f
                delete!(pointsdict_1,hs)
            else
                f = get(pointsdict_4bk,hs,nothing)
                if f!=nothing
                    c.points[i] = f
                    delete!(pointsdict_4bk,hs)
                else
                    f = get(pointsdict_3bk,hs,nothing)
                    if f!=nothing
                        c.points[i] = f
                        delete!(pointsdict_3bk,hs)
                    else
                        f = get(pointsdict_2bk,hs,nothing)
                        if f!=nothing
                            c.points[i] = f
                            delete!(pointsdict_2bk,hs)
                        else
                            f = get(pointsdict_1bk,hs,nothing)
                            if f!=nothing
                                c.points[i] = f
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
        @printf "  %5d points\n" length(mesh.points)
        @printf "  %5d total cells\n" length(mesh.cells)
        @printf "  %5d new joint cells\n" length(jcells)
        nfaces = length(mesh.faces)
        nfaces>0 && @printf("  %5d faces\n", nfaces)
        nedges = length(mesh.edges)
        nedges>0 && @printf("  %5d edges\n", nedges)
        println("  done.")
    end

    return mesh

end




