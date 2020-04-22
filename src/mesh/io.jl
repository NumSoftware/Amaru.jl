# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function save_vtk(mesh::AbstractMesh, filename::String; desc::String="")
    # Saves a UnstructuredGrid
    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)

    # Number of total connectivities
    nconns = 0
    for cell in mesh.elems
        nconns += 1 + length(cell.nodes)
    end

    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, desc)
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write nodes
    for (i,node) in enumerate(mesh.nodes)
        @printf f "%23.15e %23.15e %23.15e \n" node.coord.x node.coord.y node.coord.z
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ", ncells, " ", nconns)
    for cell in mesh.elems
        print(f, length(cell.nodes), " ")
        for node in cell.nodes
            print(f, node.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write elem types
    println(f, "CELL_TYPES ", ncells)
    for cell in mesh.elems
        println(f, Int(cell.shape.vtk_type))
    end
    println(f)

    has_node_data = !isempty(mesh.node_data)
    has_elem_data  = !isempty(mesh.elem_data)

    # Write node data
    if has_node_data
        println(f, "POINT_DATA ", npoints)
        for (field,D) in mesh.node_data
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat
            dtype = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==1
                println(f, "SCALARS $field $dtype 1")
                println(f, "LOOKUP_TABLE default")
            else
                println(f, "VECTORS ", "$field $dtype")
            end
            for i=1:npoints
                for j=1:ncomps
                    if isfloat
                        @printf f "%23.10e" Float32(D[i,j])
                    else
                        @printf f "%10d" D[i,j]
                    end
                end
            end
            println(f)
        end
    end

    # Write cell data
    if has_elem_data
        println(f, "CELL_DATA ", ncells)
        for (field,D) in mesh.elem_data
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat

            dtype = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==1
                println(f, "SCALARS $field $dtype 1")
                println(f, "LOOKUP_TABLE default")
            else
                println(f, "VECTORS ", "$field $dtype")
            end
            for i=1:ncells
                for j=1:ncomps
                    if isfloat
                        @printf f "%23.10e" Float32(D[i,j])
                    else
                        @printf f "%10d" D[i,j]
                    end
                end
            end
            println(f)
        end
    end

    close(f)

    return nothing
end

function save_vtu(mesh::AbstractMesh, filename::String; desc::String="")
    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    root = Xnode("VTKFile", Dict("type"=>"UnstructuredGrid", "version"=>"0.1", "byte_order"=>"LittleEndian"))
    ugrid = Xnode("UnstructuredGrid")
    piece = Xnode("Piece", Dict("NumberOfPoints"=>"$npoints", "NumberOfCells"=>"$ncells"))
    push!(ugrid.children, piece)
    push!(root.children, ugrid)

    io = IOBuffer()

    # Write coordinates
    xpoints = Xnode("Points")
    xcoords  = Xnode("DataArray", Dict("type"=>"Float64", "NumberOfComponents"=>"3", "format"=>"ascii"))
    for (i,node) in enumerate(mesh.nodes)
        @printf io "%17.7e %17.7e %17.7e" node.coord.x node.coord.y node.coord.z
        i<npoints && print(io, "\n")
    end
    xcoords.content = String(take!(io))
    push!(xpoints.children, xcoords)
    push!(piece.children, xpoints)

    xcells = Xnode("Cells")

    # Write connectivities
    xconn  = Xnode("DataArray", Dict("type"=>"Int32", "Name"=>"connectivity", "format"=>"ascii"))
    for cell in mesh.elems
        for node in cell.nodes
            print(io, node.id-1, "  ")
        end
    end
    xconn.content = String(take!(io))

    # Write offset
    xoffset = Xnode("DataArray", Dict("type"=>"Int32", "Name"=>"offsets", "format"=>"ascii"))
    offset = 0
    for cell in mesh.elems
        offset += length(cell.nodes)
        print(io, offset, "  ")
    end
    xoffset.content = String(take!(io))

    # Write cell types
    xtypes = Xnode("DataArray", Dict("type"=>"Int32", "Name"=>"types", "format"=>"ascii"))
    for cell in mesh.elems
        print(io, Int(cell.shape.vtk_type), "  ")
    end
    xtypes.content = String(take!(io))

    push!(xcells.children, xconn)
    push!(xcells.children, xoffset)
    push!(xcells.children, xtypes)

    push!(piece.children, xcells)

    # Node and Cell data
    has_node_data = !isempty(mesh.node_data)
    has_elem_data  = !isempty(mesh.elem_data)

    # Write node data
    if has_node_data
        xpointdata = Xnode("PointData")
        for (field,D) in mesh.node_data
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat
            dtype = isfloat ? "Float64" : "Int32"
            ncomps = size(D,2)
            xdata = Xnode("DataArray", Dict("type"=>dtype, "Name"=>"$field", "NumberOfComponents"=>"$ncomps", "format"=>"ascii"))
            for i=1:npoints
                for j=1:ncomps
                    if isfloat
                        @printf io "%17.7e" Float32(D[i,j])
                    else
                        print(io, D[i,j], "  ")
                    end
                end
            end
            xdata.content = String(take!(io))
            push!(xpointdata.children, xdata)
        end
        push!(piece.children, xpointdata)
    end

    # Write cell data
    if has_elem_data
        xcelldata = Xnode("CellData")
        for (field,D) in mesh.elem_data
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat
            dtype = isfloat ? "Float64" : "Int32"
            ncomps = size(D,2)
            xdata = Xnode("DataArray", Dict("type"=>dtype, "Name"=>"$field", "NumberOfComponents"=>"$ncomps", "format"=>"ascii"))
            for i=1:ncells
                for j=1:ncomps
                    if isfloat
                        @printf io "%20.10e" Float32(D[i,j])
                    else
                        print(io, D[i,j], "  ")
                    end
                end
            end
            xdata.content = String(take!(io))
            push!(xcelldata.children, xdata)
        end
        push!(piece.children, xcelldata)
    end

    fileatts = OrderedDict("version"=>"1.0")
    doc = Xdoc(fileatts, root)
    save(doc, filename)
end


function read_vtk(filename::String)
    # Reading file
    # ============

    # read nodal information
    alltext = read(filename, String)
    data    = split(alltext)

    #local npoints, ncells
    #local coords, connects, cell_types

    npoints = 0
    ncells  = 0
    coords  = zeros(0,0)
    connects = Array{Int,1}[]
    cell_types = Int[]

    node_data = OrderedDict{String,Array}()
    elem_data  = OrderedDict{String,Array}()

    reading_node_data = false
    reading_elem_data  = false

    TYPES = Dict("float32"=>Float32, "float64"=>Float64, "int"=>Int64)

    idx = 1
    while idx<=length(data)
        if data[idx] == "DATASET"
            gridtype = data[idx+1]
            gridtype == "UNSTRUCTURED_GRID" || error("load_VTK_unstructured_grid: this reader only support files of VTK UNSTRUCTURED_GRID")
        end

        # read coords
        if data[idx] == "POINTS"
            npoints = parse(Int64, data[idx+1]) # read number of coords
            coords  = zeros(npoints,3)
            idx += 2
            for i=1:npoints
                coords[i,1] = parse(Float64, data[idx+1])
                coords[i,2] = parse(Float64, data[idx+2])
                coords[i,3] = parse(Float64, data[idx+3])
                idx += 3
            end
        end

        # read cells connectivities
        if data[idx] == "CELLS"
            ncells = parse(Int64, data[idx+1])
            ncdata = parse(Int64, data[idx+2])
            idx += 2

            connects = Array{Int,1}[]
            for i=1:ncells
                npts = parse(Int64, data[idx+1])
                idx += 1
                conn = Int[]
                for j=1:npts
                    idx += 1
                    id = parse(Int64, data[idx]) + 1
                    push!(conn, id)
                end
                push!(connects, conn)
            end
        end

        # read type of cells
        if data[idx] == "CELL_TYPES"
            idx += 1
            cell_types = Int[]
            for i=1:ncells
                idx += 1
                vtk_shape = parse(Int64, data[idx])
                push!(cell_types, vtk_shape)
            end
        end

        if data[idx] == "POINT_DATA"
            idx += 1
            reading_node_data = true
            reading_elem_data  = false
        end

        if data[idx] == "CELL_DATA"
            idx += 1
            reading_elem_data  = true
            reading_node_data = false
        end

        if data[idx] == "VECTORS" && reading_node_data
            label = data[idx+1]
            ty = data[idx+2]
            dtype = TYPES[ty]
            idx += 2
            vectors = zeros(dtype, npoints,3)
            for i=1:npoints
                vectors[i,1] = parse(dtype, data[idx+1])
                vectors[i,2] = parse(dtype, data[idx+2])
                vectors[i,3] = parse(dtype, data[idx+3])
                idx += 3
            end
            node_data[label] = vectors
        end

        if data[idx] == "VECTORS" && reading_elem_data
            label = data[idx+1]
            ty = data[idx+2]
            dtype = TYPES[ty]
            idx += 2
            vectors = zeros(dtype, ncells,3)
            for i=1:ncells
                vectors[i,1] = parse(dtype, data[idx+1])
                vectors[i,2] = parse(dtype, data[idx+2])
                vectors[i,3] = parse(dtype, data[idx+3])
                idx += 3
            end
            elem_data[label] = vectors
        end

        if data[idx] == "SCALARS" && reading_node_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], npoints)
            for i=1:npoints
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            node_data[label] = scalars
        end

        if data[idx] == "SCALARS" && reading_elem_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], ncells)
            for i=1:ncells
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            elem_data[label] = scalars
        end

        idx += 1

    end

    return Mesh(coords, connects, cell_types, node_data, elem_data)
end


function read_vtu(filename::String)
    # Reading file
    # ============

    TYPES = Dict("Float32"=>Float32, "Float64"=>Float64, "Int32"=>Int32, "Int64"=>Int64)

    doc = Xdoc(filename)
    piece = doc.root("UnstructuredGrid", "Piece")
    npoints = parse(Int, piece.attributes["NumberOfPoints"])
    ncells  = parse(Int, piece.attributes["NumberOfCells"])
    strcoords = piece("Points", "DataArray").content
    coords = transpose(reshape(parse.(Float64, split(strcoords)), 3, :))

    xmlcells = piece("Cells")
    conn     = parse.(Int, split(xmlcells["Name"=>"connectivity"][1].content)) .+ 1
    offsets  = parse.(Int, split(xmlcells["Name"=>"offsets"][1].content))
    cell_types = parse.(Int, split(xmlcells["Name"=>"types"][1].content))

    connects = Array{Int,1}[]
    pos = 1
    for off in offsets
        push!(connects, conn[pos:off])
        pos = off+1
    end

    node_data = OrderedDict{String,Array}()
    elem_data  = OrderedDict{String,Array}()

    xpointdata = piece("PointData")
    if xpointdata!=nothing
        for arr in xpointdata.children
            ncomps = parse(Int, arr.attributes["NumberOfComponents"])
            dtype = TYPES[arr.attributes["type"]]
            label = arr.attributes["Name"]
            if ncomps==1
                node_data[label] = parse.(dtype, split(arr.content))
            else
                node_data[label] = transpose(reshape(parse.(dtype, split(arr.content)), ncomps, npoints))
            end
        end
    end

    xcelldata = piece("CellData")
    if xcelldata!=nothing
        for arr in xcelldata.children
            ncomps = parse(Int, arr.attributes["NumberOfComponents"])
            dtype = TYPES[arr.attributes["type"]]
            label = arr.attributes["Name"]
            if ncomps==1
                elem_data[label] = parse.(dtype, split(arr.content))
            else
                elem_data[label] = transpose(reshape(parse.(dtype, split(arr.content)), ncomps, ncells))
            end
        end
    end

    return Mesh(coords, connects, cell_types, node_data, elem_data)
end


# Setting a Mesh object
# =====================
function Mesh(coords, connects, vtk_types, node_data, elem_data)

    mesh = Mesh()
    npoints = size(coords,1)
    ncells  = length(connects)

    # Setting points
    for i=1:npoints
        X = coords[i,:]
        node = Node(X)
        node.id = i
        push!(mesh.nodes, node)
    end

    # Set ndim
    ndim = 1
    for node in mesh.nodes
        node.coord.y != 0.0 && (ndim=2)
        node.coord.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim

    # Setting cells
    has_polyvertex = false

    for i=1:ncells
        conn = mesh.nodes[ connects[i] ]
        vtk_shape = VTKCellType(vtk_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYV
            has_polyvertex = true
        else
            shape = get_shape_from_vtk( vtk_shape, length(conn), ndim )
        end
        cell  = Cell(shape, conn)
        push!(mesh.elems, cell)
    end

    # update mesh and get faces and edges
    fixup!(mesh, reorder=false)

    # Setting data
    mesh.node_data = merge(mesh.node_data, node_data)
    mesh.elem_data  = merge(mesh.elem_data, elem_data)

    # Fix information for 1D joints
    if haskey(mesh.elem_data, "inset-data")
        inset_data = mesh.elem_data["inset-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape.family==JOINT1D_SHAPE
                linked_ids = inset_data[i,2:3]
                cell.linked_elems = mesh.elems[linked_ids]
                cells.linked_elems[1].crossed = true # host cell is crossed
            end
        end
    end

    # Fix information for joints
    if haskey(mesh.elem_data, "joint-data")
        joint_data = mesh.elem_data["joint-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape.family==JOINT_SHAPE
                nlayers    = joint_data[i,1]
                linked_ids = joint_data[i,2:3]
                cell.linked_elems = mesh.elems[linked_ids]
                n = length(cell.nodes)
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, nlayers)
            end
        end
    end

    # remaining polyvertex cells
    #cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)

    #=
    # Fix shape for polyvertex cells
    if has_polyvertex
        # mount dictionary of cells
        cdict = Dict{UInt64, Cell}()
        for cell in mesh.elems
            hs = hash(cell)
            cdict[hs] = cell
        end

        # check cells
        for cell in mesh.elems
            if cell.shape == POLYV
                n = length(cell.nodes)
                # look for joints1D and fix shape
                if n>=5
                    # check if cell is related to a JLINK2
                    hss = hash( [ cell.nodes[i] for i=1:n-2] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK2
                        hs0   = hash( [ cell.nodes[i] for i=n-1:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_elems = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end

                    # check if cell is related to a JLINK3
                    hss = hash( [ cell.nodes[i] for i=1:n-3] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK3
                        hs0   = hash( [ cell.nodes[i] for i=n-2:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_elems = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end
                end

                # look for conventional joints and fix shape
                if n%2==0 && n>=4
                    isjoint = true
                    for i=1:div(n,2)
                        p1 = cell.nodes[i]
                        p2 = cell.nodes[div(n,2)+i]
                        if hash(p1) != hash(p2)
                            isjoint = false
                            break
                        end
                    end
                    if isjoint
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, 2)
                        continue
                    end
                end

                # look for joint elements with 3 layers and fix shape
                if n%3==0 && n>=6
                    isjoint = true
                    stride= div(n,3)
                    delta = 0.0
                    for i=1:stride
                        p1 = cell.nodes[i]
                        p2 = cell.nodes[stride+i]
                        if hash(p1) != hash(p2)
                            isjoint = false
                            break
                        end
                    end
                    if isjoint
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, 3)
                        continue
                    end
                end

                # remaining polyvertex cells
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
            end
        end

        # Update linked cells in joints

        # check if there are joints
        has_joints = any( C -> C.shape.family==JOINT_SHAPE, mesh.elems )

        if has_joints
            # generate dict of faces
            facedict = Dict{UInt64, Cell}()
            for cell in mesh.elems
                for face in get_faces(cell)
                    hs = hash(face)
                    #f  = get(facedict, hs, nothing)
                    facedict[hs] = face
                end
            end

            for cell in mesh.elems
                if cell.shape.family == JOINT_SHAPE
                    n = length(cell.nodes)
                    hs1 = hash( [ cell.nodes[i] for i=1:div(n,2)] )
                    hs2 = hash( [ cell.nodes[i] for i=div(n,2)+1:n] )
                    cell1 = facedict[hs1].oelem
                    cell2 = facedict[hs2].oelem
                    cell.linked_elems = [ cell1, cell2 ]
                end
            end
        end
    end

    # Fix linked_elems for embedded elements (line cells not connected to other elements)
    has_line = any(C->C.shape.family==LINE_SHAPE, mesh.elems)
    if has_line
        # TODO: find the owner of orphan line cells OR use parent id information
    end
    =#

    return mesh

end


function read_tetgen(filekey::String)
    # reads files .node and .ele

    local points, cells, cell_types, npoints, ncells
    node_data = Dict{String,Array}()
    #point_vector_data = Dict{String,Array}()
    elem_data  = Dict{String,Array}()

    # read nodal information

    f=open(filekey*".node")
    line = readline(f)
    npoints, ndim, _, _ = parse.(Int, split(line))
    points  = zeros(npoints,3)

    for i=1:npoints
        line = readline(f)
        items = split(line)
        points[i,1] = parse(Float64, items[2])
        points[i,2] = parse(Float64, items[3])
        points[i,3] = parse(Float64, items[4])
    end
    close(f)

    # read cell information

    f=open(filekey*".ele")
    line = readline(f)
    ncells, npts, _ = parse.(Int, split(line))
    cells = Array{Int,1}[]

    for i=1:ncells
        line = readline(f)
        items = split(line)
        pts = parse.(Int, items[2:end])
        pts .+= 1  # convert to one-based
        push!(cells, pts)
    end
    close(f)

    # cell types
    #
    cell_types = 10*ones(Int, ncells)


    #ugrid = UnstructuredGrid("Unstructured grid from tetgen", points, cells, cell_types)
    #return ugrid
    # TODO: Return a Mesh object
end


# Read a mesh from TetGen output files
# TODO: deprecate this function. Consider read_ugrid_tetgen function
function tetreader(filekey::String)
    # reading .node file
    f = open(filekey*".node")
    data = readlines(f)
    npoints, ndim, att, hasmarker = parse.(split(data[1]))
    points = Array{Node,1}(npoints)

    for i=1:npoints
        pdata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? pdata[end] : ""
        node = Node(pdata[2], pdata[3], pdata[4], tag)
        node.id = i
        points[i] = node
    end

    # reading .ele file
    f = open(filekey*".ele")
    data = readlines(f)
    ncells, ntetpts, hasatt = parse.(split(data[1]))
    celltype = ntetpts==4 ? TET4 : TET10
    cells = Array{Cell,1}(ncells)

    for i=1:ncells
        cdata = parse.(split(data[i+1]))
        tag   = hasatt==1 ? cdata[end] : 0
        cellpoints = points[cdata[2:ntetpts+1]]
        cell = Cell(celltype, cellpoints, tag)
        cell.id = i
        cells[i] = cell
    end

    # reading .face file
    f = open(filekey*".face")
    data = readlines(f)
    nfaces, hasmarker = parse.(split(data[1]))
    nfacepts = ntetpts==4 ? 3 : 6
    celltype = ntetpts==4 ? TRI3 : TRI6
    faces = Array{Cell,1}(nfaces)

    for i=1:nfaces
        fdata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? fdata[end] : 0
        facepts = points[fdata[2:nfacepts+1]]
        nfacepts==6 && (facepts = facepts[[1,2,3,6,4,5]])
        faces[i] = Cell(celltype, facepts, tag)
    end

    # reading .edge file
    f    = open(filekey*".edge")
    data = readlines(f)
    nedges, hasmarker = parse.(split(data[1]))
    nedgepts = ntetpts==4 ? 2 : 3
    celltype = ntetpts==4 ? LIN2 : LIN3
    edges = Array{Cell,1}(nedges)

    for i=1:nedges
        edata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? edata[end] : 0
        edgepts = points[edata[2:nedgepts+1]]
        edges[i] = Cell(celltype, edgepts, tag)
    end

    mesh = Mesh()
    mesh.nodes = points
    mesh.elems  = cells
    mesh.faces  = faces
    mesh.edges  = edges

    return mesh
end



"""
    save(mesh, filename, verbose=true)

Saves a mesh object into a file in VTK legacy format
"""
function save(mesh::AbstractMesh, filename::String; verbose::Bool=true, silent::Bool=false)
    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    format = split(filename, ".")[end]

    if     format=="vtk" ; save_vtk(mesh, filename, desc="File generated by Amaru")
    elseif format=="vtu"; save_vtu(mesh, filename, desc="File generated by Amaru")
    else   error("save: Cannot save $filename. Available formats are vtk and vtu.")
    end
    verbosity>0 && printstyled( "  file $filename written \033[K \n", color=:cyan)
end


"""
    Mesh(filename)

Constructs a `Mesh` object based on a file in VTK legacy format or JSON format.
"""
function Mesh(filename::String; verbose::Bool=false, silent::Bool=false, format::String="", reorder::Bool=false)

    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    verbosity>1 && printstyled("Mesh loading: filename $filename\n", bold=true, color=:cyan)

    if format==""
        basename, format = splitext(filename)
        format = format[2:end]
    end

    if format=="vtk"
        verbosity>1 && print("  Reading VTK legacy format...\n")
        mesh = read_vtk(filename)
    elseif format=="vtu"
        verbosity>1 && print("  Reading VTU format...\n")
        mesh = read_vtu(filename)
    elseif format=="json"
        verbosity>1 && print("  Reading JSON format...\n")
        mesh = read_json(filename)
    elseif format=="tetgen"
        verbosity>1 && print("  Reading tetgen output files...\n")
        mesh = read_tetgen(filename)
    else
        error("Mesh: Reading $format format is not available (file=$filename)")
    end

    verbosity>0 && printstyled( "  file $filename loaded \033[K \n", color=:cyan)

    # Reorder nodal numbering
    if reorder
        verbosity>0 && print("  reordering points...\r")
        reorder!(mesh)
    end

    if verbosity>1
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        println("  ", mesh.ndim, "d                   ")
        @printf "  %5d points\n" npoints
        @printf "  %5d cells\n" ncells
        @printf "  %5d faces\n" nfaces
        @printf "  %5d surface edges\n" nedges
    end

    return mesh
end
