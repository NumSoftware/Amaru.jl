# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

using TranscodingStreams, CodecZlib

function save_vtk(mesh::AbstractDomain, filename::String; desc::String="")
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
    has_elem_data = !isempty(mesh.elem_data)

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
            for i in 1:npoints
                for j in 1:ncomps
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
            for i in 1:ncells
                for j in 1:ncomps
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


function get_array_node!(array::AbstractArray, name::String, compressed, buf)
    dtype = string(eltype(array))
    ncomps = size(array,2)
    if compressed # appends compressed data to buf
        xdata = XmlElement("DataArray", attributes=("type"=>dtype, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"appended", "offset"=>"$(position(buf))"))
        level = 4 # compression level
        inipos = position(buf)

        array = collect(transpose(array))

        # temporary header
        write(buf, UInt64(1), UInt64(0), UInt64(0), UInt64(0))
        arr_size = length(array)*sizeof(eltype(array))

        # write compressed array
        zbuf = ZlibCompressorStream(buf, level=level)
        write(zbuf, array)
        write(zbuf, TranscodingStreams.TOKEN_END)
        flush(zbuf)
        TranscodingStreams.finalize(zbuf.codec) # finalize codec
        
        # rewrite header
        endpos = position(buf)
        comp_arr_size = endpos - inipos - 4*sizeof(UInt64(0)) # considering header size
        seek(buf, inipos)
        write(buf, UInt64(1), UInt64(arr_size), UInt64(arr_size), UInt64(comp_arr_size))
        seek(buf, endpos)
    else
        xdata = XmlElement("DataArray", attributes=("type"=>dtype, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"ascii"))
        isfloat = eltype(array)<:AbstractFloat
        io = IOBuffer()
        nrows = size(array,1)
        for i in 1:nrows
            for j in 1:ncomps
                if isfloat
                    @printf io "%20.10e" Float32(array[i,j])
                else
                    print(io, array[i,j], "  ")
                end
            end
            i<nrows && print(io, "\n")
        end
        xdata.content = String(take!(io))
    end
    return xdata
end


function save_vtu(mesh::AbstractDomain, filename::String; desc::String="", compress=false)
    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    root = XmlElement("VTKFile", attributes=("type"=>"UnstructuredGrid", "version"=>"1.0",  "byte_order"=>"LittleEndian", "header_type"=>"UInt64", "compressor"=>"vtkZLibDataCompressor"))
    ugrid = XmlElement("UnstructuredGrid")
    piece = XmlElement("Piece", attributes=("NumberOfPoints"=>"$npoints", "NumberOfCells"=>"$ncells"))
    addchild!(ugrid, piece)
    addchild!(root, ugrid)

    buf = IOBuffer() # only for compressed data

    # Write coordinates
    xpoints = XmlElement("Points")
    coords = Float64[ node.coord[i] for node in mesh.nodes, i in 1:3 ]

    addchild!(xpoints, get_array_node!(coords, "Points", compress, buf))
    addchild!(piece, xpoints)

    xcells = XmlElement("Cells")

    # Write connectivities
    conn = Int32[]
    for cell in mesh.elems
        for node in cell.nodes
            push!(conn, node.id-1)
        end
    end
    addchild!(xcells, get_array_node!(conn, "connectivity", compress, buf))

    # Write offset
    offsets = Int32[]
    offset = 0
    for cell in mesh.elems
        offset += length(cell.nodes)
        push!(offsets, offset)
    end
    # append_compressed_array!(buf, offsets, level)
    addchild!(xcells, get_array_node!(offsets, "offsets", compress, buf))

    # Write cell types
    types = Int32[]
    for cell in mesh.elems
        push!(types, Int32(cell.shape.vtk_type))
    end
    addchild!(xcells, get_array_node!(types, "types", compress, buf))
    addchild!(piece, xcells)

    # Write node data
    has_node_data = !isempty(mesh.node_data)
    if has_node_data
        xpointdata = XmlElement("PointData")
        for (field, D) in mesh.node_data
            isempty(D) && continue
            addchild!(xpointdata, get_array_node!(D, field, compress, buf))
        end
        addchild!(piece, xpointdata)
        # push!(piece.children, xpointdata)
    end

    # Write cell data
    has_elem_data = !isempty(mesh.elem_data)
    if has_elem_data
        xcelldata = XmlElement("CellData")
        for (field, D) in mesh.elem_data
            isempty(D) && continue
            addchild!(xcelldata, get_array_node!(D, field, compress, buf))
        end
        addchild!(piece, xcelldata)
    end

    if compress
        xappended = XmlElement("AppendedData", attributes=("encoding"=>"raw",))
        xappended.content = take!(buf)
        addchild!(root, xappended)
    end

    fileatts = ("version"=>"1.0",)
    doc = XmlDocument(fileatts)
    addchild!(doc, XmlComment(desc))
    addchild!(doc, root)
    save(doc, filename)
end


function add_extra_fields!(mesh::AbstractDomain)
    # Add field for joints
    if any( c.shape.family==JOINTCELL for c in mesh.elems )
        ncells = length(mesh.elems)
        joint_data = zeros(Int, ncells, 3) # nlayers, first link, second link
        for i in 1:ncells
            cell = mesh.elems[i]
            nlayers = cell.shape.name[1:2]=="J3" ? 3 : 2
            if cell.shape.family==JOINTCELL
                joint_data[i,1] = nlayers
                joint_data[i,2] = cell.linked_elems[1].id
                joint_data[i,3] = cell.linked_elems[2].id
            end
        end
        mesh.elem_data["joint-data"] = joint_data
    end

    # Add field for embedded nodes
    if any( c.shape.family==LINEJOINTCELL for c in mesh.elems )
        ncells = length(mesh.elems)
        inset_data = zeros(Int, ncells, 3) # npoints, first link id, second link id
        for i in 1:ncells
            cell = mesh.elems[i]
            if cell.shape.family==LINEJOINTCELL
                inset_data[i,1] = cell.shape.npoints
                inset_data[i,2] = cell.linked_elems[1].id
                inset_data[i,3] = cell.linked_elems[2].id
            end
        end
        mesh.elem_data["inset-data"] = inset_data
    end
end


"""
    save(mesh, filename, quiet=true)

Saves a mesh object into a file. Available formats are vtu and vtk.
"""
function save(mesh::AbstractDomain, filename::String; compress=false, quiet=true)
    formats = (".vtk", ".vtu") 
    _, format = splitext(filename)
    format in formats || error("save: Cannot save $(typeof(mesh)) to $filename. Available formats are $formats.")

    add_extra_fields!(mesh)
    desc = "File generated by Amaru Finite Element Code"

    if format==".vtk"
        save_vtk(mesh, filename, desc=desc)
    elseif format==".vtu"
        save_vtu(mesh, filename, desc=desc, compress=compress)
    end
    quiet || printstyled( "  file $filename saved \e[K \n", color=:cyan)
    return nothing
end


function read_vtk(filename::String)
    # read nodal information
    alltext = read(filename, String)
    data    = split(alltext)

    npoints = 0
    ncells  = 0
    coords  = zeros(0,0)
    connects = Array{Int,1}[]
    cell_types = Int[]

    node_data = OrderedDict{String,Array}()
    elem_data = OrderedDict{String,Array}()

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
            for i in 1:npoints
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
            for i in 1:ncells
                npts = parse(Int64, data[idx+1])
                idx += 1
                conn = Int[]
                for j in 1:npts
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
            for i in 1:ncells
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
            for i in 1:npoints
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
            for i in 1:ncells
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
            for i in 1:npoints
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
            for i in 1:ncells
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
    doc = XmlDocument(filename)

    # check if file is compressed
    if getchild(doc.root, "AppendedData")!==nothing 
        raw = Vector{UInt8}(strip(getchild(doc.root, "AppendedData").content[2:end])) # remove first character _
        nodes = getallchildren(doc.root)

        # update content on nodes with offset att
        for node in nodes
            node isa XmlElement || continue
            node.name == "DataArray" && haskey(node.attributes, "offset") || continue
            offset = parse(Int, node.attributes["offset"])

            first  = offset + 1
            last   = offset + 4*sizeof(UInt64(0))
            
            header = Int.(reinterpret(UInt64, raw[first:last]))
            len    = header[4]

            first = offset + 4*sizeof(UInt64(0)) + 1
            last  = first + len - 1

            node.content = transcode(ZlibDecompressor, raw[first:last]) # decompression
        end
    end

    piece   = doc.root("UnstructuredGrid", "Piece")
    npoints = parse(Int, piece.attributes["NumberOfPoints"])
    ncells  = parse(Int, piece.attributes["NumberOfCells"])
    coords  = get_array(piece("Points", "DataArray"))

    xcells     = piece("Cells")
    conn       = get_array(xcells["Name"=>"connectivity"][1]) .+ 1
    offsets    = get_array(xcells["Name"=>"offsets"][1])
    cell_types = get_array(xcells["Name"=>"types"][1])

    connects = Array{Int,1}[]
    pos = 1
    for off in offsets
        push!(connects, conn[pos:off])
        pos = off+1
    end

    node_data = OrderedDict{String,Array}()
    elem_data = OrderedDict{String,Array}()

    xpointdata = piece("PointData")
    if xpointdata!==nothing
        for array_data in xpointdata.children
            label = array_data.attributes["Name"]
            node_data[label] = get_array(array_data)
        end
    end

    xcelldata = piece("CellData")
    if xcelldata!==nothing
        for array_data in xcelldata.children
            label = array_data.attributes["Name"]
            elem_data[label] = get_array(array_data)
        end
    end

    return Mesh(coords, connects, cell_types, node_data, elem_data)
end


function get_array(data_array::XmlElement)
    TYPES = Dict("Float32"=>Float32, "Float64"=>Float64, "Int32"=>Int32, "Int64"=>Int64)

    ncomps = haskey(data_array.attributes, "NumberOfComponents") ? parse(Int, data_array.attributes["NumberOfComponents"]) : 1
    dtype  = TYPES[data_array.attributes["type"]]
    isbin  = data_array.content isa Vector{UInt8}
    if isbin
        if ncomps==1
            return reinterpret(dtype, data_array.content)
        else
            return transpose(reshape(reinterpret(dtype, data_array.content), ncomps, :))
        end
    else # string
        if ncomps==1
            return parse.(dtype, split(data_array.content))
        else
            return transpose(reshape(parse.(dtype, split(data_array.content)), ncomps, :))
        end
    end
end


# Setting a Mesh object
function Mesh(coords, connects, vtk_types, node_data, elem_data)

    npoints = size(coords,1)
    ncells  = length(connects)

    # Setting points
    nodes = Node[]
    for i in 1:npoints
        X = coords[i,:]
        node = Node(X)
        node.id = i
        push!(nodes, node)
    end

    # Mesh object
    ndim = getndim(nodes)

    # check for exceptional 3d cases
    if haskey(node_data, "rx") || haskey(node_data, "ry")
        ndim = 3
    end

    mesh = Mesh(ndim)
    mesh.nodes = nodes

    # Setting cells
    has_polyvertex = false

    for i in 1:ncells
        conn = mesh.nodes[ connects[i] ]
        vtk_shape = VTKCellType(vtk_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYVERTEX
            has_polyvertex = true
        else
            shape = get_shape_from_vtk( vtk_shape, length(conn), ndim )
        end
        cell  = Cell(shape, conn)
        cell.id = i
        push!(mesh.elems, cell)
    end

    # update mesh and get faces and edges
    compute_facets!(mesh)

    # Setting data
    mesh.node_data = node_data
    mesh.elem_data  = elem_data
    # mesh.node_data = merge(mesh.node_data, node_data)
    # mesh.elem_data  = merge(mesh.elem_data, elem_data)

    # Fix information for 1D joints
    if haskey(mesh.elem_data, "inset-data")
        inset_data = mesh.elem_data["inset-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape==POLYVERTEX && inset_data[i,1]>0
                linked_ids = inset_data[i,2:3]
                cell.linked_elems = mesh.elems[linked_ids]
                cell.linked_elems[1].crossed = true # host cell is crossed
                n = length(cell.linked_elems[2].nodes)
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
            end
        end
    end

    # Fix information for joints
    if haskey(mesh.elem_data, "joint-data")
        joint_data = mesh.elem_data["joint-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape==POLYVERTEX && joint_data[i,1] in (2,3)
                nlayers    = joint_data[i,1]+0
                linked_ids = joint_data[i,2:3]
                cell.linked_elems = mesh.elems[linked_ids]
                n = length(cell.nodes)
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, nlayers)
            end
        end
    end

    # Flip cells
    for cell in mesh.elems
        isinverted(cell) && flip!(cell)
    end

    syncronize!(mesh)

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
            if cell.shape == POLYVERTEX
                n = length(cell.nodes)
                # look for joints1D and fix shape
                if n>=5
                    # check if cell is related to a JLINK2
                    hss = hash( [ cell.nodes[i] for i in 1:n-2] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK2
                        hs0   = hash( [ cell.nodes[i] for i in n-1:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_elems = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end

                    # check if cell is related to a JLINK3
                    hss = hash( [ cell.nodes[i] for i in 1:n-3] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK3
                        hs0   = hash( [ cell.nodes[i] for i in n-2:n] )
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
                    for i in 1:div(n,2)
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
                    for i in 1:stride
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
        has_joints = any( C -> C.shape.family==JOINTCELL, mesh.elems )

        if has_joints
            # generate dict of faces
            facedict = Dict{UInt64, Cell}()
            for cell in mesh.elems
                for face in getfacets(cell)
                    hs = hash(face)
                    #f  = get(facedict, hs, nothing)
                    facedict[hs] = face
                end
            end

            for cell in mesh.elems
                if cell.shape.family == JOINTCELL
                    n = length(cell.nodes)
                    hs1 = hash( [ cell.nodes[i] for i in 1:div(n,2)] )
                    hs2 = hash( [ cell.nodes[i] for i in div(n,2)+1:n] )
                    cell1 = facedict[hs1].owner
                    cell2 = facedict[hs2].owner
                    cell.linked_elems = [ cell1, cell2 ]
                end
            end
        end
    end
    =#

    #=

    # Fix linked_elems for embedded elements (line cells not connected to other elements)
    has_line = any(C->C.shape.family==LINECELL, mesh.elems)
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

    for i in 1:npoints
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

    for i in 1:ncells
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

    for i in 1:npoints
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

    for i in 1:ncells
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

    for i in 1:nfaces
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

    for i in 1:nedges
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
    Mesh(filename)

Constructs a `Mesh` object based on a file.
"""
function Mesh(filename::String; reorder=false, quiet=true)
    
    formats = (".vtk", ".vtu", ".tetgen")

    quiet || printstyled("Mesh loading: filename $filename\n", bold=true, color=:cyan)

    _, format = splitext(filename)

    format in formats || error("Mesh: cannot read format \"$format\". Suitable formats are $formats.")

    if format==".vtk"
        quiet || print("  Reading VTK legacy format...\n")
        mesh = read_vtk(filename)
    elseif format==".vtu"
        quiet || print("  Reading VTU format...\n")
        mesh = read_vtu(filename)
    elseif format==".tetgen"
        quiet || print("  Reading tetgen output files...\n")
        mesh = read_tetgen(filename)
    end

    quiet || printstyled( "  file $filename loaded \e[K \n", color=:cyan)

    # Reorder nodal numbering
    if reorder
        quiet || print("  reordering points...\r")
        reorder!(mesh)
    end

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        println("  ", mesh.ctx.ndim, "d                   ")
        @printf "  %5d points\n" npoints
        @printf "  %5d cells\n" ncells
        @printf "  %5d faces\n" nfaces
        @printf "  %5d surface edges\n" nedges
    end

    return mesh
end
