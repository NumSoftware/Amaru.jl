# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

function save_vtk(mesh::Mesh, filename::String; desc::String="")
    # Saves a UnstructuredGrid
    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    # Number of total connectivities
    nconns = 0
    for cell in mesh.cells
        nconns += 1 + length(cell.points)
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
    for (i,point) in enumerate(mesh.points)
        @printf f "%23.15e %23.15e %23.15e \n" point.x point.y point.z
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ", ncells, " ", nconns)
    for cell in mesh.cells
        print(f, length(cell.points), " ")
        for point in cell.points
            print(f, point.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write elem types
    println(f, "CELL_TYPES ", ncells)
    for cell in mesh.cells
        println(f, Int(cell.shape.vtk_type))
    end
    println(f)

    has_point_data = !isempty(mesh.point_data)
    has_cell_data  = !isempty(mesh.cell_data)

    # Write point data
    if has_point_data
        println(f, "POINT_DATA ", npoints)
        for (field,D) in mesh.point_data
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
    if has_cell_data
        println(f, "CELL_DATA ", ncells)
        for (field,D) in mesh.cell_data
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

function save_vtu(mesh::Mesh, filename::String; desc::String="")
    # Saves a UnstructuredGrid
    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    tab1 = "  "
    tab2 = tab1^2
    tab3 = tab1^3
    tab4 = tab1^4
    tab5 = tab1^5

    # Open filename
    f = open(filename, "w")

    println(f, """<?xml version="1.0"?>""")
    println(f, """<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">""")
    println(f, """<UnstructuredGrid>""")
    print(f, tab1)
    println(f, """<Piece NumberOfPoints="$npoints" NumberOfCells="$ncells">""")

    # Write points
    print(f, tab2)
    println(f, """<Points>""")
    print(f, tab3)
    println(f, """<DataArray type="Float64" NumberOfComponents="3" format="ascii">""")
    for (i,point) in enumerate(mesh.points)
        print(f, tab4)
       @printf f "%23.15e %23.15e %23.15e \n" point.x point.y point.z
    end
    print(f, tab3)
    println(f, """<\\DataArray>""")
    print(f, tab2)
    println(f, """<\\Points>""")

    # Write cells
    print(f, tab2)
    println(f, """<Cells>""")
    print(f, tab3)
    println(f, """<DataArray type="Int32" Name="connectivity" format="ascii">""")

    # Write connectivities
    for cell in mesh.cells
        for point in cell.points
            print(f, tab4)
            print(f, point.id-1, " ")
        end
        println(f)
    end
    print(f, tab3)
    println(f, """<\\DataArray>""")

    # Offsets
    print(f, tab3)
    println(f, """<DataArray type="Int32" Name="offsets" format="ascii">""")
    offset = 0
    for cell in mesh.cells
        offset += length(cell.points)
        print(f, tab4)
        println(f, "        ", offset)
    end
    print(f, tab3)
    println(f, """<\\DataArray>""")

    # Offsets
    print(f, tab3)
    println(f, """      <DataArray type="Int32" Name="types" format="ascii">""")
    for cell in mesh.cells
        print(f, tab4)
        println(f, Int(cell.shape.vtk_type))
    end
    print(f, tab3)
    println(f, """<\\DataArray>""")
    print(f, tab2)
    println(f, """    <\\Cells>""")

    # Point data
    has_point_data = !isempty(mesh.point_data)
    has_cell_data  = !isempty(mesh.cell_data)

    println(f, """      <DataArray type="Float64" NumberOfComponents="3" format="ascii">""")

    # Write point data
    if has_point_data
        println(f, tab2, """<PointData>""")
        if has_point_data
            for (field,D) in mesh.point_data
                isempty(D) && continue
                isfloat = eltype(D)<:AbstractFloat
                dtype = isfloat ? "Float64" : "Int32"
                ncomps = size(D,2)
                println(f, tab3, """<DataArray type="$dtype" Name="$field" NumberOfComponents="$ncomps" format="ascii">""")
                for i=1:npoints
                    print(f, tab4)
                    for j=1:ncomps
                        @printf f "%23.10e" isfloat ? Float32(D[i,j]) : D[i,j]
                    end
                end
                println(f, tab3, """<\\DataArray>""")
            end
        end
    end
    println(f, tab2, """<\\PointData>""")

    # Write cell data
    if has_point_data
        println(f, tab2, """<CellData>""")
        if has_cell_data
            for (field,D) in mesh.cell_data
                isempty(D) && continue
                isfloat = eltype(D)<:AbstractFloat
                dtype = isfloat ? "Float64" : "Int32"
                ncomps = size(D,2)
                println(f, tab3, """<DataArray type="$dtype" Name="$field" NumberOfComponents="$ncomps" format="ascii">""")
                for i=1:npoints
                    print(f, tab4)
                    for j=1:ncomps
                        @printf f "%23.10e" isfloat ? Float32(D[i,j]) : D[i,j]
                    end
                end
                println(f, tab3, """<\\DataArray>""")
            end
        end
    end
    println(f, tab2, """<\\PointData>""")

    println(f, """<\\UnstructuredGrid>""")
    println(f, """<\\VTKFile""")

end

function read_vtk(filename::String)
    # Reading file
    # ============

    # read nodal information
    alltext = read(filename, String)
    data    = split(alltext)

    local coords, connects, cell_types, npoints, ncells
    point_data = OrderedDict{String,Array}()
    cell_data  = OrderedDict{String,Array}()
    reading_point_data = false
    reading_cell_data  = false

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
            reading_point_data = true
            reading_cell_data  = false
        end

        if data[idx] == "CELL_DATA"
            idx += 1
            reading_cell_data  = true
            reading_point_data = false
        end

        if data[idx] == "VECTORS" && reading_point_data
            label = data[idx+1]
            idx += 2
            vectors = zeros(npoints,3)
            for i=1:npoints
                vectors[i,1] = parse(Float64, data[idx+1])
                vectors[i,2] = parse(Float64, data[idx+2])
                vectors[i,3] = parse(Float64, data[idx+3])
                idx += 3
            end
            point_data[label] = vectors
        end

        if data[idx] == "SCALARS" && reading_point_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], npoints)
            for i=1:npoints
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            point_data[label] = scalars
        end

        if data[idx] == "SCALARS" && reading_cell_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], ncells)
            for i=1:ncells
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            cell_data[label] = scalars
        end

        idx += 1

    end

    # Setting a Mesh object
    # =====================

    mesh = Mesh()

    # Setting points
    for i=1:npoints
        X = coords[i,:]
        point = Point(X)
        point.id = i
        push!(mesh.points, point)
    end

    # Set ndim
    ndim = 1
    for point in mesh.points
        point.y != 0.0 && (ndim=2)
        point.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim

    # Setting cells
    has_polyvertex = false

    for i=1:ncells
        conn = mesh.points[ connects[i] ]
        vtk_shape = VTKCellType(cell_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYV
            has_polyvertex = true
        else
            shape = get_shape_from_vtk( vtk_shape, length(conn), ndim )
        end
        cell  = Cell(shape, conn)
        push!(mesh.cells, cell)
    end

    # update mesh and get faces and edges
    fixup!(mesh, reorder=false)

    # Setting data
    mesh.point_data = merge(mesh.point_data, point_data)
    mesh.cell_data  = merge(mesh.cell_data, cell_data)

    # Fix shape for polyvertex cells
    if has_polyvertex
        # mount dictionary of cells
        cdict = Dict{UInt64, Cell}()
        for cell in mesh.cells
            hs = hash(cell)
            cdict[hs] = cell
        end

        # check cells
        for cell in mesh.cells
            if cell.shape == POLYV
                n = length(cell.points)
                # look for joints1D and fix shape
                if n>=5
                    # check if cell is related to a JLINK2
                    hss = hash( [ cell.points[i] for i=1:n-2] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK2
                        hs0   = hash( [ cell.points[i] for i=n-1:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_cells = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end

                    # check if cell is related to a JLINK3
                    hss = hash( [ cell.points[i] for i=1:n-3] )
                    if haskey(cdict, hss)
                        cell.shape = JLINK3
                        hs0   = hash( [ cell.points[i] for i=n-2:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_cells = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end
                end

                # look for conventional joints and fix shape
                if n%2==0 && n>=4
                    isjoint = true
                    for i=1:div(n,2)
                        p1 = cell.points[i]
                        p2 = cell.points[div(n,2)+i]
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
                    stride= div(n,3)
                    delta = 0.0
                    for i=1:stride
                        p1 = cell.points[i]
                        p2 = cell.points[stride+i]
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
        has_joints = any( C -> C.shape.family==JOINT_SHAPE, mesh.cells )

        if has_joints
            # generate dict of faces
            facedict = Dict{UInt64, Cell}()
            for cell in mesh.cells
                for face in get_faces(cell)
                    hs = hash(face)
                    f  = get(facedict, hs, nothing)
                    facedict[hs] = face
                end
            end

            for cell in mesh.cells
                if cell.shape.family == JOINT_SHAPE
                    n = length(cell.points)
                    hs1 = hash( [ cell.points[i] for i=1:div(n,2)] )
                    hs2 = hash( [ cell.points[i] for i=div(n,2)+1:n] )
                    cell1 = facedict[hs1].ocell
                    cell2 = facedict[hs2].ocell
                    cell.linked_cells = [ cell1, cell2 ]
                end
            end
        end
    end

    # Fix linked_cells for embedded elements (line cells not connected to other elements)
    has_line = any(C->C.shape.family==LINE_SHAPE, mesh.cells)
    if has_line
        # TODO: find the owner of orphan line cells OR use parent id information
    end

    return mesh

end


function read_tetgen(filekey::String)
    # reads files .node and .ele

    local points, cells, cell_types, npoints, ncells
    point_data = Dict{String,Array}()
    #point_vector_data = Dict{String,Array}()
    cell_data  = Dict{String,Array}()

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
    points = Array{Point,1}(npoints)

    for i=1:npoints
        pdata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? pdata[end] : ""
        point = Point(pdata[2], pdata[3], pdata[4], tag)
        point.id = i
        points[i] = point
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
    mesh.points = points
    mesh.cells  = cells
    mesh.faces  = faces
    mesh.edges  = edges

    return mesh
end



"""
    save(mesh, filename, verbose=true)

Saves a mesh object into a file in VTK legacy format
"""
function save(mesh::Mesh, filename::String; verbose::Bool=true)
    save_vtk(mesh, filename, desc="File generated by FemMesh")
    verbose && printstyled( "  file $filename written (Mesh)\n", color=:cyan)
end


"""
    Mesh(filename)

Constructs a `Mesh` object based on a file in VTK legacy format or JSON format.
"""
function Mesh(filename::String; verbose::Bool=true, silent::Bool=false, format::String="", reorder::Bool=false)

    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    verbosity>0 && printstyled("Mesh loading: filename $filename\n", bold=true, color=:cyan)

    if format==""
        basename, format = splitext(filename)
        format = format[2:end]
    end

    if format=="vtk"
        verbosity>0 && print("  Reading VTK legacy format...\n")
        mesh = read_vtk(filename)
    elseif format=="json"
        verbosity>0 && print("  Reading JSON format...\n")
        mesh = read_json(filename)
    elseif format=="tetgen"
        verbosity>0 && print("  Reading tetgen output files...\n")
        mesh = read_tetgen(filename)
    else
        error("Could not read $format format file")
    end

    # Reorder nodal numbering
    if reorder
        verbosity>0 && print("  reordering points...\r")
        reorder!(mesh)
    end

    if verbosity>1
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
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
