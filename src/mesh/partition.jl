# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh


mutable struct CellPartition
    cellpartition::Array{Array{Cell,1},3}
    bbox::Array{Float64,2}
    lbin::Float64
    function CellPartition(nx=0, ny=0, nz=0, bbox=nothing)
        this = new()
        this.cellpartition = Array{Array{Cell,1}}(undef, nx, nx, nx)
        this.bbox = zeros(0,0)
        # this.lbin = ..
        return this
    end
end


function get_bin_cells(cellpartition::CellPartition, p::Point) # returns an array of cells: TODO: untested function
    minx, miny, minz = cellpartition.bbox[1,:]
    lbin = cellpartition.lbin
    ix = floor(Int, (p.x - minx)/lbin) + 1
    iy = floor(Int, (p.y - miny)/lbin) + 1
    iz = floor(Int, (p.z - minz)/lbin) + 1
    return cellpartition.cellpartition[ix, iy, iz]
end

function add_cell(cellpartition::CellPartition, cell::Cell) # TODO: untested function
    minx, miny, minz = cellpartition.bbox[1,:]
    lbin = cellpartition.lbin
    bbox = bounding_box(cell)
    X, Y, Z  = bbox[:,1], bbox[:,2], bbox[:,3]
    verts    = [ [x y z] for z in Z, y in Y, x in X ]
    cell_pos = Set()
    for V in verts
        x, y, z = V
        ix = floor(Int, (x - minx)/lbin) + 1
        iy = floor(Int, (y - miny)/lbin) + 1
        iz = floor(Int, (z - minz)/lbin) + 1
        push!(cell_pos, (ix, iy, iz))
    end

    for (ix, iy, iz) in cell_pos
        push!(cellpartition.cellpartition[ix, iy, iz], cell)
    end
end

# Grepares a group of cellpartition that contain given cells
function build_bins(cells::Array{Cell,1}, cellpartition::CellPartition)
    # Get all points
    cellpartition.bbox = bounding_box(cells)
    minx, miny, minz = cellpartition.bbox[1,:]

    # Get max global lengths
    Lx, Ly, Lz = cellpartition.bbox[2,:] - cellpartition.bbox[1,:]
    max_L = max(Lx, Ly, Lz)

    # Get max cell lengths
    max_l = 0.0
    for cell in cells
        if cell.shape.family!=SOLID_SHAPE continue end
        bbox = bounding_box(cell)
        l    = maximum(bbox[2,:] - bbox[1,:])
        if l>max_l; max_l = l end
    end

    # Get number of divisions
    ndiv = min(50, 1*floor(Int, max_L/max_l)) # calibrate for cellpartition efficiency

    lbin = max_L/ndiv     # Get bin length
    cellpartition.lbin = lbin

    nx = floor(Int, Lx/lbin) + 1
    ny = floor(Int, Ly/lbin) + 1
    nz = floor(Int, Lz/lbin) + 1

    # Allocate cellpartition
    cellpartition.cellpartition = Array{Array{Cell,1}}(undef, nx, ny, nz)
    for k=1:nz, j=1:ny, i=1:nx
        cellpartition.cellpartition[i,j,k] = Cell[]
    end

    # Fill cellpartition
    for cell in cells
        bbox = bounding_box(cell)
        X, Y, Z  = bbox[:,1], bbox[:,2], bbox[:,3]
        verts    = [ (x, y, z) for z in Z, y in Y, x in X ]
        cell_loc = Set() # cells can be in more than one bin
        for (x, y, z) in verts
            ix = floor(Int, (x - minx)/lbin) + 1
            iy = floor(Int, (y - miny)/lbin) + 1
            iz = floor(Int, (z - minz)/lbin) + 1
            push!(cell_loc, (ix, iy, iz))
        end

        for (ix, iy, iz) in cell_loc
            push!(cellpartition.cellpartition[ix, iy, iz], cell)
        end
    end
end

# Find the cell that contains a given point
function find_cell(X::Array{Float64,1}, cells::Array{Cell,1}, cellpartition::CellPartition, tol::Float64, exc_cells::Array{Cell,1})
    # Point coordinates
    x, y, z = vcat(X, 0)[1:3]
    lbin = cellpartition.lbin

    # Build cellpartition if empty
    length(cellpartition.cellpartition) == 0 && build_bins(cells, cellpartition)

    for attempt=1:2
        Cmin = reshape(cellpartition.bbox[1,:],3)
        Cmax = reshape(cellpartition.bbox[2,:],3)
        lbin = cellpartition.lbin

        if any(X .< Cmin .- tol) || any(X .> Cmax .+ tol)
            error("find_cell: point outside bounding box")
        end

        # Find bin index
        ix = floor(Int, (x - Cmin[1])/lbin) + 1
        iy = floor(Int, (y - Cmin[2])/lbin) + 1
        iz = floor(Int, (z - Cmin[3])/lbin) + 1

        # Search cell in bin
        bin = cellpartition.cellpartition[ix, iy, iz]
        for cell in bin
            coords = getcoords(cell)
            if is_inside(cell.shape, coords, X, tol) && !(cell in exc_cells)
                return cell
            end
        end

        # If not found in the first attempt try then rebuild cellpartition
        if attempt==1
            @warn "find_cell: Bin search failed. Rebuilding cellpartition..."
            build_bins(cells, cellpartition)
        end
    end

    @warn "find_cell: Bin search failed"
    return nothing
end
