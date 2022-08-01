# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


mutable struct ElemPartition
    partition::Array{Array{AbstractCell,1},3}
    bbox::Array{Float64,2}
    lbin::Float64
    function ElemPartition(nx=0, ny=0, nz=0, bbox=nothing)
        this = new()
        this.partition = Array{Array{AbstractCell,1}}(undef, nx, nx, nx)
        this.bbox = zeros(0,0)
        # this.lbin = ..
        return this
    end
end


function get_bin_cells(cellpartition::ElemPartition, p::Node) # returns an array of cells: TODO: untested function
    minx, miny, minz = cellpartition.bbox[1,:]
    lbin = cellpartition.lbin
    ix = floor(Int, (p.coord.x - minx)/lbin) + 1
    iy = floor(Int, (p.coord.y - miny)/lbin) + 1
    iz = floor(Int, (p.coord.z - minz)/lbin) + 1
    return cellpartition.partition[ix, iy, iz]
end

function add_cell(cellpartition::ElemPartition, cell::AbstractCell) # TODO: untested function
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
        push!(cellpartition.partition[ix, iy, iz], cell)
    end
end

# Grepares a group of cellpartition that contain given cells
function build_bins(cells::Array{<:AbstractCell,1}, cellpartition::ElemPartition)
    # Get all points
    cellpartition.bbox = bounding_box(cells)
    minx, miny, minz = cellpartition.bbox[1,:]

    # Get max global lengths
    Lx, Ly, Lz = cellpartition.bbox[2,:] - cellpartition.bbox[1,:]
    max_L = max(Lx, Ly, Lz)

    # Get max cell lengths
    max_l = 0.0
    for cell in cells
        cell.shape.family==SOLID_CELL || continue
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
    cellpartition.partition = Array{Array{AbstractCell,1}}(undef, nx, ny, nz)
    for k in 1:nz, j=1:ny, i=1:nx
        cellpartition.partition[i,j,k] = AbstractCell[]
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
            push!(cellpartition.partition[ix, iy, iz], cell)
        end
    end
end

# Find the cell that contains a given point
function find_elem(X::Array{Float64,1}, cells::Array{<:AbstractCell,1}, cellpartition::ElemPartition, tol::Float64=1e-7; exclude::Array{<:AbstractCell,1}=AbstractCell[])
    # Node coordinates
    x, y, z = vcat(X, 0)[1:3]
    lbin = cellpartition.lbin

    # Build cellpartition if empty
    length(cellpartition.partition) == 0 && build_bins(cells, cellpartition)

    for attempt=1:2
        Cmin = reshape(cellpartition.bbox[1,:],3)
        Cmax = reshape(cellpartition.bbox[2,:],3)
        lbin = cellpartition.lbin

        if any(X .< Cmin .- tol) || any(X .> Cmax .+ tol)
            error("find_elem: point outside bounding box $X")
        end

        # Find bin index
        ix = floor(Int, (x - Cmin[1])/lbin) + 1
        iy = floor(Int, (y - Cmin[2])/lbin) + 1
        iz = floor(Int, (z - Cmin[3])/lbin) + 1

        # Search cell in bin
        bin = cellpartition.partition[ix, iy, iz]
        for cell in bin
            coords =getcoords(cell)
            if is_inside(cell.shape, coords, X, tol) && !(cell in exclude)
                return cell
            end
        end

        # If not found in the first attempt try then rebuild cellpartition
        if attempt==1
            notify("find_elem: Bin search failed. Rebuilding cellpartition...")
            build_bins(cells, cellpartition)
        end
    end

    alert("find_elem: Bin search failed")
    return nothing
end
