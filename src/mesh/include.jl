include("vtk.jl")
include("quadrature.jl")

include("shape.jl")
export ShapeType, ALL_ISO_SHAPES, ShapeFamily
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

include("point.jl")
export Point, Cell, hash, get_x, get_y, get_z
include("cell.jl")
export getcoords, get_point, get_points, get_faces, get_patches, cell_extent, cell_quality
export tag!, setquadrature!
include("partition.jl")

include("mesh.jl")
include("io.jl")
export Mesh, fixup!, quality!, reorder!, save, get_surface, get_neighbors, threshold, datafields

include("block.jl")
export Block, Block2D, Block3D, BlockTruss, BlockCoords, BlockCylinder

include("block_inset.jl")
export BlockInset

include("refine.jl")
export hrefine, prefine

include("operators.jl")
export move!, array, copy, mirror, rotate!, polar, rollaxes!, scale!

include("extrude.jl")
export extrude

include("smooth.jl")
export smooth!, laplacian_smooth!, fast_smooth!

include("split.jl")

include("embedded.jl")
export generate_embedded_cells!

# show function for Amaru types
#for datatype in (:ShapeType, :Point, :Cell, :Block, :Mesh, :UnstructuredGrid )
for datatype in (:ShapeType, :Point, :Cell, :Block, :Mesh)
    eval( quote
        function Base.show(io::IO, obj::$datatype)
            print_field_values(io, obj)
        end

        function Base.show(io::IO, array::Array{<:($datatype),1})
            print_array_values(io, array)
        end
    end )
end
