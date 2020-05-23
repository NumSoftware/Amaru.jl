include("vtk.jl")
include("quadrature.jl")

include("shape.jl")
export ShapeType, ALL_ISO_SHAPES, ShapeFamily
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

#include("point.jl")
export Node, Cell, hash, get_x, get_y, get_z
include("cell.jl")
include("collapse.jl")

export get_coords, get_node, get_nodes, get_faces, get_patches, cell_extent, cell_quality
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
export move!, array, copy, mirror, rotate!, polar, scale!, rollaxes!, changeaxes!

include("extrude.jl")
export extrude

include("revolve.jl")
export revolve

include("smooth.jl")
export smooth!, laplacian_smooth!, fast_smooth!

include("split.jl")

include("embedded.jl")
export generate_embedded_cells!

# show function for mesh related types
Base.show(io::IO, obj::ShapeType) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Block) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Mesh) = custom_dump(io, obj, 2, "")
