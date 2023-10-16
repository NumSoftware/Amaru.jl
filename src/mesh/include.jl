include("geo-model.jl")

include("vtk.jl")
include("quadrature.jl")

include("shape.jl")
export CellShape, ALL_ISO_SHAPES, CellFamily
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

include("mesh-env.jl")

#include("point.jl")
export Node, Cell, hash, get_x, get_y, get_z
include("cell.jl")
include("collapse.jl")

export get_coords, get_node, getnodes, getfaces, get_patches, cell_extent, cell_quality, cell_aspect_ratio
include("partition.jl")

include("mesh.jl")
include("io.jl")
export Mesh, fixup!, quality!, reorder!, save, get_surface, get_neighbors, threshold, datafields

include("block.jl")
export Block, BlockTruss, BlockCoords, BlockCylinder, BlockGrid

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

include("convex_hull.jl")

include("slice.jl")
export slice

include("smooth.jl")
export smooth!, laplacian_smooth!, fast_smooth!

include("split.jl")
export generate_joints!, cracksmesh

include("embedded.jl")
export generate_embedded_cells!

include("outline.jl")

# show function for mesh related types
Base.show(io::IO, obj::CellShape) = _show(io, obj, 2)
Base.show(io::IO, obj::Block) = _show(io, obj, 2)
Base.show(io::IO, obj::Mesh) = _show(io, obj, 2)
