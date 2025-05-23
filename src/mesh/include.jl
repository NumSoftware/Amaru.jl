

include("block.jl")
export Block, BlockTruss, BlockCoords, BlockCylinder, BlockGrid

include("mesh-env.jl")

#include("point.jl")
export Node, Cell, hash, get_x, get_y, get_z
include("cell.jl")
include("collapse.jl")

export get_coords, get_node, getnodes, getfacets, getfaces, getedges, get_patches, cell_extent, cell_quality, cell_aspect_ratio
include("partition.jl")

include("mesh.jl")
include("structured.jl")
include("unstructured.jl")
include("gen_insets.jl")
include("genmesh.jl")

include("io.jl")
export Mesh, fixup!, quality!, sortnodes!, save, get_outer_facets, get_neighbors, threshold, datafields

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
export insert_cohesive_elements!, cracksmesh

include("embedded.jl")
export generate_embedded_cells!

include("outline.jl")

# show function for mesh related types
Base.show(io::IO, obj::CellShape) = _show(io, obj, 2)
Base.show(io::IO, obj::Block) = _show(io, obj, 2)
Base.show(io::IO, obj::Mesh) = _show(io, obj, 2)
