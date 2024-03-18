# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

Mesh_Geo_params = [
    FunInfo(:Mesh, "Creates a `Mesh` structure from a `GeoModel` object."),
    ArgInfo(:geo, "GeoModel object"),
    KwArgInfo(:size, "Characteristic length for meshing", 0.1),
    KwArgInfo(:quadratic, "Flag for quadratic cells (gmsh)", false),
    KwArgInfo(:algorithm, "Algorithm for triangulation (gmsh)", :delaunay),
    KwArgInfo(:recombine, "Flag for cell recombination (gmsh)", false),
    KwArgInfo(:quiet, "Flag for quiet run", false),
]
@doc docstring(Mesh_Geo_params) Mesh(geo::GeoModel; kwargs...)


function Mesh(geo::GeoModel; kwargs...)
    args = checkargs([geo], kwargs, Mesh_Geo_params)
    quiet = args.quiet

    # check for blocks
    blocks  = [ b for b in geo.blocks if b isa Block ]
    iblocks = [ b for b in geo.blocks if b isa BlockInset ]

    if length(geo.surfaces)>0 || length(geo.volumes)>0
        !quiet && printstyled("Unstructured mesh generation:\n", bold=true, color=:cyan)

        length(blocks)>0 && warn("Mesh: Blocks are being ignored")
        mesh = mesh_unstructured(geo; args...)
        if length(iblocks)>0
            mesh = Mesh(mesh, iblocks; args...)
        end
    elseif length(blocks)>0
        !quiet && printstyled("Structured mesh generation:\n", bold=true, color=:cyan)
        mesh = mesh_structured(geo; args...)
    else
        error("Mesh: No blocks or surfaces/volumes found")
    end

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        @printf "  %4dd mesh\033[K\n" mesh.env.ndim
        @printf "  %5d nodes\n" npoints
        @printf "  %5d cells\n" ncells
        nfaces  = length(mesh.faces)
        @printf "  %5d faces\n" nfaces
        nedges  = length(mesh.edges)
        nedges>0 && @printf "  %5d surface edges\n" nedges
    end

    return mesh
end


function Mesh(
    items     ::Union{Mesh, AbstractBlock, Array{<:Union{AbstractBlock, Array},1}}...;
    ndim      ::Int = 0,
    quiet     ::Bool = false,
)

    # Flatten items list
    fitems = flatten(items)

    # Get list of blocks and update mesh objects
    blocks = AbstractBlock[]
    meshes = Mesh[]
    blocks = filter(b->b isa AbstractBlock, fitems)
    meshes = filter(b->b isa Mesh, fitems)

    # Find ndim
    # ndim = max(ndim, maximum((b.ndim for b in blocks), init=0), maximum((m.env.ndim for m in meshes), init=0) )
    
    nmeshes = length(meshes)
    nblocks = length(blocks)
    if !quiet
        printstyled("Mesh generation:\n", bold=true, color=:cyan)
        nmeshes>0 && @printf "  %5d meshes\n" nmeshes
        @printf "  %5d blocks\n" nblocks
    end

    # Join meshes
    mesh = Mesh(ndim)
    for m in meshes
        join_mesh!(mesh, m)
    end

    # Split blocks: generates nodes and cells
    for (i,b) in enumerate(blocks)
        split_block!(mesh, b)
        quiet || print("  spliting block ", i, "...    \r")
    end

    # Updates numbering, quality, facets and edges
    syncronize!(mesh)

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        @printf "  %4dd mesh\033[K\n" mesh.env.ndim
        @printf "  %5d nodes\n" npoints
        @printf "  %5d cells\n" ncells
        nfaces  = length(mesh.faces)
        @printf "  %5d faces\n" nfaces
        nedges  = length(mesh.edges)
        nedges>0 && @printf "  %5d surface edges\n" nedges
    end

    return mesh
end