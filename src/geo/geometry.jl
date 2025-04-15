

GeoModel_params = [
    FunInfo(:GeoModel, "A geometry model"),
    KwArgInfo(:size, "Characteristic length for meshing", 0.1),
    KwArgInfo(:quiet, "Flag for quiet run", false),
]
@doc docstring(GeoModel_params) GeoModel


mutable struct GeoModel
    points::Vector{Point}
    edges::Vector{Edge}
    loops::Vector{Loop}
    faces::Vector{Face}
    volumes::Vector{Volume}
    subpaths::Vector{SubPath}
    blocks::Vector{AbstractBlock}
    size::Float64
    ndim::Int
    quiet::Bool
    _id::Int
    _volume_detection::Bool

    function GeoModel(; kwargs...)
        args  = checkargs(kwargs, Mesh_Geo_params)
        quiet = args.quiet
        size  = args.size
        quiet = args.quiet
        !quiet && printstyled("Geometry model:\n", bold=true, color=:cyan)
        return new( [], [], [], [], [], [], [], size, 0, quiet, 0, true )
    end
end


# Show functions
Base.show(io::IO, geo::GeoModel) = _show(io, geo, 3, "")




function Base.copy(geo::GeoModel)
    _geo = GeoModel(quiet=true)
    _geo.size = geo.size
    _geo.ndim = geo.ndim
    _geo._id = geo._id

    # copy points
    for point in geo.points
        _point = Point(point.coord, id=point.id, size=point.size, tag=point.tag)
        push!(_geo.points, _point)
    end

    # copy edges
    for edge in geo.edges
        ids = [ p.id for p in edge.points ]
        _points = [ _geo[id] for id in ids ]
        if edge isa Line
            _edge = Line(_points..., n=edge.n, tag=edge.tag, id=edge.id)
        elseif edge isa Arc
            _edge = Arc(_points..., n=edge.n, tag=edge.tag, id=edge.id)

            # extra points
            for p in edge.extrapoints
                push!(_edge.extrapoints, copy(p))
            end
        end

        push!(_geo.edges, _edge)
    end

    # update point edges
    for point in geo.points
        _point = _geo[point.id]
        _edge_ids = [ edge.id for edge in point.edges ]
        _edges = Edge[ _geo[id] for id in _edge_ids ]
        _point.edges = _edges
    end

    # copy Loops
    for loop in geo.loops
        _edge_ids = [ dart.edge.id for dart in loop.darts ]
        _edges = [ _geo[id] for id in _edge_ids ]

        _darts = [ Dart(_edge, forward=dart.forward) for (dart, _edge) in zip(loop.darts, _edges) ]

        _loop = Loop(_darts, id=loop.id, flat=loop.flat)
        push!(_geo.loops, _loop)
    end

    # copy faces
    for face in geo.faces
        _loop_ids = [ loop.id for loop in face.loops ]
        _loops = [ _geo[id] for id in _loop_ids ]

        _face = Face(_loops..., id=face.id, tag=face.tag)
        _face.flat = face.flat

        push!(_geo.faces, _face)
    end

    # update edge faces
    for edge in geo.edges
        _edge = _geo[edge.id]
        _face_ids = [ face.id for face in edge.faces ]
        _faces = [ _geo[id] for id in _face_ids ]
        _edge.faces = _faces
    end

    # copy volumes
    for volume in geo.volumes
        _face_ids = [ spin.face.id for spin in volume.spins ]
        _normals =  [ copy(spin.normal) for spin in volume.spins ]
        _spins = [ FaceSpin(_geo[id], _normal) for (id, _normal) in zip(_face_ids, _normals) ]
        _volume = Volume(_spins, id=volume.id, tag=volume.tag )
        push!(_geo.volumes, _volume)
    end

    # update volumes in faces
    for face in geo.faces
        _face = _geo[face.id]
        _volume_ids = [ volume.id for volume in face.volumes ]
        _volumes = [ geo[id] for id in _volume_ids ]
        _face.volumes = _volumes
    end

    return _geo

end


function Base.getindex(geo::GeoModel, idx::Int)
    for p in geo.points
        p.id==idx && return p
    end
    for l in geo.edges
        l.id==idx && return l
    end
    for lo in geo.loops
        lo.id==idx && return lo
    end
    for s in geo.faces
        s.id==idx && return s
    end
    for v in geo.volumes
        v.id==idx && return v
    end
    for sp in geo.subpaths
        sp.id==idx && return sp
    end
    for b in geo.blocks
        b.id==idx && return b
    end

    error("GeoModel: No entity with id $idx")
end


function getpoint(geo::GeoModel, p::Point)
    for pp in geo.points
        p==pp && return pp
    end

    return nothing
end


function getline(geo::GeoModel, line::Edge)
    idx = findfirst(==(line), geo.edges)
    idx === nothing && (idx=0) 
    return get(geo.edges, idx, nothing)
end


function getline(geo::GeoModel, p1::Point, p2::Point)
    l = Line(p1, p2)
    return getline(geo, l)
end


function getloop(geo::GeoModel, loop::Loop)
    idx = findfirst(==(loop), geo.loops)
    idx === nothing && return idx
    return geo.loops[idx]
end


function getface(geo::GeoModel, face::Face)
    idx = findfirst(==(face), geo.faces)
    idx === nothing && return idx
    return geo.faces[idx]
end


function getface(geo::GeoModel, loop::Loop)
    loops = [ face.loops[1] for face in geo.faces ]
    idx = findfirst(==(loop), loops)
    idx === nothing && return idx
    return geo.faces[idx]
end


function getvolume(geo::GeoModel, vol::Volume)
    idx = findfirst(==(vol), geo.volumes)
    idx === nothing && return idx
    return geo.volumes[idx]
end