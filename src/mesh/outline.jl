
function get_facet_normal(face::AbstractCell)
    ndim = 1 + face.shape.ndim
    C = getcoords(face, ndim)

    if ndim==2
        C .+= [pi pi^1.1]
    else
        C .+= [pi pi^1.1 pi^1.2]
    end

    # calculate the normal
    I = ones(size(C,1))
    N = pinv(C)*I # best fit normal
    normalize!(N) # get unitary vector

    return N
end


mutable struct EdgeFaces
    edge::CellEdge
    faces::Array{CellFace,1}
    function EdgeFaces()
    end
end

export get_feature_edges, get_feature_mesh


function get_feature_edges(cells::Vector{<:AbstractCell}; angle=150)
    
    faces_dict = Dict{UInt64, Cell}()

    # Get faces
    for cell in cells
        cell.shape.family == BULKCELL || continue # only bulk cells
        if cell.shape.ndim==2
            hs = hash(cell)
            faces_dict[hs] = cell
            continue
        end

        for face in getedges(cell)
            hs = hash(face)
            if haskey(faces_dict, hs)
                delete!(faces_dict, hs)
            else
                faces_dict[hs] = face
            end
        end
    end

    faces = values(faces_dict)
    # @show faces

    # Get normals
    normals = Dict{CellFace,Array{Float64,1}}( f => get_facet_normal(f) for f in faces ) 
    face_edge_d = Dict{UInt64,Cell}()
    outline = CellEdge[]

    # Get edges with non-coplanar adjacent faces
    for face in faces
        n1 = normals[face]
        for edge in getedges(face)
            hs = hash(edge)
            edge0 = get(face_edge_d, hs, nothing)
            if edge0===nothing
                face_edge_d[hs] = edge
            else
                delete!(face_edge_d, hs)
                n2 = normals[edge0.owner]
                α = 180 - acos( abs(clamp(dot(n1,n2),-1,1)) )*180/pi
                α = round(α, digits=2)
                α<=angle && push!(outline, edge)
            end
        end
    end
    # @show length(face_edge_d)
    outline = [ outline; collect(values(face_edge_d)) ]

    return outline
end


function get_feature_mesh(mesh::Mesh; angle=150)
    if mesh.ndim==2
        return Mesh(mesh.edges)
    else
        return Mesh(get_feature_edges(mesh.faces; angle=angle))
    end
end
