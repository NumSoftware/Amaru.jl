
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
    edge::Edge
    faces::Array{Face,1}
    function EdgeFaces()
    end
end

export get_outline_edges, get_outline

function get_outline_edges(cells::Array{<:AbstractCell,1}; angle=150)
    faces_dict = Dict{UInt64, Cell}()

    # Get faces
    for cell in cells
        cell.shape.family == SOLID_SHAPE || continue # only bulk cells
        if cell.shape.ndim==2
            hs = hash(cell)
            faces_dict[hs] = cell
            continue
        end

        for face in getfaces(cell)
            hs = hash(face)
            if haskey(faces_dict, hs)
                delete!(faces_dict, hs)
            else
                faces_dict[hs] = face
            end
        end
    end

    faces = values(faces_dict)

    # Get normals
    normals = Dict{Face,Array{Float64,1}}( f => get_facet_normal(f) for f in faces ) 
    edge_dict = Dict{UInt64,Cell}()
    outline = Edge[]

    # Get edges with non-coplanar adjacent faces
    for face in faces
        n1 = normals[face]
        for edge in getedges(face)
            hs = hash(edge)
            edge0 = get(edge_dict, hs, nothing)
            if edge0===nothing
                edge_dict[hs] = edge
            else
                delete!(edge_dict, hs)
                n2 = normals[edge0.owner]
                α = 180 - acos( abs(clamp(dot(n1,n2),-1,1)) )*180/pi
                α = round(α, digits=2)
                α<=angle && push!(outline, edge)
            end
        end
    end

    return outline
end


function get_outline(mesh::Mesh; angle=150)
    if mesh.ndim==2
        return Mesh(mesh.edges)
    else
        return Mesh(get_outline_edges(mesh.faces; angle=angle))
    end

    return Mesh(outline)
end
