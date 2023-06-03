"""
    getpolygon(nodes)

Orders an array of `nodes` to get a closed polygon.
"""
function getpolygon(nodes::Array{Node,1})
    # the nodes might be in placed an arbitrary 3D plane

    # find the best fit plane
    C = getcoords(nodes, 3)
    nnodes = length(nodes)

    # move the coordinates to avoid singular case
    # when the regression line/planes crosses the origin
    Cm =  C .+ [pi 2*pi 3*pi]

    I = ones(nnodes)
    N = pinv(Cm)*I # best fit normal
    N = round.(N, digits=15) # to avoid almost zero values

    # rotation matrix
    if N[1]==0.0
        V1 = [ 1.0, 0.0, 0.0 ]
    elseif N[2]==0.0
        V1 = [ 0.0, 1.0, 0.0 ]
    else
        V1 = [ 0.0, 0.0, 1.0 ]
    end

    V2 = cross(N, V1)
    V3 = cross(V1, V2)

    normalize!(V2)
    normalize!(V3)

    R = [ V1 V2 V3 ]
    C = C*R

    # get center
    Xc = vec(sum(C, dims=1))./nnodes

    # get polar angles
    angles = zeros(nnodes)
    for i in 1:nnodes
        x, y, _ = C[i,:]-Xc
        angles[i] = atan(y,x)
    end

    perm = sortperm(angles)
    C = C[perm,:]

    # get internal angles
    minα = 2*pi
    mini = 0
    for i in 1:nnodes
        prev = i==1 ? nnodes : i-1
        next = i==nnodes ? 1 : i+1
        X = C[i,:]
        V1 = normalize(C[prev,:]-X)
        V2 = normalize(C[next,:]-X)
        α = acos(dot(V1,V2))
        if α < minα
            minα=α
            mini = i
        end
    end

    # put point with smaller internal angle as first
    circshift(perm, 1-mini)

    # Amaru.@showm getcoords(nodes[perm])*R

    # build contour
    return nodes[perm]
end


function makecells(nodes::Array{Node,1}; ndim=2, quadratic=false)
    cells = Cell[]
    nnodes = length(nodes)
    # @assert nnodes==shape.npoints

    if ndim==2
        nodes = getpolygon(nodes)
        if nnodes==3
            cell = Cell(TRI3, nodes)
            push!(cells, cell)
        elseif nnodes==4
            cell = Cell(QUAD4, nodes)
            push!(cells, cell)
        elseif nnodes==6
            cell1 = Cell(QUAD4, nodes[1:4])
            cell2 = Cell(QUAD4, nodes[[4,5,6,1]])
            push!(cells, cell1)
            push!(cells, cell2)
        elseif nnodes==5
            cell1 = Cell(TRI3, nodes[1:3])
            cell2 = Cell(QUAD4, nodes[[3,4,5,1]])
            push!(cells, cell1)
            push!(cells, cell2)
        end

        return cells
    end
    error("Not implemented")

end


"""
    slice(mesh; base, axis)

Generates a planar mesh by slicing a 3D `mesh` using a plane defined by a `base` point and an `axis`.
The original nodal data is interpolated to the nodes of the resulting mesh.
"""
function slice(
    mesh::AbstractDomain; 
    base::AbstractArray{<:Real,1} = Float64[],
    axis::AbstractArray{<:Real,1} = Float64[],
    )

    @check length(base)==3
    @check length(axis)==3
    @check mesh.ndim==3

    axis = Vec3(axis)
    base = Vec3(base)

    # get all mesh edges
    edgedict = Dict{UInt, Edge}()
    for cell in mesh.elems
        for edge in getedges(cell)
            hs = hash(edge)
            edgedict[hs] = edge
        end
    end

    edges = collect(values(edgedict))

    # find intersection points between edges and the plane
    edgeint = Dict{UInt, Node}()
    tolξ = 1e-9
    for edge in edges
        
        Xa = edge.nodes[1].coord
        Xb = edge.nodes[2].coord
        
        fa = dot(Xa-base, axis)
        fb = dot(Xb-base, axis)

        if fa==0.0
            edgeint[hash(edge)] = Node(Xa)
            continue
        end
        if fb==0.0
            edgeint[hash(edge)] = Node(Xb)
            continue
        end
        fa*fb > 0.0 && continue

        Xi = (Xa+Xb)/2
        n  = floor(Int, log(2, 2/tolξ)) + 1
        Xi = Vec3(0,0,0)
        C  = getcoords(edge)
        a  = -1.0
        b  = +1.0
        for i in 1:n
            ξ = (a+b)/2
            N = edge.shape.func([ξ])
            Xi = C'*N

            fi = dot(Xi-base, axis)
            if fa*fi < 0.0
                fb = fi
                b = ξ
            else
                fa = fi
                a = ξ
            end
            fi == 0 && break
        end

        Xi = round.(Xi, digits=8)
        edgeint[hash(edge)] = Node(Xi)
    end

    # mount slice elements
    newcells = Cell[]
    for cell in mesh.elems
        nodedict = Dict{UInt,Node}()
        for edge in cell.edges
            Xi = get(edgeint, hash(edge), nothing)
            Xi === nothing && continue
            nodedict[hash(Xi)] = Xi
        end
        
        nodes = collect(values(nodedict))
        nnodes = length(nodes) 
        nnodes == 0 && continue

        cells = makecells(nodes, ndim=2)
        for c in cells
            c.tag = cell.tag
            c.owner = cell
        end
        append!(newcells, cells)

    end

    # generate new mesh
    newmesh = Mesh()
    newmesh.nodes = getnodes(newcells)
    newmesh.elems = newcells
    fixup!(newmesh)

    # interpolate node data
    nnodes = length(newmesh.nodes)
    for (key,data) in mesh.node_data
        haskey(newmesh.node_data, key) && continue
        count = zeros(Int, nnodes)
        sz    = size(data)
        dim   = length(sz)
        if dim==1
            newdata = zeros(eltype(data), nnodes)
        else
            newdata = zeros(eltype(data), nnodes, sz[2])
        end

        for cell in newmesh.elems
            ocell = cell.owner
            ids = [ node.id for node in ocell.nodes ]
            V = data[ids,:]
            coords = getcoords(ocell)
            for node in cell.nodes
                Ξ = inverse_map(ocell.shape, coords, node.coord)
                N = ocell.shape.func(Ξ)
                val = V'*N
                newdata[node.id, :] .+= val
                count[node.id] += 1
            end
        end

        newmesh.node_data[key] = newdata./count
    end
    
    # update cell data
    ncelss = length(newmesh.elems)
    for (key,data) in mesh.elem_data
        haskey(newmesh.elem_data, key) && continue
        sz = size(data)
        dim = length(sz)
        if dim==1
            newdata = zeros(eltype(data), ncelss)
        else
            newdata = zeros(eltype(data), ncelss, size[2])
        end

        for cell in newmesh.elems
            ocell = cell.owner
            val = data[ocell.id, :]
            newdata[cell.id, :] .= val 
        end

        newmesh.elem_data[key] = newdata
    end

    return newmesh
end

