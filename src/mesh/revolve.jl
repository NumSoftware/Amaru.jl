
"""
    revolve(mesh; base=[0,0,0], axis=[0,0,1], minangle=0, maxangle=360, angle=360, n=8, collapse=true)

Generates a 3D mesh by revolving a 2D `mesh` using a `base` point, 
a direction `axis`, a rotation `angle` and a number of divisions `n`.
If `collapse` is true, elements with concident nodes are collapsed
into simpler elements with fewer nodes.
"""
function revolve(mesh::Mesh; 
                 base::AbstractArray{<:Real,1} = [0.0,0.0,0.0],
                 axis::AbstractArray{<:Real,1} = [0.0,1.0,0.0],
                 minangle::Real = 0,
                 maxangle::Real = 360,
                 angle::Real = NaN,
                 n::Int=8,
                 collapse::Bool=true
                )
    @assert length(axis)==3
    @assert length(base)==3
    @assert n>0

    axis = Vec3(normalize(axis))
    base = Vec3(base)

    if !isnan(angle)
        minangle = 0
        maxangle = angle
    end

    @assert 0<=minangle<360
    @assert 0<maxangle<=360
    @assert minangle<maxangle

    cells = Cell[]
    Δθ = (maxangle-minangle)/n*pi/180

    if !collapse
        # move points that lie at the axis
        normal = Vec3(axis[2], -axis[1], 0.0) # TODO: improve this strategy
        for node in mesh.nodes
            X = node.coord
            if norm( dot(X-base, axis)*axis - X ) < 1e-8
                node.coord = node.coord + normal*1e-6
            end
        end
    end

    for θ in range(minangle,step=Δθ,length=n)
        θhalf = θ + Δθ/2
        θend  = θ + Δθ

        Rini  = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
        Rhalf = Quaternion(cos(θhalf/2), axis[1]*sin(θhalf/2), axis[2]*sin(θhalf/2), axis[3]*sin(θhalf/2))
        Rend  = Quaternion(cos(θend/2), axis[1]*sin(θend/2), axis[2]*sin(θend/2), axis[3]*sin(θend/2))

        for cell in mesh.elems
            nodes = Node[]

            if cell.shape==LIN2
                newshape = QUAD4
                n1, n2 = cell.nodes
                nnodes = length(cell.nodes)

                for (n, R) in zip([n1, n2, n2, n1], [Rend, Rend, Rini, Rini])
                    coord = base + R*(n.coord-base)*conj(R)
                    coord = round.(coord, digits=8)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape==LIN3
                newshape = QUAD8
                n1, n2, n3 = cell.nodes
                nnodes = length(cell.nodes)

                for (n, R) in zip([n1, n2, n2, n1, n3, n2, n3, n1], [Rend, Rend, Rini, Rini, Rend, Rhalf, Rini, Rhalf])
                    coord = base + R*(n.coord-base)*conj(R)
                    coord = round.(coord, digits=8)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape in (TRI3, QUAD4)
                newshape = cell.shape==TRI3 ? WED6 : HEX8
                nnodes = length(cell.nodes)

                for R in (Rend, Rini)
                    for n in cell.nodes[1:nnodes]
                        coord = base + R*(n.coord-base)*conj(R)
                        coord = round.(coord, digits=8)
                        push!(nodes, Node(coord))
                    end
                end
            elseif cell.shape==TRI6
                newshape = WED15

                for R in (Rend, Rini)
                    for n in cell.nodes[1:3]
                        coord = base + R*(n.coord-base)*conj(R)
                        coord = round.(coord, digits=8)
                        push!(nodes, Node(coord))
                    end
                end

                for R in (Rend, Rini)
                    for n in cell.nodes[4:6]
                        coord = base + R*(n.coord-base)*conj(R)
                        coord = round.(coord, digits=8)
                        push!(nodes, Node(coord))
                    end
                end

                for n in cell.nodes[1:3]
                    coord = base + Rhalf*(n.coord-base)*conj(Rhalf)
                        coord = round.(coord, digits=8)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape==QUAD8
                newshape = HEX20

                for R in (Rend, Rini)
                    for n in cell.nodes[1:4]
                        coord = base + R*(n.coord-base)*conj(R)
                        coord = round.(coord, digits=8)
                        push!(nodes, Node(coord))
                    end
                end

                for R in (Rend, Rini)
                    for n in cell.nodes[5:8]
                        coord = base + R*(n.coord-base)*conj(R)
                        coord = round.(coord, digits=8)
                        push!(nodes, Node(coord))
                    end
                end

                for n in cell.nodes[1:4]
                    coord = base + Rhalf*(n.coord-base)*conj(Rhalf)
                        coord = round.(coord, digits=8)
                    push!(nodes, Node(coord))
                end
            else
                error("revolve: Cell shape $(cell.shape.name) is not supported")
            end

            newcell = Cell(newshape, nodes, tag=cell.tag)
            push!(cells, newcell)
        end
    end

    # Merge points
    point_dict = Dict{UInt64,Node}( hash(n) => n for c in cells for n in c.nodes )

    for cell in cells
        cell.nodes = Node[ point_dict[hash(n)] for n in cell.nodes ]
    end
    nodes = collect(values(point_dict))


    if collapse
        for elem in cells
            inaxis = false
            for node in elem.nodes
                X = node.coord
                if norm( dot(X-base, axis)*axis - X ) < 1e-8
                    inaxis = true
                    break
                end
            end
            inaxis && Amaru.collapse!(elem)
        end
    else
        # move back points that lied at the axis
        for node in mesh.nodes
            X = node.coord
            if norm( dot(X-base, axis)*axis - X ) < 1.00001e-6
                node.coord = round.( dot(X-base, axis)*axis , digits=8)
            end
        end
    end

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    fixup!(newmesh, reorder=true)

    return newmesh
end

