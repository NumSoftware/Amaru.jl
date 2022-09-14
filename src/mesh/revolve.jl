
"""
    revolve(mesh; base=[0,0,0], axis=[0,0,1], minangle=0, maxangle=360, angle=360, n=8, collapse=true, lagrangian=false)

Generates a 3D mesh by revolving a 1D o3 2D `mesh` using a `base` point, 
an `axis`, a rotation `angle` and a number of divisions `n`.
If `collapse` is true, each element with concident nodes is collapsed
into simpler element with fewer nodes. If lagrangian is true, lagrangian elements
are generated when posible.
"""
function revolve(mesh::Mesh; 
                 base::AbstractArray{<:Real,1} = [0,0,0],
                 axis::AbstractArray{<:Real,1} = [0,0,1],
                 minangle  ::Real = 0,
                 maxangle  ::Real = 360,
                 angle     ::Real = NaN,
                 n         ::Int  = 8,
                 collapse  ::Bool = true,
                 lagrangian::Bool = false
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
    minθ = minangle*pi/180
    maxθ = maxangle*pi/180
    Δθ = (maxθ-minθ)/n

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

    for θ in range(minθ,step=Δθ,length=n)
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
                for (n, R) in zip([n1, n2, n2, n1], [Rend, Rend, Rini, Rini])
                    coord = base + R*(n.coord-base)*conj(R)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape==LIN3
                n1, n2, n3 = cell.nodes
                newnodes  = [n1, n2, n2, n1, n3, n2, n3, n1]
                angles = [Rend, Rend, Rini, Rini, Rend, Rhalf, Rini, Rhalf]

                if lagrangian
                    newshape = QUAD9
                    push!(newnodes, n3)
                    push!(angles, Rhalf)
                else
                    newshape = QUAD8
                end

                for (n, R) in zip(newnodes, angles)
                    coord = base + R*(n.coord-base)*conj(R)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape==LIN4
                newshape = QUAD12
                n1, n2, n3, n4 = cell.nodes
                θ13 = θ + 1/3*Δθ
                θ23 = θ + 2/3*Δθ
                R13 = Quaternion(cos(θ13/2), axis[1]*sin(θ13/2), axis[2]*sin(θ13/2), axis[3]*sin(θ13/2))
                R23 = Quaternion(cos(θ23/2), axis[1]*sin(θ23/2), axis[2]*sin(θ23/2), axis[3]*sin(θ23/2))

                for (n, R) in zip([n1, n2, n2, n1, n3, n4, n2, n2, n4, n3, n1, n1], [Rend, Rend, Rini, Rini, Rend, Rend, R23, R13, Rini, Rini, R13, R23])
                    coord = base + R*(n.coord-base)*conj(R)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape in (TRI3, QUAD4)
                newshape = cell.shape==TRI3 ? WED6 : HEX8
                nnodes = length(cell.nodes)

                for R in (Rend, Rini)
                    for n in cell.nodes[1:nnodes]
                        coord = base + R*(n.coord-base)*conj(R)
                        push!(nodes, Node(coord))
                    end
                end
            elseif cell.shape==TRI6
                newshape = WED15

                for R in (Rend, Rini)
                    for n in cell.nodes[1:3]
                        coord = base + R*(n.coord-base)*conj(R)
                        push!(nodes, Node(coord))
                    end
                end

                for R in (Rend, Rini)
                    for n in cell.nodes[4:6]
                        coord = base + R*(n.coord-base)*conj(R)
                        push!(nodes, Node(coord))
                    end
                end

                for n in cell.nodes[1:3]
                    coord = base + Rhalf*(n.coord-base)*conj(Rhalf)
                    push!(nodes, Node(coord))
                end
            elseif cell.shape==QUAD8
                newshape = HEX20

                for R in (Rend, Rini)
                    for n in cell.nodes[1:4]
                        coord = base + R*(n.coord-base)*conj(R)
                        push!(nodes, Node(coord))
                    end
                end

                for R in (Rend, Rini)
                    for n in cell.nodes[5:8]
                        coord = base + R*(n.coord-base)*conj(R)
                        push!(nodes, Node(coord))
                    end
                end

                for n in cell.nodes[1:4]
                    coord = base + Rhalf*(n.coord-base)*conj(Rhalf)
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

    for elem in cells
        isinverted(elem) && flip!(elem)
    end

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    fixup!(newmesh, reorder=true)

    return newmesh
end

"""
    revolve(node; base=[0,0,0], axis=[0,0,1], minangle=0, maxangle=360, angle=360, n=8)

Generates a mesh composed by LIN3 elements by revolving a `node` using a `base` point, 
an `axis`, a rotation `angle` and a number of divisions `n`.
"""
function revolve(node::Node; 
    base    ::AbstractArray{<:Real,1} = [0,0,0],
    axis    ::AbstractArray{<:Real,1} = [0,0,1],
    minangle::Real = 0,
    maxangle::Real = 360,
    angle   ::Real = NaN,
    n       ::Int  = 8,
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
    minθ = minangle*pi/180
    maxθ = maxangle*pi/180
    Δθ = (maxθ-minθ)/n

    for θ in range(minθ,step=Δθ,length=n)
        θhalf = θ + Δθ/2
        θend  = θ + Δθ

        Rini  = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
        Rhalf = Quaternion(cos(θhalf/2), axis[1]*sin(θhalf/2), axis[2]*sin(θhalf/2), axis[3]*sin(θhalf/2))
        Rend  = Quaternion(cos(θend/2), axis[1]*sin(θend/2), axis[2]*sin(θend/2), axis[3]*sin(θend/2))

        nodes = Node[]
        for R in (Rend, Rini, Rhalf)
            coord = base + R*(node.coord-base)*conj(R)
            push!(nodes, Node(coord))
        end

        newcell = Cell(LIN3, nodes)
        push!(cells, newcell)
    end

    # Merge points
    point_dict = Dict{UInt64,Node}( hash(n) => n for c in cells for n in c.nodes )

    for cell in cells
        cell.nodes = Node[ point_dict[hash(n)] for n in cell.nodes ]
    end
    nodes = collect(values(point_dict))

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    fixup!(newmesh, reorder=true)

    return newmesh

end