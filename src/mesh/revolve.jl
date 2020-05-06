
"""
    revolve(mesh; base=[0,0,0], axis=[0,0,1], minangle=0, maxangle=360, angle=360, n=8)

Generates a 3D mesh by revolving a 2D `mesh` using a `base` point, 
a direction `axis`, a rotation `angle` a number `n` of divisions.
"""
function revolve(mesh::Mesh; 
                 base::AbstractArray{<:Real,1} = [0.0,0.0,0.0],
                 axis::AbstractArray{<:Real,1} = [0.0,1.0,0.0],
                 minangle::Real = 0,
                 maxangle::Real = 360,
                 angle::Real = NaN,
                 n::Int=4
                )
    @assert length(axis)==3
    @assert length(base)==3
    #@assert n>1

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

    # move points that lie at the axis
    #normal = Vec3(axis[2], -axis[1], 0.0) # TODO: improve
    #for node in mesh.nodes
        #X = node.coord
        #if norm( dot(X-base, axis)*axis - X ) < 1e-8
            #node.coord = node.coord + normal*1e-6
        #end
    #end

    for θ in range(minangle,step=Δθ,length=n)
        θhalf = θ + Δθ/2
        θend  = θ + Δθ

        Rini  = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
        Rhalf = Quaternion(cos(θhalf/2), axis[1]*sin(θhalf/2), axis[2]*sin(θhalf/2), axis[3]*sin(θhalf/2))
        Rend  = Quaternion(cos(θend/2), axis[1]*sin(θend/2), axis[2]*sin(θend/2), axis[3]*sin(θend/2))

        for cell in mesh.elems
            nodes = Node[]

            if cell.shape in (TRI3, QUAD4)
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
                error("revolve: Cell $(cell.shape.name) is not supported")
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

    # move back points that lied at the axis
    #for node in mesh.nodes
        #X = node.coord
        #if norm( dot(X-base, axis)*axis - X ) < 1.00001e-6
            #node.coord = round.( dot(X-base, axis)*axis , digits=8)
        #end
    #end

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells

    for elem in newmesh.elems
        inaxis = false
        for node in elem.nodes
            X = node.coord
            if norm( dot(X-base, axis)*axis - X ) < 1e-8
                inaxis = true
                break
                node.coord = node.coord + normal*1e-6
            end
        end
        inaxis && collapse!(elem)
        #collapse!(elem)
    end

    fixup!(newmesh, reorder=true)

    return newmesh
end

