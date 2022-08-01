function hrefine(mesh::Mesh; n=2, printlog=false)
    n==1 && return copy(mesh)
    newmesh = Mesh()

    for cell in mesh.elems
        coords =getcoords(cell.nodes)

        if cell.shape==TRI3
            p_arr = Array{Node}(undef, n+1, n+1)
            for j = 1:n+1
                for i = 1:n+1
                    i+j > n+2 && continue
                    r = (1.0/n)*(i-1)
                    s = (1.0/n)*(j-1)

                    N = cell.shape.func([r, s])
                    C = round.(N'*coords, digits=8)
                    if i==1 || j==1 || i+j==n+2
                        p = get_node(newmesh._pointdict, C)
                        if p===nothing
                            p = Node(C); push!(newmesh.nodes, p)
                            newmesh._pointdict[hash(p)] = p
                        end
                    else
                        p = Node(C); push!(newmesh.nodes, p)
                    end
                    p_arr[i,j] = p
                end
            end

            for j = 1:n
                for i = 1:n
                    i+j >= n+2 && continue
                    p1 = p_arr[i  , j  ]
                    p2 = p_arr[i+1, j  ]
                    p3 = p_arr[i  , j+1]

                    cell1 = Cell(cell.shape, [p1, p2, p3], tag=cell.tag)
                    push!(newmesh.elems, cell1)

                    if i+j < n+1
                        p4 = p_arr[i+1, j+1]
                        cell2 = Cell(cell.shape, [p2, p4, p3], tag=cell.tag)
                        push!(newmesh.elems, cell2)
                    end
                end
            end
            continue
        end

        if cell.shape==HEX8
            p_arr = Array{Node}(undef, n+1, n+1, n+1)
            for k = 1:n+1
                for j = 1:n+1
                    for i = 1:n+1
                        r = (2.0/n)*(i-1) - 1.0
                        s = (2.0/n)*(j-1) - 1.0
                        t = (2.0/n)*(k-1) - 1.0
                        N = cell.shape.func([r, s, t])
                        C = round.(N'*coords, digits=8)
                        if i in (1, n+1) || j in (1, n+1) || k in (1, n+1)
                            p =get_node(newmesh._pointdict, C)
                            if p===nothing
                                p = Node(C); push!(newmesh.nodes, p)
                                newmesh._pointdict[hash(p)] = p
                            end
                        else
                            p = Node(C); push!(newmesh.nodes, p)
                        end
                        p_arr[i,j,k] = p
                    end
                end
            end

            for k = 1:n
                for j = 1:n
                    for i = 1:n
                        p1 = p_arr[i  , j  , k  ]
                        p2 = p_arr[i+1, j  , k  ]
                        p3 = p_arr[i+1, j+1, k  ]
                        p4 = p_arr[i  , j+1, k  ]
                        p5 = p_arr[i  , j  , k+1]
                        p6 = p_arr[i+1, j  , k+1]
                        p7 = p_arr[i+1, j+1, k+1]
                        p8 = p_arr[i  , j+1, k+1]

                        cell = Cell(cell.shape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=cell.tag)
                        push!(newmesh.elems, cell)
                    end
                end
            end
            continue
        end

        error("hrefine: Cannot refine mesh containing elements of type $(cell.shape.name)")

    end

    fixup!(newmesh, reorder=true)
    return newmesh
end


function prefine(mesh::Mesh; n=2, printlog=false)
    #newmesh = Mesh()

    NS = Dict{CellShape,CellShape}(TRI3=>TRI6, TET4=>TET10)

    # Generate new cells
    cells = Cell[]

    for cell in mesh.elems
        newshape = NS[cell.shape]
        coords = getcoords(cell.nodes)
        points = Node[]
        for i in 1:newshape.npoints
            R = newshape.nat_coords[i,:]
            N = cell.shape.func(R)
            C = coords'*N
            p = Node(C);
            push!(points, p)
        end
        newcell = Cell(newshape, points, tag=cell.tag)
        push!(cells, newcell)
    end

    # Merge points
    node_dict = Dict{UInt64,Node}( hash(n) => n for c in cells for n in c.nodes )

    for cell in cells
        cell.nodes = Node[ node_dict[hash(n)] for n in cell.nodes ]
    end
    nodes = collect(values(node_dict))

    # New mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    fixup!(newmesh, reorder=true)

    return newmesh
end
