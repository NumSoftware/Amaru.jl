function hrefine(mesh::Mesh; n=2, verbose=true)
    msh = Mesh()

    for cell in mesh.cells
        if cell.shape==TRI3
            coords = getcoords(cell.points)

            p_arr = Array{Point}(undef, n+1, n+1)
            for j = 1:n+1
                for i = 1:n+1
                    i+j > n+2 && continue
                    r = (1.0/n)*(i-1)
                    s = (1.0/n)*(j-1)

                    N = cell.shape.func([r, s])
                    C = N'*coords
                    if i==1 || j==1 || i+j==n+2
                        C = round.(C, digits=8)
                        p = get_point(msh.pointdict, C)
                        if p==nothing
                            p = Point(C); push!(msh.points, p)
                            msh.pointdict[hash(p)] = p
                        end
                    else
                        p = Point(C); push!(msh.points, p)
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
                    push!(msh.cells, cell1)

                    if i+j < n+1
                        p4 = p_arr[i+1, j+1]
                        cell2 = Cell(cell.shape, [p2, p4, p3], tag=cell.tag)
                        push!(msh.cells, cell2)
                    end
                end
            end
        end
    end
    fixup!(msh, reorder=true)
    return msh
end


function prefine(mesh::Mesh; n=2, verbose=true)
    msh = Mesh()

    for cell in mesh.cells
        if cell.shape==TRI3
            NS = TRI6 # new shape
            coords = getcoords(cell.points)
            points = Point[]
            for i=1:NS.npoints
                R = NS.nat_coords[i,:]
                N = cell.shape.func(R)
                C = coords'*N
                C = round.(C, digits=8)
                p = get_point(msh.pointdict, C)
                if p==nothing
                    p = Point(C);
                    push!(msh.points, p)
                    msh.pointdict[hash(p)] = p
                end
                push!(points, p)
            end
            newcell = Cell(NS, points, tag=cell.tag)
            push!(msh.cells, newcell)
        end
    end
    fixup!(msh, reorder=true)
    return msh
end
