# flattens an array of nested arrays
function flatten(x, y)
    ty = typeof(x)
    if ty <: Tuple || ty <: Array
        for item in x
            flatten(item, y)
        end
    else
        push!(y, x)
    end
    return y
end
flatten(x)=flatten(x, [])


function mesh_structured(geo::GeoModel; args...)
    mesh = Mesh() # empty mesh
    return mesh_structured(mesh, geo.blocks; args...)
end


function mesh_structured(mesh::Mesh, block_or_arr...; args...)
    blocks = flatten(block_or_arr)

    for b in blocks
        b isa block || error("mesh_structured: expected Block or Array{Block}, got $(typeof(b))")
        split_block!(mesh, b)
    end

    fixup!(mesh)

    return mesh
end


# Splits a block
function split_block!(mesh::Mesh, bl::Block)
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    rx, ry, rz = bl.rx, bl.ry, bl.rz
    coords = getcoords(bl.points)
    shape = bl.cellshape

    if shape==LIN2
        p_arr = Array{Node}(undef, nx+1)
        for i = 1:nx+1
            # r = (2.0/nx)*(i-1) - 1.0
            r = -1.0 + 2.0*(rx==1 ? (1/nx)*(i-1) : (1-rx^(i-1))/(1-rx^nx))
            N = bl.shape.func([r])
            C = N'*coords
            C = round.(C, digits=8)
            p = get_node(mesh._pointdict, C)
            if p===nothing
                p = Node(C);
                push!(mesh.nodes, p)
                mesh._pointdict[hash(p)] = p
            end
            p_arr[i] = p
        end

        for i = 1:nx
            p1 = p_arr[i  ]
            p2 = p_arr[i+1]

            cell = Cell(shape, [p1, p2], tag=bl.tag, env=mesh.env)
            push!(mesh.elems, cell)
        end
        return
    end

    if shape==LIN3
        p_arr = Array{Node}(undef, 2*nx+1)
            for i = 1:2*nx+1
                # r = (1.0/nx)*(i-1) - 1.0
                r = -1.0 + 2.0*(rx==1 ? (1/(2*nx))*(i-1) : (1-rx^(i-1))/(1-rx^(2*nx)))
                N = bl.shape.func([r])
                C = N'*coords
                C = round.(C, digits=8)
                p = get_node(mesh._pointdict, C)
                if p===nothing
                    p = Node(C); 
                    push!(mesh.nodes, p)
                    mesh._pointdict[hash(p)] = p
                end
                p_arr[i] = p
            end

            for i = 1:2:2*nx
                p1 = p_arr[i  ]
                p2 = p_arr[i+2]
                p3 = p_arr[i+1]

                cell = Cell(shape, [p1, p2, p3], tag=bl.tag, env=mesh.env)
                push!(mesh.elems, cell)
            end
        return
    end

    if shape==LIN4
        p_arr = Array{Node}(undef, 3*nx+1)
            for i = 1:3*nx+1
                r = -1.0 + 2.0*(rx==1 ? (1/(3*nx))*(i-1) : (1-rx^(i-1))/(1-rx^(3*nx)))
                N = bl.shape.func([r])
                C = N'*coords
                C = round.(C, digits=8)
                p = get_node(mesh._pointdict, C)
                if p===nothing
                    p = Node(C); 
                    push!(mesh.nodes, p)
                    mesh._pointdict[hash(p)] = p
                end
                p_arr[i] = p
            end

            for i = 1:3:3*nx
                p1 = p_arr[i  ]
                p2 = p_arr[i+3]
                p3 = p_arr[i+1]
                p4 = p_arr[i+2]

                cell = Cell(shape, [p1, p2, p3, p4], tag=bl.tag, env=mesh.env)
                push!(mesh.elems, cell)
            end
        return
    end

    if shape in (QUAD4, TRI3)
        p_arr = Array{Node}(undef, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                # r = (2.0/nx)*(i-1) - 1.0
                # s = (2.0/ny)*(j-1) - 1.0 
                r = -1.0 + 2.0*(rx==1 ? (1/nx)*(i-1) : (1-rx^(i-1))/(1-rx^nx))
                s = -1.0 + 2.0*(ry==1 ? (1/ny)*(j-1) : (1-ry^(j-1))/(1-ry^ny))
                
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, digits=8)
                    p =get_node(mesh._pointdict, C)
                    if p===nothing
                        p = Node(C); push!(mesh.nodes, p)
                        mesh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C);
                    push!(mesh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                if shape==QUAD4
                    cell = Cell(shape, [p1, p2, p3, p4], tag=bl.tag, env=mesh.env)
                    push!(mesh.elems, cell)
                else
                    C = (p1.coord+p2.coord+p3.coord+p4.coord)/4
                    p5 = Node(C); push!(mesh.nodes, p5)
                    mesh._pointdict[hash(p5)] = p5
                    
                    cell1 = Cell(shape, [p1, p2, p5], tag=bl.tag, env=mesh.env)
                    cell2 = Cell(shape, [p2, p3, p5], tag=bl.tag, env=mesh.env)
                    cell3 = Cell(shape, [p3, p4, p5], tag=bl.tag, env=mesh.env)
                    cell4 = Cell(shape, [p4, p1, p5], tag=bl.tag, env=mesh.env)
                    push!(mesh.elems, cell1)
                    push!(mesh.elems, cell2)
                    push!(mesh.elems, cell3)
                    push!(mesh.elems, cell4)
                end
            end
        end
        return
    end

    if shape in (QUAD8, QUAD9)
        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                shape==QUAD8 && iseven(i) && iseven(j) && continue

                # r = (1.0/nx)*(i-1) - 1.0
                # s = (1.0/ny)*(j-1) - 1.0
                r = -1.0 + 2.0*(rx==1 ? (1/(2*nx))*(i-1) : (1-rx^(i-1))/(1-rx^(2*nx)))
                s = -1.0 + 2.0*(ry==1 ? (1/(2*ny))*(j-1) : (1-ry^(j-1))/(1-ry^(2*ny)))
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(mesh._pointdict, C)
                    if p===nothing
                        p = Node(C); push!(mesh.nodes, p)
                        mesh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(mesh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                if shape==QUAD8
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag, env=mesh.env)
                else
                    p9   = p_arr[i+1, j+1]
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9], tag=bl.tag, env=mesh.env)
                end
                push!(mesh.elems, cell)
            end
        end
        return
    end

    if shape == QUAD12
        p_arr = Array{Node}(undef, 3*nx+1, 3*ny+1)
        for j = 1:3*ny+1
            for i = 1:3*nx+1
                shape==QUAD12 && (i-1)%3>0 && (j-1)%3>0 && continue

                # r = ((2/3)/nx)*(i-1) - 1.0
                # s = ((2/3)/ny)*(j-1) - 1.0
                r = -1.0 + 2.0*(rx==1 ? (1/(3*nx))*(i-1) : (1-rx^(i-1))/(1-rx^(3*nx)))
                s = -1.0 + 2.0*(ry==1 ? (1/(3*ny))*(j-1) : (1-ry^(j-1))/(1-ry^(3*ny)))

                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 3*nx+1) || j in (1, 3*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(mesh._pointdict, C)
                    if p===nothing
                        p = Node(C); push!(mesh.nodes, p)
                        mesh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(mesh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:3:3*ny
            for i = 1:3:3*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+3, j  ]
                p3 = p_arr[i+3, j+3]
                p4 = p_arr[i  , j+3]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j  ]
                p7 = p_arr[i+3, j+1]
                p8 = p_arr[i+3, j+2]
                p9 = p_arr[i+2, j+3]
                p10 = p_arr[i+1, j+3]
                p11 = p_arr[i  , j+2]
                p12 = p_arr[i  , j+1]

                cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], tag=bl.tag, env=mesh.env)
                push!(mesh.elems, cell)
            end
        end
        return
    end

    # if shape == TRI3
    #     p_arr = Array{Node}(undef, nx+1, ny+1)
    #     for j = 1:ny+1
    #         for i = 1:nx+1

    #             # r = (2.0/nx)*(i-1) - 1.0
    #             # s = (2.0/ny)*(j-1) - 1.0
    #             r = -1.0 + 2.0*(rx==1 ? (1/nx)*(i-1) : (1-rx^(i-1))/(1-rx^nx))
    #             s = -1.0 + 2.0*(ry==1 ? (1/ny)*(j-1) : (1-ry^(j-1))/(1-ry^ny))
                
    #             N = bl.shape.func([r, s])
    #             C = N'*coords
    #             p::Any = nothing
    #             if i in (1, nx+1) || j in (1, ny+1)
    #                 C = round.(C, digits=8)
    #                 p =get_node(mesh._pointdict, C)
    #                 if p===nothing
    #                     p = Node(C); push!(mesh.nodes, p)
    #                     mesh._pointdict[hash(p)] = p
    #                 end
    #             else
    #                 p = Node(C); push!(mesh.nodes, p)
    #             end
    #             p_arr[i,j] = p
    #         end
    #     end

    #     for j = 1:ny
    #         for i = 1:nx
    #             p1 = p_arr[i  , j  ]
    #             p2 = p_arr[i+1, j  ]
    #             p3 = p_arr[i+1, j+1]
    #             p4 = p_arr[i  , j+1]

    #             cell1 = Cell(shape, [p1, p2, p3], tag=bl.tag, env=mesh.env)
    #             cell2 = Cell(shape, [p4, p1, p3], tag=bl.tag, env=mesh.env)
    #             push!(mesh.elems, cell1)
    #             push!(mesh.elems, cell2)
    #         end
    #     end
    #     return
    # end

    if shape == TRI6

        #=   4       7       3
               @-----@-----@
               |         / |
               |       /   |
             8 @     @     @ 6
               |   /  9    |
               | /         |
               @-----@-----@
             1       5       2     =#

        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
                    p =get_node(mesh._pointdict, C)
                    if p===nothing
                        p = Node(C); push!(mesh.nodes, p)
                        mesh._pointdict[hash(p)] = p
                    end
                else
                    p = Node(C); push!(mesh.nodes, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                p9   = p_arr[i+1, j+1]

                cell1 = Cell(shape, [p1, p2, p3, p5, p6, p9], tag=bl.tag, env=mesh.env)
                cell2 = Cell(shape, [p4, p1, p3, p8, p9, p7], tag=bl.tag, env=mesh.env)
                push!(mesh.elems, cell1)
                push!(mesh.elems, cell2)
            end
        end
        return
    end

    if shape in (HEX8, TET4, PYR5)
        p_arr = Array{Node}(undef, nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    # r = (2.0/nx)*(i-1) - 1.0
                    # s = (2.0/ny)*(j-1) - 1.0
                    # t = (2.0/nz)*(k-1) - 1.0
                    r = -1.0 + 2.0*(rx==1 ? (1/nx)*(i-1) : (1-rx^(i-1))/(1-rx^nx))
                    s = -1.0 + 2.0*(ry==1 ? (1/ny)*(j-1) : (1-ry^(j-1))/(1-ry^ny))
                    t = -1.0 + 2.0*(rz==1 ? (1/nz)*(k-1) : (1-rz^(k-1))/(1-rz^nz))
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        C = round.(C, digits=8)
                        p =get_node(mesh._pointdict, C)
                        if p===nothing
                            p = Node(C); push!(mesh.nodes, p)
                            mesh._pointdict[hash(p)] = p
                        end
                    else
                        p = Node(C); push!(mesh.nodes, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    p1 = p_arr[i  , j  , k  ]
                    p2 = p_arr[i+1, j  , k  ]
                    p3 = p_arr[i+1, j+1, k  ]
                    p4 = p_arr[i  , j+1, k  ]
                    p5 = p_arr[i  , j  , k+1]
                    p6 = p_arr[i+1, j  , k+1]
                    p7 = p_arr[i+1, j+1, k+1]
                    p8 = p_arr[i  , j+1, k+1]

                    if shape==HEX8
                        cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag, env=mesh.env)
                        push!(mesh.elems, cell)
                    end
                    if shape==TET4
                        push!( mesh.elems, Cell(shape, [p2, p4, p1, p8], tag=bl.tag, env=mesh.env) )
                        push!( mesh.elems, Cell(shape, [p2, p1, p5, p8], tag=bl.tag, env=mesh.env) )
                        push!( mesh.elems, Cell(shape, [p2, p5, p6, p8], tag=bl.tag, env=mesh.env) )
                        push!( mesh.elems, Cell(shape, [p2, p6, p7, p8], tag=bl.tag, env=mesh.env) )
                        push!( mesh.elems, Cell(shape, [p2, p3, p4, p8], tag=bl.tag, env=mesh.env) )
                        push!( mesh.elems, Cell(shape, [p2, p7, p3, p8], tag=bl.tag, env=mesh.env) )
                    end
                    if shape==PYR5
                        C = (p1.coord+p2.coord+p3.coord+p4.coord+p5.coord+p6.coord+p7.coord+p8.coord)/8
                        p9 = Node(C); push!(mesh.nodes, p9)
                        mesh._pointdict[hash(p9)] = p9
                        
                        cell1 = Cell(shape, [p1, p2, p3, p4, p9], tag=bl.tag, env=mesh.env)
                        cell2 = Cell(shape, [p2, p6, p7, p3, p9], tag=bl.tag, env=mesh.env)
                        cell3 = Cell(shape, [p4, p3, p7, p8, p9], tag=bl.tag, env=mesh.env)
                        cell4 = Cell(shape, [p1, p4, p8, p5, p9], tag=bl.tag, env=mesh.env)
                        cell5 = Cell(shape, [p2, p1, p5, p6, p9], tag=bl.tag, env=mesh.env)
                        cell6 = Cell(shape, [p6, p5, p8, p7, p9], tag=bl.tag, env=mesh.env)
                        push!(mesh.elems, cell1)
                        push!(mesh.elems, cell2)
                        push!(mesh.elems, cell3)
                        push!(mesh.elems, cell4)
                        push!(mesh.elems, cell5)
                        push!(mesh.elems, cell6)
                    end
                end
            end
        end
        return
    end

    if shape in (TET10, HEX20, HEX27)
        p_arr = Array{Node}(undef, 2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if shape==HEX20
                        if iseven(i) && iseven(j) continue end
                        if iseven(j) && iseven(k) continue end
                        if iseven(k) && iseven(i) continue end
                    end

                    # r = (1.0/nx)*(i-1) - 1.0
                    # s = (1.0/ny)*(j-1) - 1.0
                    # t = (1.0/nz)*(k-1) - 1.0
                    r = -1.0 + 2.0*(rx==1 ? (1/(2*nx))*(i-1) : (1-rx^(i-1))/(1-rx^(2*nx)))
                    s = -1.0 + 2.0*(ry==1 ? (1/(2*ny))*(j-1) : (1-ry^(j-1))/(1-ry^(2*ny)))
                    t = -1.0 + 2.0*(rz==1 ? (1/(2*nz))*(k-1) : (1-rz^(k-1))/(1-rz^(2*nz)))
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        C = round.(C, digits=8)
                        p =get_node(mesh._pointdict, C)
                        if p===nothing
                            p = Node(C); push!(mesh.nodes, p)
                            mesh._pointdict[hash(p)] = p
                        end
                    else
                        p = Node(C); push!(mesh.nodes, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:2:2*nz
            for j = 1:2:2*ny
                for i = 1:2:2*nx
                    p1  = p_arr[i  , j  , k  ]
                    p2  = p_arr[i+2, j  , k  ]
                    p3  = p_arr[i+2, j+2, k  ]
                    p4  = p_arr[i  , j+2, k  ]
                    p5  = p_arr[i  , j  , k+2]
                    p6  = p_arr[i+2, j  , k+2]
                    p7  = p_arr[i+2, j+2, k+2]
                    p8  = p_arr[i  , j+2, k+2]

                    p9  = p_arr[i+1, j  , k  ]
                    p10 = p_arr[i+2, j+1, k  ]
                    p11 = p_arr[i+1, j+2, k  ]
                    p12 = p_arr[i  , j+1, k  ]
                    p13 = p_arr[i+1, j  , k+2]
                    p14 = p_arr[i+2, j+1, k+2]
                    p15 = p_arr[i+1, j+2, k+2]
                    p16 = p_arr[i  , j+1, k+2]

                    p17 = p_arr[i  , j  , k+1]
                    p18 = p_arr[i+2, j  , k+1]
                    p19 = p_arr[i+2, j+2, k+1]
                    p20 = p_arr[i  , j+2, k+1]

                    if shape == HEX20
                        cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20], tag=bl.tag, env=mesh.env)
                        push!(mesh.elems, cell)
                    end
                    if shape in (TET10, HEX27)

                        p21 = p_arr[i  , j+1, k+1]
                        p22 = p_arr[i+2, j+1, k+1]
                        p23 = p_arr[i+1, j  , k+1]
                        p24 = p_arr[i+1, j+2, k+1]
                        p25 = p_arr[i+1, j+1, k  ]
                        p26 = p_arr[i+1, j+1, k+2]
                        p27 = p_arr[i+1, j+1, k+1]

                        if shape==TET10
                            push!( mesh.elems, Cell(shape, [p2, p4, p1, p8, p25, p12, p9,  p27, p20, p21], tag=bl.tag, env=mesh.env) )
                            push!( mesh.elems, Cell(shape, [p2, p1, p5, p8, p9,  p17, p23, p27, p21, p16], tag=bl.tag, env=mesh.env) )
                            push!( mesh.elems, Cell(shape, [p2, p5, p6, p8, p23, p13, p18, p27, p16, p26], tag=bl.tag, env=mesh.env) )
                            push!( mesh.elems, Cell(shape, [p2, p6, p7, p8, p18, p14, p22, p27, p26, p15], tag=bl.tag, env=mesh.env) )
                            push!( mesh.elems, Cell(shape, [p2, p3, p4, p8, p10, p11, p25, p27, p24, p20], tag=bl.tag, env=mesh.env) )
                            push!( mesh.elems, Cell(shape, [p2, p7, p3, p8, p22, p19, p10, p27, p15, p24], tag=bl.tag, env=mesh.env) )
                        else
                            cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27], tag=bl.tag, env=mesh.env)
                            push!(mesh.elems, cell)
                        end
                    end
                end
            end
        end
        return
    end
    error("block: Can not discretize using shape $(shape.name)")
end
