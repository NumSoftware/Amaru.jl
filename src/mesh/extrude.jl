# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

"""
`extrude(block, [axis=[0,0,1],] [len=1.0,] [n=1])`

Extrudes a 2D `block` generating a 3D block based on a direction `axis`, a lenght `len` and a number `n` of divisions.
"""
function extrude(block::Block2D; axis=[0,0,1], len::Number=1.0, n::Int=1)::Block3D

    V = axis/norm(axis)
    δ = len/n
    #coords = block.coords
    coords = getcoords(block.points)

    # check numbering order
    v1 = vec(coords[1,:])
    v2 = vec(coords[2,:])
    v3 = vec(coords[3,:])
    normal_order = dot(cross(v2-v1, v3-v2), axis) > 0.0

    # generate extra points
    npoints = size(coords,1)
    topcoords = coords .+ V'*len

    local newcoords::Array{Float64,2}

    if npoints==4
        if normal_order
            newcoords = vcat(coords, topcoords)
        else
            newcoords = vcat(topcoords, coords)
        end
    end

    if npoints==8
        midcoords = coords[1:4,:] .+ V'*len*0.5
        if normal_order
            newcoords = vcat(coords[1:4,:], topcoords[1:4,:], coords[5:8,:], topcoords[5:8,:], midcoords)
        else
            newcoords = vcat(coords[5:8,:], topcoords[5:8,:], coords[1:4,:], topcoords[1:4,:], midcoords)
        end
    end

    #shape = HEX8
    #if block.shape==TRI3
        #shape = WED6
    #elseif block.shape==TRI6
        #shape = WED15
    #elseif block.shape==QUAD8
        #shape = HEX20
    #end

    shape = block.cellshape==QUAD8 ? HEX20 : HEX8

    return Block3D( newcoords, nx=block.nx, ny=block.ny, nz=n, cellshape=shape)

end

function extrude(blocks::Array; axis=[0,0,1], len=1.0::Number, n=1::Int)::Array{Block3D,1}
    blocks3D = []

    for bl in blocks
        bl3D = extrude(bl, axis=axis, len=len, n=n)
        push!(blocks3D, bl3D)
    end

    return blocks3D
end

# Generates a new mesh obtained by extrusion of a 2D mesh
"""
`extrude(mesh, [axis=[0,0,1],] [len=1.0,] [n=1,] [verbose=true,] [genedges=false])`

Extrudes a 2D `mesh` generating a 3D mesh based on a direction `axis`, a lenght `len` and a number `n` of divisions.
Only meshes with cell shapes TRI3, TRI6, QUAD4 and QUAD8 are allowed.
"""
function extrude(mesh::Mesh; len::Number=1.0, n::Int=1, verbose::Bool=true)
    verbose && printstyled("Mesh extrude:\n", bold=true, color=:cyan)
    Δz = len/n

    # Generate new cells
    cells = Cell[]

    for cell in mesh.cells
        if cell.shape in (TRI3, QUAD4)
            newshape = cell.shape==TRI3 ? WED6 : HEX8
            for i=1:n
                zi = (i-1)*Δz
                points  = [ Point(p.x, p.y, z) for z in (zi, zi+Δz) for p in cell.points ]
                newcell = Cell(newshape, points, tag=cell.tag)
                push!(cells, newcell)
            end
        elseif cell.shape==TRI6
            for i=1:n
                zi = (i-1)*Δz
                points1 = [ Point(p.x, p.y, z) for z in (zi, zi+Δz) for p in view(cell.points, 1:3) ]
                points2 = [ Point(p.x, p.y, z) for z in (zi, zi+Δz) for p in view(cell.points, 4:6) ]
                points3 = [ Point(p.x, p.y, zi+Δz/2) for p in view(cell.points, 1:3) ]
                newcell = Cell(WED15, [ points1; points2; points3 ], tag=cell.tag)
                push!(cells, newcell)
            end
        elseif cell.shape==QUAD8
            for i=1:n
                zi = (i-1)*Δz
                points1 = [ Point(p.x, p.y, z) for z in (zi, zi+Δz) for p in view(cell.points, 1:4) ]
                points2 = [ Point(p.x, p.y, z) for z in (zi, zi+Δz) for p in view(cell.points, 5:8) ]
                points3 = [ Point(p.x, p.y, zi+Δz/2) for p in view(cell.points, 1:4) ]
                newcell = Cell(WED6, [ points1; points2; points3 ], tag=cell.tag)
                push!(cells, newcell)
            end
        else
            error("extrude: Cell $(cell.shape.name) is not supported")
        end
    end

    # Merge points
    pointsD = Dict{UInt64,Point}( hash(p) => p for c in cells for p in c.points )

    for cell in cells
        cell.points = [ pointsD[hash(p)] for p in cell.points ]
    end

    # New mesh
    newmesh = Mesh()
    newmesh.points = collect(values(pointsD))
    newmesh.cells  = cells
    fixup!(newmesh, reorder=true)

    if verbose
        @printf "  %5d points obtained\n" length(newmesh.points)
        @printf "  %5d cells obtained\n" length(newmesh.cells)
        @printf "  %5d faces obtained\n" length(newmesh.faces)
        @printf "  %5d surface edges obtained\n" length(newmesh.edges)
        println("  done.")
    end

    return newmesh

end



function extrude2(mesh::Mesh; axis=[0.0,0.0,1.], len::Number=1.0, n::Int=1, verbose::Bool=true, genedges::Bool=false)
    verbose && printstyled("Mesh extrude:\n", bold=true, color=:cyan)

    V = axis/norm(axis)
    δ = len/n
    inipoints = mesh.points
    inicells  = mesh.cells
    newmesh   = Mesh()

    length(inicells)>0 || error("Extrude: Cannot extrude mesh with no cells.")

    # check if all cells are the same
    shape = inicells[1].shape
    any(Bool[ c.shape!=shape for c in inicells ]) && error("Extrude: Input mesh shoud have same type of cells.")

    # check if cells are QUAD4 or QUAD8
    (shape == QUAD4 || shape == QUAD8) || error("Error: Only can extrude meshes of QUAD4 and QUAD8 cells.")

    # check numbering order
    p1, p2, p3 = inicells[1].points[1:3]
    v1 = [ p1.x, p1.y, p1.z ]
    v2 = [ p2.x, p2.y, p2.z ]
    v3 = [ p3.x, p3.y, p3.z ]
    normal_order = dot(cross(v2-v1, v3-v2), axis) > 0.0

    # generate extra nodes and cells
    if shape == QUAD4
        # generate new points
        p_arr = [ Array{Point}(n+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                C = X + V*δ*(j-1)
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array{Cell}( n) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:n
                cpoints = Array{Point}( 8)
                cpoints[1] = p_arr[ c.points[1].id ][j]
                cpoints[2] = p_arr[ c.points[2].id ][j]
                cpoints[3] = p_arr[ c.points[3].id ][j]
                cpoints[4] = p_arr[ c.points[4].id ][j]
                cpoints[5] = p_arr[ c.points[1].id ][j+1]
                cpoints[6] = p_arr[ c.points[2].id ][j+1]
                cpoints[7] = p_arr[ c.points[3].id ][j+1]
                cpoints[8] = p_arr[ c.points[4].id ][j+1]

                if !normal_order
                    cpoints = vcat(cpoints[5:8], cpoints[1:4])
                end

                cell = Cell(HEX8, cpoints, tag=c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    # generate extra nodes and cells
    if shape == QUAD8

        # flags for middle points
        middle = trues(length(inipoints))
        for c in inicells
            for i=1:4
                middle[ c.points[i].id ] = false
            end
        end

        # generate new points
        p_arr = [ Array{Point}(2*n+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                (middle[i] && iseven(j)) && continue
                C = X + V*δ*(j-1)*0.5
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array{Cell}( n) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:2:2*n
                cpoints = Array{Point}( 20)
                cpoints[1 ] = p_arr[ c.points[1].id ][j]
                cpoints[2 ] = p_arr[ c.points[2].id ][j]
                cpoints[3 ] = p_arr[ c.points[3].id ][j]
                cpoints[4 ] = p_arr[ c.points[4].id ][j]
                cpoints[5 ] = p_arr[ c.points[1].id ][j+2]
                cpoints[6 ] = p_arr[ c.points[2].id ][j+2]
                cpoints[7 ] = p_arr[ c.points[3].id ][j+2]
                cpoints[8 ] = p_arr[ c.points[4].id ][j+2]
                cpoints[9 ] = p_arr[ c.points[5].id ][j]
                cpoints[10] = p_arr[ c.points[6].id ][j]
                cpoints[11] = p_arr[ c.points[7].id ][j]
                cpoints[12] = p_arr[ c.points[8].id ][j]
                cpoints[13] = p_arr[ c.points[5].id ][j+2]
                cpoints[14] = p_arr[ c.points[6].id ][j+2]
                cpoints[15] = p_arr[ c.points[7].id ][j+2]
                cpoints[16] = p_arr[ c.points[8].id ][j+2]
                cpoints[17] = p_arr[ c.points[1].id ][j+1]
                cpoints[18] = p_arr[ c.points[2].id ][j+1]
                cpoints[19] = p_arr[ c.points[3].id ][j+1]
                cpoints[20] = p_arr[ c.points[4].id ][j+1]

                if !normal_order
                    cpoints = vcat(cpoints[5:8], cpoints[1:4], cpoints[13:16], cpoints[9:12], cpoints[17:20])
                end

                cell = Cell(HEX20, cpoints, tag=c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    fixup!(newmesh, reorder=true, genfacets=true, genedges=genedges)

    if verbose
        @printf "  %5d points obtained\n" length(newmesh.points)
        @printf "  %5d cells obtained\n" length(newmesh.cells)
        @printf "  %5d faces obtained\n" length(newmesh.faces)
        if genedges
            @printf "  %5d surface edges obtained\n" length(newmesh.edges)
        end
        println("  done.")
    end

    return newmesh
end


function extrude1(mesh::Mesh; axis=[0.0,0.0,1.], len::Number=1.0, n::Int=1, verbose::Bool=true, genedges::Bool=false)
    verbose && printstyled("Mesh extrude:\n", bold=true, color=:cyan)

    length(inicells)>0 || error("Extrude: Cannot extrude mesh with no cells.")

    V = axis/norm(axis)
    δ = len/n
    inipoints = mesh.points
    inicells  = mesh.cells
    newmesh   = Mesh()
    newpoints = Dict{UInt64, Point}()

    # Replicate points
    X = zeros(3)
    for j=1:n
        for point in mesh.points
            X.x = point.x
            X.y = point.y
            X.z = point.z
            Xn = X + δ*V
            newpoints[hash(p)] = p
        end
    end


    # loop along cells
    for cell in inicells
        shape = cell.shape
        if shape == QUAD4
            for j=1:n
                for point in cell.points
                    p
                end
            end
        elseif shape == QUAD8
        elseif shape == TRI3
        elseif shape == TRI6
        else
            error("Extrude: Only QUAD4, QUAD8, TRI3 and TRI6 cells can be extruded.")
        end
    end

    # check if all cells are the same
    shape = inicells[1].shape
    any(Bool[ c.shape!=shape for c in inicells ]) && error("Extrude: Input mesh shoud have same type of cells.")

    # check if cells are QUAD4 or QUAD8
    (shape == QUAD4 || shape == QUAD8) || error("Error: Only can extrude meshes of QUAD4 and QUAD8 cells.")

    # check numbering order
    p1, p2, p3 = inicells[1].points[1:3]
    v1 = [ p1.x, p1.y, p1.z ]
    v2 = [ p2.x, p2.y, p2.z ]
    v3 = [ p3.x, p3.y, p3.z ]
    normal_order = dot(cross(v2-v1, v3-v2), axis) > 0.0

    # generate extra nodes and cells
    if shape == QUAD4
        # generate new points
        p_arr = [ Array{Point}(n+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                C = X + V*δ*(j-1)
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array{Cell}( n) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:n
                cpoints = Array{Point}( 8)
                cpoints[1] = p_arr[ c.points[1].id ][j]
                cpoints[2] = p_arr[ c.points[2].id ][j]
                cpoints[3] = p_arr[ c.points[3].id ][j]
                cpoints[4] = p_arr[ c.points[4].id ][j]
                cpoints[5] = p_arr[ c.points[1].id ][j+1]
                cpoints[6] = p_arr[ c.points[2].id ][j+1]
                cpoints[7] = p_arr[ c.points[3].id ][j+1]
                cpoints[8] = p_arr[ c.points[4].id ][j+1]

                if !normal_order
                    cpoints = vcat(cpoints[5:8], cpoints[1:4])
                end

                cell = Cell(HEX8, cpoints, tag=c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    # generate extra nodes and cells
    if shape == QUAD8

        # flags for middle points
        middle = trues(length(inipoints))
        for c in inicells
            for i=1:4
                middle[ c.points[i].id ] = false
            end
        end

        # generate new points
        p_arr = [ Array{Point}(2*n+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                (middle[i] && iseven(j)) && continue
                C = X + V*δ*(j-1)*0.5
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array{Cell}( n) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:2:2*n
                cpoints = Array{Point}( 20)
                cpoints[1 ] = p_arr[ c.points[1].id ][j]
                cpoints[2 ] = p_arr[ c.points[2].id ][j]
                cpoints[3 ] = p_arr[ c.points[3].id ][j]
                cpoints[4 ] = p_arr[ c.points[4].id ][j]
                cpoints[5 ] = p_arr[ c.points[1].id ][j+2]
                cpoints[6 ] = p_arr[ c.points[2].id ][j+2]
                cpoints[7 ] = p_arr[ c.points[3].id ][j+2]
                cpoints[8 ] = p_arr[ c.points[4].id ][j+2]
                cpoints[9 ] = p_arr[ c.points[5].id ][j]
                cpoints[10] = p_arr[ c.points[6].id ][j]
                cpoints[11] = p_arr[ c.points[7].id ][j]
                cpoints[12] = p_arr[ c.points[8].id ][j]
                cpoints[13] = p_arr[ c.points[5].id ][j+2]
                cpoints[14] = p_arr[ c.points[6].id ][j+2]
                cpoints[15] = p_arr[ c.points[7].id ][j+2]
                cpoints[16] = p_arr[ c.points[8].id ][j+2]
                cpoints[17] = p_arr[ c.points[1].id ][j+1]
                cpoints[18] = p_arr[ c.points[2].id ][j+1]
                cpoints[19] = p_arr[ c.points[3].id ][j+1]
                cpoints[20] = p_arr[ c.points[4].id ][j+1]

                if !normal_order
                    cpoints = vcat(cpoints[5:8], cpoints[1:4], cpoints[13:16], cpoints[9:12], cpoints[17:20])
                end

                cell = Cell(HEX20, cpoints, tag=c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    fixup!(newmesh, genfacets=true, genedges=genedgesm, reorder=true)

    if verbose
        @printf "  %5d points obtained\n" length(newmesh.points)
        @printf "  %5d cells obtained\n" length(newmesh.cells)
        @printf "  %5d faces obtained\n" length(newmesh.faces)
        if genedges
            @printf "  %5d surface edges obtained\n" length(newmesh.edges)
        end
        println("  done.")
    end

    return newmesh
end
