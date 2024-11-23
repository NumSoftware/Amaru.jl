# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function gen_insets!(mesh::Mesh, subpaths::Vector{SubPath})
    for subpath in subpaths
        gen_insets!(mesh, subpath)
    end
end


function gen_insets!(mesh::Mesh, subpath::SubPath)
    cmds        = subpath.path.cmds
    closed      = subpath.closed
    tag         = subpath.tag
    jointtag    = subpath.jointtag
    tiptag = subpath.tiptag
    tips    = subpath.tips
    embedded    = subpath.embedded
    shape       = subpath.shape

    # TODO: add option: merge_points
    # TODO: check closed when hole 

    # tolerances
    ε  = 1e-9  # for bissection
    εn = 1e-4  # for checking first cell and end of path
    εc = 1e-9  # for checking if point is inside a cell
    λ  = 1.0   # step for checking multiple intersections in one cell (0 < λ < 1)

    npoints = subpath.shape==LIN2 ? 2 : 3
    jntshape = subpath.shape==LIN2 ? JLINK2 : JLINK3

    # Initial conditions
    len = 1.0

    firstnode = Node(cmds[1].points[1].coord)
    lastnode = nothing
    for cmd in cmds
        cmd.key==:M && continue

        # find the initial element
        X0    = cmd(εn) # a little bit ahead from 0.0
        ccell = find_elem(X0, mesh.elems, mesh._elempartition, εc) # the first tresspased cell

        first_segment = true
        last_segment  = false
        s  = 0.0
        s1 = 0.0

        # Splitting cmd
        k = 0
        while true
            k +=1

            if ccell !== nothing
                ccell_coords = getcoords(ccell)
                # Default step
                step  = 0.5*(1.0-s)

                # Finding step (introduced to catch elements in curved paths)
                str = s     # trial point
                nits = round(Int, 1.0/λ)
                for i in 1:nits
                    str += λ
                    str>1.0 && break
                    X = cmd(str)
                    if !is_inside(ccell.shape, ccell_coords, X, εc)
                        step = 0.5*(str-s)
                        break
                    end
                end

                s += step
                X  = cmd(s)
                n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

                for i in 1:n
                    step *= 0.5
                    if is_inside(ccell.shape, ccell_coords, X, εc)
                        s += step
                    else
                        s -= step
                    end

                    X = cmd(s)
                end
            else # ccell is nothing (hole or gap)
                step  = 0.5*(1.0-s)
                s += step
                X  = cmd(s)
                n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

                for i in 1:n
                    step *= 0.5
                    ccell = find_elem(X, mesh.elems, mesh._elempartition, εc)
                    if ccell === nothing
                        s += step
                    else
                        s -= step
                    end

                    X = cmd(s)
                end
            end

            # Check if end was reached
            if s > len - εn
                last_segment = true
            end

            # Getting line cell nodes
            if lastnode===nothing
                P1 = firstnode
                push!(mesh.nodes, P1)
            else
                P1 = lastnode
            end

            if last_segment && closed
                P2 = firstnode
            else 
                P2 = Node(X)
                push!(mesh.nodes, P2)
            end

            if npoints==2
                points = [P1, P2]
            else
                P3 = Node(cmd((s1+s)/2))
                push!(mesh.nodes, P3)
                points = [P1, P2, P3]
            end

            # Saving line cell
            lcell = Cell(shape, points, tag=tag)
            push!(mesh.elems, lcell)

            if ccell!== nothing
                if embedded
                    # Set line as embedded
                    lcell.embedded = true
                    lcell.linked_elems = [ ccell ]
                else
                    # Generate a continuous joint element
                    jntpts  = vcat( ccell.nodes, lcell.nodes )
                    jntcell = Cell(jntshape, jntpts, tag=jointtag)
                    push!(mesh.elems, jntcell)
                    jntcell.linked_elems = [ccell, lcell]

                    # generate tip joints
                    if first_segment && tips in (:front, :both)
                        tip = P1
                        tipjointnodes = vcat(ccell.nodes, tip)
                        tipjointcell = Cell(tips, tipjointnodes, tag=tiptag)
                        tipjointcell.linked_elems = jntcell.linked_elems
                        push!(mesh.elems, tipjointcell)
                    end
                    if last_segment && tips in (:end, :both)
                        tip = P2
                        tipjointnodes = vcat(ccell.nodes, tip)
                        tipjointcell = Cell(tips, tipjointnodes, tag=tiptag)
                        tipjointcell.linked_elems = jntcell.linked_elems
                        push!(mesh.elems, tipjointcell)
                    end
                end
                ccell.crossed = true
            end

            # update first and last points
            if firstnode === nothing;
                firstnode = P1 
            end
            lastnode = P2
            first_segment = false

            last_segment && break

            # Preparing for the next iteration
            ccell = find_elem(cmd(s + εn), mesh.elems, mesh._elempartition, εc, exclude=[ccell])
            s1    = s
            s     = s + εn
        end

    end

end

