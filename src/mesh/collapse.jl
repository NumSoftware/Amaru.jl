


"""
    collapse!(element)

Collapses an `element` with coincident nodes into a simpler element.

"""
function collapse!(elem::AbstractCell)

    shape = elem.shape
    uniquenodes = unique(elem.nodes)
    nunodes = length(uniquenodes)
    length(elem.nodes)==length(uniquenodes) && return

    collapse_error(shape, nnodes) = error("collapse: cannot collapse $(shape.name) element to $nnodes nodes.")
    shape.ndim==1 && error("collapse: cannot collapse $(shape.name) element.")

    if shape.ndim==2
        if shape==QUAD4
            if nunodes==3
                elem.nodes = uniquenodes
                elem.shape = TRI3
            else
                collapse_error(shape, nunodes)
            end
        elseif shape==QUAD8
            if nunodes==6
                validfaces = filter(face -> face.nodes[1]!=face.nodes[2], get_faces(elem))
                elem.nodes = [ 
                            validfaces[1].nodes[1],
                            validfaces[2].nodes[1],
                            validfaces[3].nodes[1],
                            validfaces[1].nodes[3],
                            validfaces[2].nodes[3],
                            validfaces[3].nodes[3]
                           ]
                elem.shape = TRI6
            else
                collapse_error(shape, nunodes)
            end

        else
            error("collapse: cannot collapse $(shape.name) element.")
        end

    elseif shape.ndim==3
        faces = get_faces(elem)
        validfaces = Face[]
        for face in faces
            funiquenodes = unique(face.nodes)
            nfunodes = length(funiquenodes)
            length(face.nodes)==nfunodes && push!(validfaces, face)
            if (face.shape==QUAD4 && nfunodes==3) || (face.shape==QUAD8 && nfunodes==6)
                collapse!(face)
                push!(validfaces, face)
            end
        end

        if shape==HEX8 && nunodes==6 # => WED6
            bottom = filter(f -> length(f.nodes)==3, validfaces)[1]
            bottom.nodes = bottom.nodes[[1,3,2]]
            edges = get_edges(elem)

            sideedges = Edge[]
            for node in bottom.nodes
                for edge in edges
                    n1, n2 = edge.nodes
                    if node==n1 && !(n2 in bottom.nodes)
                        push!(sideedges, edge)
                        break
                    end
                    if node==n2 && !(n1 in bottom.nodes)
                        edge.nodes[1], edge.nodes[2] = edge.nodes[2], edge.nodes[1]
                        push!(sideedges, edge)
                        break
                    end
                end
            end

            elem.nodes = [
                          bottom.nodes;
                          sideedges[1].nodes[2];
                          sideedges[2].nodes[2];
                          sideedges[3].nodes[2]
                         ]
            elem.shape = WED6
        elseif shape==HEX20 && nunodes==15 # => WED15
            bases  = filter(f -> length(f.nodes)==6, validfaces)
            bottom = bases[1]
            top    = bases[2]
            bottom.nodes = bottom.nodes[[1,3,2,6,5,4]]
            edges  = get_edges(elem)

            sideedges = Edge[]
            for node in bottom.nodes[1:3]
                for edge in edges
                    n1, n2 = edge.nodes
                    if node==n1 && !(n2 in bottom.nodes)
                        push!(sideedges, edge)
                        break
                    end
                    if node==n2 && !(n1 in bottom.nodes)
                        edge.nodes[1], edge.nodes[2] = edge.nodes[2], edge.nodes[1]
                        push!(sideedges, edge)
                        break
                    end
                end
            end

            topcornernodes = [ sideedges[1].nodes[2],
                               sideedges[2].nodes[2],
                               sideedges[3].nodes[2] ]
            topedges = get_edges(top)
            topmiddlenodes = Node[]
            for (i,j) in ((1,2),(2,3),(3,1))
                node1 = topcornernodes[i]
                node2 = topcornernodes[j]
                for edge in topedges
                    n1, n2 = edge.nodes
                    if (node1,node2)==(n1,n2) || (node1,node2)==(n2,n1)
                        push!(topmiddlenodes, edge.nodes[3])
                    end
                end
            end

            elem.nodes = [
                          bottom.nodes[1],
                          bottom.nodes[2],
                          bottom.nodes[3],
                          topcornernodes[1],
                          topcornernodes[2],
                          topcornernodes[3],
                          bottom.nodes[4],
                          bottom.nodes[5],
                          bottom.nodes[6],
                          topmiddlenodes[1],
                          topmiddlenodes[2],
                          topmiddlenodes[3],
                          sideedges[1].nodes[3],
                          sideedges[2].nodes[3],
                          sideedges[3].nodes[3],
                         ]
            elem.shape = WED15
        elseif shape==WED6 && nunodes==4 # => TET4
            bottom = filter(f -> length(f.nodes)==3, validfaces)[1]
            reverse!(bottom.nodes)
            tip = setdiff(uniquenodes, bottom.nodes)[1]

            elem.nodes = [ bottom.nodes; tip ]
            elem.shape = TET4
        elseif shape==WED15 && nunodes==10 # => TET10
            bottom = filter(f -> length(f.nodes)==6, validfaces)[1]
            bottom.nodes = bottom.nodes[[1,3,2,6,5,4]]
            edges = get_edges(elem)

            sideedges = Edge[]
            for node in bottom.nodes[1:3]
                for edge in edges
                    n1, n2 = edge.nodes
                    if node==n1 && !(n2 in bottom.nodes)
                        push!(sideedges, edge)
                        break
                    end
                    if node==n2 && !(n1 in bottom.nodes)
                        edge.nodes[1], edge.nodes[2] = edge.nodes[2], edge.nodes[1]
                        push!(sideedges, edge)
                        break
                    end
                end
            end

            tip = sideedges[1].nodes[2]

            elem.nodes = [
                          bottom.nodes[1],
                          bottom.nodes[2],
                          bottom.nodes[3],
                          tip,
                          bottom.nodes[4],
                          bottom.nodes[5],
                          bottom.nodes[6],
                          sideedges[1].nodes[3],
                          sideedges[2].nodes[3],
                          sideedges[3].nodes[3]
                         ]
            elem.shape = TET10
        elseif shape==WED6 && nunodes==5 # => PYR5
            bottom = filter(f -> length(f.nodes)==4, validfaces)[1]
            reverse!(bottom.nodes)
            tip = setdiff(uniquenodes, bottom.nodes)[1]

            elem.nodes = [
                          bottom.nodes[1],
                          bottom.nodes[2],
                          bottom.nodes[3],
                          bottom.nodes[4],
                          tip
                         ]
            elem.shape = PYR5
        elseif shape==WED15 && nunodes==13 # => PYR13
            bottom = filter(f -> length(f.nodes)==8, validfaces)[1]
            bottom.nodes = bottom.nodes[[1,4,3,2,8,7,6,5]]
            edges = get_edges(elem)

            sideedges = Edge[]
            for node in bottom.nodes[1:4]
                for edge in edges
                    n1, n2 = edge.nodes
                    if node==n1 && !(n2 in bottom.nodes)
                        push!(sideedges, edge)
                        break
                    end
                    if node==n2 && !(n1 in bottom.nodes)
                        edge.nodes[1], edge.nodes[2] = edge.nodes[2], edge.nodes[1]
                        push!(sideedges, edge)
                        break
                    end
                end
            end

            tip = sideedges[1].nodes[2]

            elem.nodes = [
                          bottom.nodes[1],
                          bottom.nodes[2],
                          bottom.nodes[3],
                          bottom.nodes[4],
                          tip,
                          bottom.nodes[5],
                          bottom.nodes[6],
                          bottom.nodes[7],
                          bottom.nodes[8],
                          sideedges[1].nodes[3],
                          sideedges[2].nodes[3],
                          sideedges[3].nodes[3],
                          sideedges[4].nodes[3]
                         ]
            elem.shape = PYR13
        else
            collapse_error(shape, nunodes)
        end
    else
        error("collapse: cannot collapse $(shape.name) element.")
    end

    return elem

end
