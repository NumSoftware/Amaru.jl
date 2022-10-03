"""
    convex_hull(nodes)

Computes the convex hull of an array of `nodes`.
"""
function convex_hull(nodes::Array{Node,1})
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

    @show V1
    V2 = cross(N, V1)
    V3 = cross(V1, V2)

    normalize!(V2)
    normalize!(V3)

    R = [ V1 V2 V3 ]

    Amaru.@showm R

    # C = C*R'
    C = C*R

    # Jarvis march or gift wrapping algorithm

    # find leftmost point
    _, idx1 = findmax(C[:,1])

    perm = Int[]
    while true
        push!(perm, idx1)
        idx2 = 1
        P1 = C[idx1, :]
        P2 = C[idx2, :]
        for j in 1:nnodes
            if idx1==idx2
                idx2 = j
                continue
            end
            
            P2 = C[idx2,:]
            Pj = C[j,:]

            # check if point j is on right of line idx1 to idx2
            if (P2[1] - P1[1])*(Pj[2] - P1[2]) - (P2[2] - P1[2])*(Pj[1] - P1[1]) < 0
                idx2 = j
            end
        end

        idx1 = idx2
        idx2==perm[1] && break
    end

    Amaru.@showm C
    Amaru.@showm C[perm, :]

    # build contour
    return nodes[perm]
end
