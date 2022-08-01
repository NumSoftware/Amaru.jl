# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
`BlockInset(coords, [curvetype=0,] [cellshape=LIN3,] [closed=false,] [tag="",] [toln=1e-4,] [tolc=1e-9,] [lam=1.0,])`

Generates an inset block object from a matrix of coordinates `coords`. `curvetype` is an integer value
that represents the inset curve and can be:

0:polyline, 2:lagrangian, 3:cubic Bezier.

`cellshape` represents the shape for 1D elements used in the final mesh. Possible values are LIN2 and LIN3.
`closed=true` can be used if a closed inset curve is required.
"""
mutable struct BlockInset <: AbstractBlock
    nodes    ::Array{Node,1}
    curvetype::Union{Int,AbstractString} # 0:polyline, 1:closed polyline, 2: lagrangian, 3:cubic Bezier with inner points
    closed   ::Bool
    embedded ::Bool
    shape    ::CellShape
    cellshape::CellShape
    tag      ::String
    jointtag ::String
    tipjointtag::String
    tipjoint::Symbol
    ε        ::Float64 # bisection tolerance
    εn       ::Float64 # increment to find next cell
    εc       ::Float64 # tolerance to find cells
    λ        ::Float64 # jump distance to find multiple intersections in one cell
    id       ::Int64
    icount   ::Int64
    _endpoint  ::Union{Node, Nothing}
    _startpoint::Union{Node, Nothing}

    function BlockInset(
        coords::Array{<:Real,2}; 
        curvetype = 0,
        closed::Bool    = false,
        embedded::Bool  = false,
        cellshape = LIN3,
        tag     ::String  = "",
        jointtag::String  = "",
        tipjointtag::String  = "",
        tipjoint::Symbol = :none,
        tol     ::Float64 = 1e-9,
        toln    ::Float64 = 1e-4,
        tolc    ::Float64 = 1e-9,
        lam     ::Float64 = 1.0,
        id      ::Int     = -1,
    )
        # TODO: add option: merge_points
        # TODO: add case: endpoints outside mesh
        if typeof(curvetype)<:Integer
            if !(0<=curvetype<=3); error("Wrong curve type") end
            ctype = curvetype
        else
            cases = Dict("polyline"=>0, "closed polyline"=>1, "lagrangian"=>2, "Bezier"=>3, "bezier"=>3)
            ctype = get(cases, curvetype, -1)
            if ctype==-1; error("Wrong curve type") end
        end

        nrows  = size(coords,1)
        nodes = [ Node(coords[i,:]) for i in 1:nrows ]

        closed && (tipjoint=:none)
        embedded && (tipjoint=:none)
        this = new(nodes, ctype, closed, embedded, LIN2, cellshape, tag, jointtag, tipjointtag, tipjoint, tol, toln, tolc, lam, id)
        this.icount = 0
        this.ε  = tol
        this.εn = toln
        this.εc = tolc
        this.λ  = lam
        this._endpoint   = nothing
        this._startpoint = nothing
        return this
    end
end


function Base.copy(bl::BlockInset; dx=0.0, dy=0.0, dz=0.0)
    BlockInset(getcoords(bl.nodes) .+ [dx dy dz], curvetype=bl.curvetype, closed=bl.closed,
                       embedded=bl.embedded, cellshape=bl.cellshape, tag=bl.tag,
                       jointtag=bl.jointtag, tipjointtag=bl.tipjointtag, tipjoint=bl.tipjoint)
end


function cubicBezier(s::Float64, coords::Array{Float64,2}, isQ::Bool=false)
    # check
    if size(coords,1) < 4
        println("coords = ", coords)
        error("List of points must have at least 4 points for cubic Bezier.")
    end

    # input data
    ndim = size(coords,2)     # space dimension
    np   = size(coords,1)     # number of points
    ns   = np - 1             # number of spacings
    nb   = floor(Int64, ns/3) # number of bezier curves
    ds   = 1.0 / nb           # spacing between curves

    # find index of Bezier and local coordinate t
    ib = floor(Int64, s/ds) + 1    # index of Bezier
    if ib > nb; ib = nb end # fix index if s ~= 1+eps
    s0 = (ib-1) * ds        # s @ left point
    t  = (s - s0) / ds      # local t for current Bezier
    if t > 1.0; t = 1.0 end # clean rubbish. e.g. 1.00000000002

    # collect control points
    Q = zeros(4, ndim)    # control points
    k = 1 + (ib-1) * 3    # position of first point of bezier
    if isQ
        for i in 1:4
            Q[i,:] = coords[k+i-1]
        end
    else
        PQ1 = coords[k  ,:]
        PQ2 = coords[k+1,:]
        PQ3 = coords[k+2,:]
        PQ4 = coords[k+3,:]
        Q[1,:] =         PQ1
        Q[2,:] = (-5.0 * PQ1 + 18.0 * PQ2 -  9.0 * PQ3 + 2.0 * PQ4) / 6.0
        Q[3,:] = ( 2.0 * PQ1 -  9.0 * PQ2 + 18.0 * PQ3 - 5.0 * PQ4) / 6.0
        Q[4,:] =                                               PQ4
    end

    # compute Bezier
    Q1, Q2, Q3, Q4 = Q[1,:], Q[2,:], Q[3,:], Q[4,:]
    a =       Q4 - 3.0 * Q3 + 3.0 * Q2 - Q1
    b = 3.0 * Q3 - 6.0 * Q2 + 3.0 * Q1
    c = 3.0 * Q2 - 3.0 * Q1
    d =       Q1
    return vec(a*t*t*t + b*t*t + c*t + d)
end

function interLagrange(s::Float64, coords::Array{Float64,2})
    #  Interpolates coordinates for s between 0 and 1
    #
    #  0                   +1
    #  0---1---2---3---..---n  -->s
    n = size(coords, 1)
    S = zeros(n)

    for i in 1:n
        f  = 1.0
        si = 1.0*(i-1)/(n-1.0)
        for j in 1:n
            sj = 1.0*(j-1)/(n-1.0)
            if i != j
                f *= (s - sj)/(si - sj)
            end
        end
        S[i] = f
    end
    return coords'*S
end


function split_block(bl::BlockInset, mesh::Mesh)
    coords = getcoords(bl.nodes)
    n, ndim = size(coords)

    if n<2; error("At list two points are required in BlockInset") end
    # 0:polyline, 1:closed polyline, 2: lagrangian, 3:cubic Bezier, 4:Bezier with control points

    ncells = length(mesh.elems) # initial number of cells in mesh

    if bl.curvetype in (2,3) # Lagrangian or Bezier with inner points
        split_curve(coords, bl, false, mesh)
    else # Polyline
        if bl.closed && n<3; error("At least three points are required for closed polyline in BlockInset") end

        bl._endpoint = nothing
        for i in 1:n-1
            coordsi = coords[i:i+1,:]
            split_curve(coordsi, bl, false, mesh)
        end
        if bl.closed
            coordsi = [ coords[n:n,:] ; coords[1:1,:] ]
            split_curve(coordsi, bl, true, mesh)
        end
    end

    newjoints = mesh.elems[ncells+1:end].linejoints

    if bl.tipjoint in (:front, :both)
        joint = newjoints[1]
        tip = joint.linked_elems[2].nodes[1]
        tipjointpts  = vcat(joint.linked_elems[1].nodes, tip)
        tipjointcell = Cell(TIPJOINT, tipjointpts, tag=bl.tipjointtag)
        tipjointcell.linked_elems = joint.linked_elems
        push!(mesh.elems, tipjointcell)
    end
    if bl.tipjoint in (:end, :both)
        joint = newjoints[end]
        tip = joint.linked_elems[2].nodes[2]
        tipjointpts  = vcat(joint.linked_elems[1].nodes, tip )
        tipjointcell = Cell(TIPJOINT, tipjointpts, tag=bl.tipjointtag)
        tipjointcell.linked_elems = joint.linked_elems
        push!(mesh.elems, tipjointcell)
    end

    bl._startpoint = nothing
    bl._endpoint = nothing
end

function get_point(s::Float64, coords::Array{Float64,2}, curvetype::Int)
    s = s>1.0 ? 1.0 : s
    if curvetype<=2;
        return interLagrange(s, coords)
    else
        return cubicBezier(s, coords)
    end
end

function split_curve(coords::Array{Float64,2}, bl::BlockInset, closed::Bool, msh::Mesh)
    # Tolerances
    ε  = bl.ε
    εn = bl.εn
    εc = bl.εc
    λ  = bl.λ

    itcount=0 ##

    # Constants
    shape   = bl.cellshape
    npoints = shape==LIN2 ? 2 : 3
    jntshape = shape==LIN2 ? JLINK2 : JLINK3
    curvetype = bl.curvetype

    # Initial conditions
    bdist = 0.0  # boundary function initial value
    len   = 1.0

    # Defining required vectors
    X1 = vec(coords[  1,:])
    Xn = vec(coords[end,:])

    # Find the initial and final element
    ecells = Cell[]
    s0 = get_point(εn, coords, curvetype)
    icell = find_elem(s0, msh.elems, msh._elempartition, εc, exclude=ecells) # The first tresspased cell

    if icell === nothing
        error("Inset point $(s0) outside the mesh")
    end

    # Initializing more variables
    ccell  = icell
    # nodes = Array{Node}(undef, npoints)

    # Do not set _endpoint to nothing ( bl._endpoint = nothing ) to allow connectivity between segments!

    end_reached  = false
    s  = 0.0
    sp = 0.0
    nits = round(Int, 1.0/λ)

    # Splitting inset
    k = 0
    while true
        k +=1
        ccell_coords = getcoords(ccell)
        # Default step
        step  = 0.50*(1.0-s)

        # Finding step
        st = s     # trial point
        for i in 1:nits
            st += λ
            if st>1.0; break end
            X =get_point(st, coords, curvetype)
            is_in = is_inside(ccell.shape, ccell_coords, X, ε)
            if !is_in
                step  = 0.5*(st-s)
                break
            end
        end

        s += step
        X  = get_point(s, coords, curvetype)
        n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection

        itcount+=n ##

        for i in 1:n

            step *= 0.5
            is_in = is_inside(ccell.shape, ccell_coords, X, εc)
            if is_in
                s += step
            else
                s -= step
            end

            X =get_point(s, coords, curvetype)

            R     = inverse_map(ccell.shape, ccell_coords, X)
            bdist = bdistance(ccell.shape, R)
        end

        # Check if end was reached
        if s > len - εn
            end_reached = true
            if curvetype<=1; X=Xn end
            # TODO test (check also incomplete Bezier...)
        end

        # Counter
        bl.icount += end_reached ? 0 : 1

        # Getting line cell nodes
        if bl._endpoint===nothing
            P1 = Node(X1)
            push!(msh.nodes, P1)
        else
            P1 = bl._endpoint
        end

        if !(closed && end_reached)
            P2 = Node(X)
            push!(msh.nodes, P2)
        else
            P2 = bl._startpoint
        end

        if npoints==2
            Ps = [P1, P2]
        else
            P3 = Node(get_point( (sp+s)/2.0, coords, curvetype))
            push!(msh.nodes, P3)
            Ps = [P1, P2, P3]
        end

        # Saving line cell
        lcell = Cell(shape, Ps, tag=bl.tag)
        push!(msh.elems, lcell)

        if bl._startpoint === nothing; 
            bl._startpoint = P1 
        end
        bl._endpoint = P2

        if bl.embedded
            # Set line as embedded
            lcell.embedded = true
            lcell.linked_elems = [ ccell ]
        else
            # Generate a continuous joint element
            jntpts  = vcat( ccell.nodes, lcell.nodes )
            jntcell = Cell(jntshape, jntpts, tag=bl.jointtag)
            push!(msh.elems, jntcell)
            jntcell.linked_elems = [ccell, lcell]
        end

        ccell.crossed = true

        if end_reached
            return
        end

        # Preparing for the next iteration
        ncell  = find_elem(get_point(s + εn, coords, curvetype), msh.elems, msh._elempartition, εc, exclude=[ccell])
        if ncell === nothing
            error("Hole found while searching for next tresspassed cell")
        end

        ccell = ncell
        sp = s
        s = s+εn
    end
end
