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
    path::Path
    embedded ::Bool
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
        path::Path;
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


        if path.closed || embedded
            tipjoint = :none
        end
        this = new(path, embedded, cellshape, tag, jointtag, tipjointtag, tipjoint, tol, toln, tolc, lam, id)
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

const PathInset = BlockInset
export PathInset


function Base.copy(bl::BlockInset; dx=0.0, dy=0.0, dz=0.0)
    BlockInset(getcoords(bl.nodes) .+ [dx dy dz], curvetype=bl.curvetype, closed=bl.closed,
                       embedded=bl.embedded, cellshape=bl.cellshape, tag=bl.tag,
                       jointtag=bl.jointtag, tipjointtag=bl.tipjointtag, tipjoint=bl.tipjoint)
end

function Base.copy(blocks::Array{BlockInset,1}; dx=0.0, dy=0.0, dz=0.0)
    return [ copy(obj; dx=dx, dy=dy, dz=dz) for obj in blocks ]
end


function split_block(bl::BlockInset, mesh::Mesh)
    ncells = length(mesh.elems) # initial number of cells in mesh
    pointdict = Dict{UInt64, Node}()
    joints = Array{Cell,1}()
    
    for cmd in bl.path.cmds
        p = endpoint(cmd)
        node = Node(p.coord)
        push!(mesh.nodes, node)
        pointdict[hash(node)] = node
        cmd isa MoveCmd && continue
        
        split!(cmd, mesh, bl.cellshape, pointdict, joints)
    end

    newjoints = mesh.elems[ncells+1:end].linejoints

    if bl.tipjoint in (:front, :both)
        joint = newjoints[1]
        tip = joint.linked_elems[2].nodes[1]
        tipjointnodes = vcat(joint.linked_elems[1].nodes, tip)
        tipjointcell = Cell(tipjoint, tipjointnodes, tag=bl.tipjointtag)
        tipjointcell.linked_elems = joint.linked_elems
        push!(mesh.elems, tipjointcell)
    end
    if bl.tipjoint in (:end, :both)
        joint = newjoints[end]
        tip = joint.linked_elems[2].nodes[2]
        tipjointnodes = vcat(joint.linked_elems[1].nodes, tip )
        tipjointcell = Cell(tipjoint, tipjointnodes, tag=bl.tipjointtag)
        tipjointcell.linked_elems = joint.linked_elems
        push!(mesh.elems, tipjointcell)
    end

    bl._startpoint = nothing
    bl._endpoint = nothing
end


# function get_point(s::Float64, coords::Array{Float64,2}, curvetype::Int)
#     s = s>1.0 ? 1.0 : s
#     if curvetype<=2;
#         return interLagrange(s, coords)
#     else
#         return cubicBezier(s, coords)
#     end
# end


function split_path(cmd::PathCmd,  msh::Mesh, cellshape, pointdict, joints)
    # Tolerances
    ε  = bl.ε
    εn = bl.εn
    εc = bl.εc
    λ  = bl.λ

    # itcount=0 ##

    # Constants
    shape   = cellshape
    npoints = shape==LIN2 ? 2 : 3
    jntshape = shape==LIN2 ? JLINK2 : JLINK3

    # Initial conditions
    len   = 1.0

    # Defining required vectors
    X1 = cmd.p1.coord

    # Find the initial and final element
    ecells = Cell[]
    s0 = cmd(εn)
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
            X = cmd(st)
            is_in = is_inside(ccell.shape, ccell_coords, X, ε)
            if !is_in
                step  = 0.5*(st-s)
                break
            end
        end

        s += step
        X  = cmd(s)
        n  = floor(Int, log(2, step/ε)) + 1  # number of required iterations to find intersection


        for i in 1:n
            is_in = is_inside(ccell.shape, ccell_coords, X, εc)
            step *= 0.5
            if is_in
                s += step
            else
                s -= step
            end

            X = cmd(s)
        end

        # Check if end was reached
        if s > len - εn
            end_reached = true
        end

        # Getting line cell nodes
        P1 = Node(X1)
        P1 = get(pointdict, hash(P1), P1)
        P2 = Node(X)
        end_reached && (P2 = get(pointdict, hash(P2), P2)) # makes sure it picks the first point in case of closed
        
        if npoints==2
            Ps = [P1, P2]
        else
            P3 = Node(cmd((sp+s)/2))
            push!(msh.nodes, P3)
            Ps = [P1, P2, P3]
        end

        # Saving line cell
        lcell = Cell(shape, Ps, tag=bl.tag)
        push!(msh.elems, lcell)

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
            push!(joints, jntcell)
        end

        ccell.crossed = true
        end_reached && return

        # Preparing for the next iteration
        ncell  = find_elem(cmd(s + εn), msh.elems, msh._elempartition, εc, exclude=[ccell])
        ncell === nothing && error("Hole found while searching for next tresspassed cell")

        ccell = ncell
        sp = s
        s = s+εn
    end
end
