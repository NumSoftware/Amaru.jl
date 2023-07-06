# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechTipJointElem<:MechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    MechTipJoint() = new()
end

matching_shape_family(::Type{MechTipJointElem}) = TIPJOINTCELL

function mountB(elem::MechTipJointElem, Ch, Ct)
    # Calculates the matrix that relates nodal displacements with relative displacements
    # ==================================================================================

    # B = T* [-MM'  I]      ndim x ndim*(m+n)

    # where
    # T is a untari vector pointing outwards the rod
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the tip node coords

    # MM' = [ M_1*I M_2*I ... M_m*I]
    # I is a ndim x ndim identity matrix


    ndim = elem.env.ndim
    tip = elem.nodes[end]
    bulk = elem.linked_elems[1]
    rod  = elem.linked_elems[2]
    nsnodes = length(bulk.nodes)

    if hash(tip)==hash(rod.nodes[1])
        R = [-1.0, 0, 0]
        D = rod.shape.deriv(R)
        J = Ct'*D
        T = -normalize(J)
    else
        R = [+1.0, 0, 0]
        D = rod.shape.deriv(R)
        J = Ct'*D
        T = normalize(J)
    end

    # Mount MM matrix
    stack = Array{Float64,2}[]
    X = tip.coord
    R = inverse_map(bulk.shape, Ch, X)
    M = bulk.shape.func(R)
    for Mi in M
        push!(stack, Mi*eye(ndim))
    end

    MM = hvcat(nsnodes, stack...)
    B = T'*[ -MM  eye(ndim) ]
    return B
end

function elem_stiffness(elem::MechTipJointElem)
    ndim = elem.env.ndim
    bulk = elem.linked_elems[1]
    rod  = elem.linked_elems[2]
    Ch = getcoords(bulk)
    Ct = getcoords(rod)

    B = mountB(elem, Ch, Ct)
    k = calcD(elem.matparams, elem.ips[1].state)
    coef = k
    K = coef*B'*B

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end

function update_elem!(elem::MechTipJointElem, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    bulk = elem.linked_elems[1]
    rod  = elem.linked_elems[2]
    Ch = getcoords(bulk)
    Ct = getcoords(rod)

    # nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    B  = mountB(elem, Ch, Ct)
    Δw = dot(B, dU)
    
    Δf, _ = update_state(elem.matparams, elem.ips[1].state, Δw)
    coef = Δf
    dF = coef*B'

     return dF, map, success()
end
