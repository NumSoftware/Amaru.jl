# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechRSJoint

struct MechRSJointProps<:ElemProperties
    p::Float64

    function MechRSJointProps(;p=NaN, A=NaN, diameter=NaN)
        # A : section area
        # p : section perimeter
        @check (p>0 || A>0 || diameter>0)

        if isnan(p)
            if A>0
                p = 2.0*(A*pi)^0.5
            else
                p = pi*diameter
            end
        end
        @check p>0

        this = new(p)
        return this
    end
    
end

MechRSJoint = MechRSJointProps


"""
    MechRodSolidJoint

An interface element to link `MechRod` elements to bulk elements
in mechanical equilibrium analyses.
"""
mutable struct MechRSJointElem<:MechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    props ::MechRSJointProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    # specific fields
    cache_B   ::Array{Array{Float64,2}}
    cache_detJ::Array{Float64}

    function MechRSJointElem(props=MechRSJointProps())
        this = new()
        this.props = props
        return this
    end
end

matching_shape_family(::Type{MechRSJointElem}) = LINEJOINTCELL
matching_elem_type(::Type{MechRSJointProps}) = MechRSJointElem
matching_props_type(::Type{MechRSJointElem}) = MechRSJointProps


function elem_init(elem::MechRSJointElem)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    Ch = getcoords(hook)
    Ct = getcoords(bar)
    elem.cache_B = []
    elem.cache_detJ = []
    for ip in elem.ips
        B, detJ = mountB(elem, ip.R, Ch, Ct)
        push!(elem.cache_B, B)
        push!(elem.cache_detJ, detJ)
    end
    return nothing
end

function mount_T(J::Matx)
    ndim = length(J)
    nJ   = norm(J)
    L1   = vec(J/nJ)

    if ndim==2
        L2 = [ -L1[2],  L1[1] ]
        return hcat(L1,L2)'
    end

    # Finding second vector
    if     abs(L1[1]) == 1.0; L2 = [0.0, 1.0, 0.0]
    elseif abs(L1[2]) == 1.0; L2 = [0.0, 0.0, 1.0]
    elseif abs(L1[3]) == 1.0; L2 = [1.0, 0.0, 0.0]
    else
        # Auxiliar vector L which must be different from L1
        L = [1.0, 0.0, 0.0]
        if norm(L-L1) < 1.0e-4; L = [0.0, 1.0, 0.0] end
        # Performing cross product to obtain a second vector
        L2  = normalize(cross(L1, L))
    end

    # Finding third vector
    L3 = normalize(cross(L1, L2))

    return [ L1'; L2'; L3' ]

    #return hcat(L1, L2, L3)'
end


function mountB(elem::MechRSJointElem, R, Ch, Ct)
    # Calculates the matrix that relates nodal displacements with relative displacements
    # ==================================================================================

    # B = T* [NN*MM  -NN]      ndim x ndim*(m+n)

    # where
    # T is a direction cosines matrix
    # NN is a matrix containing truss element shape functions
    # evaluated at the point of interest R.
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the n truss nodal points.

    #          [ M_11*I M_21*I ... M_m1*I]
    #     MM = [ M_12*I M_22*I ... M_m2*I]
    #          [ M_13*I M_23*I ... M_mn*I] ndim*n x ndim*m

    # where
    # M_12 is the first shape function from the tresspased solid element
    # evaluated at the second node of the truss element.
    # I is a ndim x ndim identity matrix


    ndim = elem.env.ndim
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    # nnodes  = length(elem.nodes)
    nbnodes = length(bar.nodes)
    nsnodes = length(hook.nodes)
    D = bar.shape.deriv(R)
    J = Ct'*D
    T = mount_T(J)

    # Mount NN matrix
    N = bar.shape.func(R)
    NN = hcat([ Ni*eye(ndim) for Ni in N  ]...)

    # Mount MM matrix
    stack = Array{Float64,2}[]
    for i in 1:nbnodes
        Xj = bar.nodes[i].coord
        R  = inverse_map(hook.shape, Ch, Xj)
        M  = hook.shape.func(R)
        for Mi in M
            push!(stack, Mi*eye(ndim))
        end
    end
    MM = hvcat(nsnodes, stack...)

    B = T*[ -NN*MM  NN ]
    detJ = norm(J)
    return B, detJ
end

function elem_stiffness(elem::MechRSJointElem)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    p = elem.props.p

    K  = zeros(nnodes*ndim, nnodes*ndim)
    DB = zeros(ndim, nnodes*ndim)

    for (i,ip) in enumerate(elem.ips)
        B    = elem.cache_B[i]
        detJ = elem.cache_detJ[i]
        D    = calcD(elem.matparams, ip.state)
        coef = p*detJ*ip.w
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end

function update_elem!(elem::MechRSJointElem, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    p = elem.props.p

    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    Δu = zeros(ndim)

    # hook = elem.linked_elems[1]
    # bar  = elem.linked_elems[2]
    for (i,ip) in enumerate(elem.ips)
        B    = elem.cache_B[i]
        detJ = elem.cache_detJ[i]
        # D    = calcD(matparams, ip.state)
        @gemv Δu = B*dU
        Δσ, _ = update_state(elem.matparams, ip.state, Δu)
        coef = p*detJ*ip.w
        # Δσ[1]  *= p
        @gemv dF += coef*B'*Δσ
    end

    return dF, map, success()
end


function elem_extrapolated_node_vals(elem::MechRSJointElem)
    all_ip_vals = [ ip_state_vals(elem.matparams, ip.state) for ip in elem.ips ]
    nips        = length(elem.ips)
    fields      = keys(all_ip_vals[1])
    nfields     = length(fields)

    # matrix with all ip values (nip x nvals)
    W = mapreduce(transpose, vcat, collect.(values.(all_ip_vals)))

    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]

    E = extrapolator(bar.shape, nips)
    N = E*W # (nbnodes x nfields)

    nhnodes = length(hook.nodes)
    N = [ zeros(nhnodes, nfields); N ]

    # Filling nodal and elem vals
    node_vals = OrderedDict{Symbol, Array{Float64,1}}(field => N[:,i] for (i,field) in enumerate(fields))

    return node_vals
end
