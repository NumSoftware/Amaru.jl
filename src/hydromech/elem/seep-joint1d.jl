# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct SeepJoint1D<:Hydromechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    # specific fields
    cache_B   ::Array{Array{Float64,2}}
    cache_detJ::Array{Float64}

    SeepJoint1D() = new()
end

matching_shape_family(::Type{SeepJoint1D}) = JOINT1D_SHAPE

function elem_config_dofs(elem::SeepJoint1D)
    nnodes = length(elem.nodes)
    for (i, node) in enumerate(elem.nodes)
        add_dof(node, :uw, :fw)
    end
end

function elem_init(elem::SeepJoint1D)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    Ch = elem_coords(hook)
    Ct = elem_coords(bar)
    elem.cache_B = []
    elem.cache_detJ = []
    for ip in elem.ips
        B, detJ = mountB(elem, ip.R, Ch, Ct)
        push!(elem.cache_B, B)
        push!(elem.cache_detJ, detJ)
    end
    return nothing
end

function mountB(elem::SeepJoint1D, R, Ch, Ct)
    ndim = elem.env.ndim
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    nnodes  = length(elem.nodes)
    nbnodes = length(bar.nodes)
    nsnodes = length(hook.nodes)
    D = bar.shape.deriv(R)
    J = D*Ct
    @show J 
    @show D 
    @show Ct

    # Mount NN matrix
    N = bar.shape.func(R)
    NN = hcat([ Ni for Ni in N  ]...)

    # Mount MM matrix
    stack = Array{Float64,1}[]
    for i=1:nbnodes
        Xj = bar.nodes[i].X
        R  = inverse_map(hook.shape, Ch, Xj)
        M  = hook.shape.func(R)
        for Mi in M
            push!(stack, Mi*ones(1))
        end
    end

    MM = hvcat(nsnodes, stack...)
    B = [ NN*MM  -NN ]

    detJ = norm(J)
    @show detJ
    return B, detJ
end

function elem_conductivity_matrix(elem::SeepJoint1D)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    h    = elem.mat.h

    H  = zeros(nnodes, nnodes)

    for (i,ip) in enumerate(elem.ips)
        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        coef = detJ*ip.w*(elem.mat.k/elem.mat.γw)*h
        H -= coef*Bp'*Bp 
    end

	map_p = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map_p, map_p
end

function elem_internal_forces(elem::SeepJoint1D, F::Array{Float64,1})
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    dFw  = zeros(nnodes)
    h    = elem.mat.h

    map_p = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    for (i,ip) in enumerate(elem.ips)
        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        # internal volumes dFw    
        V    = ip.data.V
        coef = detJ*ip.w*h
        dFw += coef*Bp'*V
    end

    F[map_p] += dFw
end

function elem_update!(elem::SeepJoint1D, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    h    = elem.mat.h

    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]
    
    dUw = DU[map_p] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ]
    Uw += dUw # nodal pore-pressure at step n+1

    dFw = zeros(nnodes)

    for (i,ip) in enumerate(elem.ips)

        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        # flow gradient
        G  = dot(Bp,Uw)/elem.mat.γw # flow gradient
        Δuw = dot(Bp,dUw)# interpolation to the integ. point

        V = update_state!(elem.mat, ip.data, Δuw, G)

        coef = Δt*detJ*ip.w*h
        dFw += coef*Bp'*V
    end

    DF[map_p] += dFw
end

function elem_extrapolated_node_vals(elem::SeepJoint1D)
    all_ip_vals = [ ip_state_vals(elem.mat, ip.data) for ip in elem.ips ]
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