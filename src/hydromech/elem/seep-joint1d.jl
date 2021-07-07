# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct SeepJoint1D<:Hydromechanical
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    # specific fields
    uw_dofs   ::Array{Dof,1} # list of pore-pressure dofs
    cache_B   ::Array{Array{Float64,2}}
    cache_detJ::Array{Float64}

    SeepJoint1D() = new()
end

matching_shape_family(::Type{SeepJoint1D}) = LINEJOINT_SHAPE

function elem_config_dofs(elem::SeepJoint1D)
    # No need to add dofs since they will be added by linked elements
end

function elem_init(elem::SeepJoint1D)
    # Only pore-pressure dofs
    elem.uw_dofs = get_dofs(elem, :uw)

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

function mountB(elem::SeepJoint1D, R, Ch, Ct)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    nnodes  = length(elem.nodes)
    nbnodes = length(bar.nodes)
    D = bar.shape.deriv(R)
    J = D*Ct

    # Mount NN matrix
    N = bar.shape.func(R)
    NN = transpose(N)

    # Babuška-Brezzi condition
    hookshape = hook.shape
    hook_uw_dofs = get_dofs(hook, :uw)
    if length(hook.nodes)!=length(hook_uw_dofs) 
        hookshape = hook.shape.basic_shape
    end

    # Mount MM matrix
    stack = Array{Float64,1}[]
    for i=1:nbnodes
        Xj = bar.nodes[i].coord
        R  = inverse_map(hook.shape, Ch, Xj)
        M  = hookshape.func(R)
        push!(stack, M)
    end

    MM = transpose(hcat(stack...))
    B = [ NN*MM  -NN ]

    detJ = norm(J)

    return B, detJ
end

function elem_conductivity_matrix(elem::SeepJoint1D)
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    h    = elem.mat.h
    map_w = [ dof.eq_id for dof in elem.uw_dofs ]
    nuwnodes = length(map_w)

    H  = zeros(nuwnodes, nuwnodes)

    for (i,ip) in enumerate(elem.ips)
        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        coef = detJ*ip.w*(elem.mat.k/elem.mat.γw)*h
        H -= coef*Bp'*Bp
    end

    return H, map_w, map_w, elem.nodes
end

function elem_internal_forces(elem::SeepJoint1D, F::Array{Float64,1})
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    h    = elem.mat.h

    map_w = [ dof.eq_id for dof in elem.uw_dofs ]
    nuwnodes = length(map_w)
    dFw  = zeros(nuwnodes)

    for (i,ip) in enumerate(elem.ips)
        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        # internal volumes dFw
        D    = ip.state.D
        coef = detJ*ip.w*h
        dFw += coef*Bp'*D
    end

    F[map_w] += dFw
end

function elem_update!(elem::SeepJoint1D, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    h    = elem.mat.h

    map_w = [ dof.eq_id for dof in elem.uw_dofs ]
    nuwnodes = length(map_w)

    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ dof.vals[:uw] for dof in elem.uw_dofs ]
    Uw += dUw # nodal pore-pressure at step n+1

    dFw = zeros(nuwnodes)

    for (i,ip) in enumerate(elem.ips)

        Bp   = elem.cache_B[i]
        detJ = elem.cache_detJ[i]

        # poropression difference between solid and drain
        ΔFw = dot(Bp,Uw)/elem.mat.γw
        V = update_state!(elem.mat, ip.state, ΔFw, Δt)

        coef = Δt*detJ*ip.w*h
        dFw += coef*Bp'*V
    end

    DF[map_w] += dFw
    return success()
end

function elem_extrapolated_node_vals(elem::SeepJoint1D)
    all_ip_vals = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
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
