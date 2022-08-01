# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechEmbRod<:Mechanical
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
    cache_NN::Array{Float64,2}

    MechEmbRod() = new()
end

matching_shape_family(::Type{MechEmbRod}) = LINE_CELL

function elem_config_dofs(elem::MechEmbRod)
    # No-op function.
    # The solid linked element will set the required dofs.
end

function elem_map(elem::MechEmbRod)
    ndim = elem.env.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    solid = elem.linked_elems[1]
    return [ node.dofdict[key].eq_id for node in solid.nodes for key in keys ]
end

function _mountNN(elem::MechEmbRod)
    ndim = elem.env.ndim
    solid = elem.linked_elems[1]
    n  = length(solid.nodes)
    m  = length(elem.nodes)
    NN = zeros(ndim*n, ndim*m)
    Cs = getcoords(solid)

    for j in 1:m
        R = inverse_map(solid.shape, Cs, elem.nodes[j].coord)
        N = solid.shape.func(R)
        for i in 1:n
            for k in 1:ndim
                NN[(i-1)*ndim+k, (j-1)*ndim+k] = N[i]
            end
        end
    end
    return NN
end


function elem_init(elem::MechEmbRod)
    elem.cache_NN = _mountNN(elem)
    return nothing
end

function elem_displacements(elem::MechEmbRod)
    ndim   = elem.env.ndim
    NN     = elem.cache_NN
    keys   = (:ux, :uy, :uz)[1:ndim]
    Ubulk  = [ node.dofdict[key].vals[key] for node in elem.linked_elems[1].nodes for key in keys ]
    nodemap = [ node.id for node in elem.nodes for key in keys ]
    dimmap  = [ i for node in elem.nodes for i in 1:ndim ]

    return NN'*Ubulk, nodemap, dimmap
end

function elem_stiffness(elem::MechEmbRod)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A = elem.mat.A
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        E    = calcD(elem.mat, ip.state)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end

    NN = elem.cache_NN
    map = elem_map(elem)
    return NN*K*NN', map, map
end

function elem_update!(elem::MechEmbRod, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A      = elem.mat.A
    NN     = elem.cache_NN


    map = elem_map(elem)
    dU  = U[map]
    dUr = NN'*dU

    dF = zeros(nnodes*ndim)
    C  = getcoords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, 1)
    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        deps = (B*dUr)[1]
        dsig, _ = stress_update(elem.mat, ip.state, deps)
        coef = A*detJ*ip.w
        dF  += coef*B'*dsig
    end

    F[map] += NN*dF
    return success()
end


function elem_vals(elem::MechEmbRod)
    # get area and average stress and axial force
    vals = OrderedDict(:A => elem.mat.A )
    mean_sa = mean( ip_state_vals(elem.mat, ip.state)[:sa] for ip in elem.ips )
    vals[:sa] = mean_sa
    vals[:fa] = elem.mat.A*mean_sa
    return vals
end
