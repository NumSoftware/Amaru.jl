# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


mutable struct MechEmbRodElem<:MechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    props ::MechRodProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    # specific fields
    cache_NN::Array{Float64,2}

    function MechEmbRodElem(props=MechRodProps())
        this = new()
        this.props = props
        return this
    end
end

matching_shape_family(::Type{MechEmbRodElem}) = LINECELL



function elem_config_dofs(::MechEmbRodElem)
    # No-op function.
    # The solid linked element will set the required dofs.
end


function elem_map(elem::MechEmbRodElem)
    ndim = elem.env.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    solid = elem.linked_elems[1]
    return [ node.dofdict[key].eq_id for node in solid.nodes for key in keys ]
end


function _mountNN(elem::MechEmbRodElem)
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


function elem_init(elem::MechEmbRodElem)
    elem.cache_NN = _mountNN(elem)
    return nothing
end


function elem_displacements(elem::MechEmbRodElem)
    ndim    = elem.env.ndim
    NN      = elem.cache_NN
    keys    = (:ux, :uy, :uz)[1:ndim]
    Ubulk   = [ node.dofdict[key].vals[key] for node in elem.linked_elems[1].nodes for key in keys ]
    Urod    = NN'*Ubulk
    nodemap = [ node.id for node in elem.nodes for key in keys ]
    dimmap  = [ i for node in elem.nodes for i in 1:ndim ]

    return Urod, nodemap, dimmap
end


function elem_stiffness(elem::MechEmbRodElem)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A = elem.props.A
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

        E    = calcD(elem.matparams, ip.state)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end

    NN = elem.cache_NN
    map = elem_map(elem)
    return NN*K*NN', map, map
end


function update_elem!(elem::MechEmbRodElem, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A      = elem.props.A
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
        dsig, _ = update_state(elem.matparams, ip.state, deps)
        coef = A*detJ*ip.w
        dF  += coef*B'*dsig
    end

    return NN*dF, map, success()
end


function elem_vals(elem::MechEmbRodElem)
    # get area and average stress and axial force
    vals = OrderedDict(:A => elem.props.A )
    mean_sa = mean( ip_state_vals(elem.matparams, ip.state)[:sa] for ip in elem.ips )
    vals[:sa] = mean_sa
    vals[:fa] = elem.props.A*mean_sa
    return vals
end
