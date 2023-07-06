# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DrainPipe

struct DrainPipeProps<:ElemProperties
    A::Float64

    function DrainPipeProps(;A=0.0, diameter=0.0)
        @check A>0.0 || diameter>0.0
        A>0 || (A=π*diameter^2/4)
        return new(A)
    end    
end

DrainPipe = DrainPipeProps

mutable struct DrainPipeElem<:HydromechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    props ::DrainPipeProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function DrainPipeElem(props=DrainPipeProps())
        this = new()
        this.props = props
        return this
    end
end

matching_shape_family(::Type{DrainPipeElem}) = LINECELL
matching_elem_type(::Type{DrainPipeProps}) = DrainPipeElem
matching_props_type(::Type{DrainPipeElem}) = DrainPipeProps

function elem_config_dofs(elem::DrainPipeElem)
    for node in elem.nodes
        add_dof(node, :uw, :fw)
    end
end

function elem_conductivity_matrix(elem::DrainPipeElem)
local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    H  = zeros(nnodes, nnodes)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # mount Bw
        Bw .= 0.0
        Bw = dNdR'/detJ

        # compute H
        coef = detJ*ip.w*(elem.matparams.k/elem.env.anaprops.γw)*A
        @gemm H -= coef*Bw'*Bw
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map, elem.nodes
end

function elem_RHS_vector(elem::DrainPipeElem)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    Q  = zeros(nnodes)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        Jvert = J[end]/detJ

        # mount Bw
        Bw = dNdR'/detJ

        # compute Q 
        coef = detJ*ip.w*elem.matparams.k*A*Jvert
        Q += coef*Bw'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end

function elem_internal_forces(elem::DrainPipeElem, F::Array{Float64,1})
    local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)
    dFw = zeros(nnodes)

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # mount Bw
        Bw = dNdR'/detJ

        # internal volumes dFw
        D    = ip.state.D
        coef = detJ*ip.w*A
        dFw += coef*Bw'*D

    end

    F[map_w] += dFw
end

function update_elem!(elem::DrainPipeElem, DU::Array{Float64,1}, Δt::Float64)
    local A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    Bw = zeros(1, nnodes)

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ]
    Uw += dUw # nodal pore-pressure at step n+1

    J  = Array{Float64}(undef, ndim, 1)
    dFw = zeros(nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        Jvert = J[end]/detJ

        # mount Bw
        Bw = dNdR'/detJ

        # compute Nw vector
        N  = elem.shape.func(ip.R)
        Nw = N'

        # flow gradient
        G  = dot(Bw,Uw)/(elem.env.anaprops.γw) # flow gradient
        G += Jvert; # gradient due to gravity
        Δuw = dot(Nw,dUw) # interpolation to the integ. point

        V = update_state!(elem.matparams, ip.state, Δuw, G, Δt)

        coef = Δt*detJ*ip.w*A
        dFw += coef*Bw'*V
    end

    return dFw, map_w, success()
end


function elem_vals(elem::DrainPipeElem)
    # get area and average fluid velocity and flow
    vals = OrderedDict(:A => elem.props.A )
    mean_va = mean( ip_state_vals(elem.matparams, ip.state)[:va] for ip in elem.ips )
    vals[:va] = mean_va
    vals[:Qa] = elem.props.A*mean_va
    return vals
end
