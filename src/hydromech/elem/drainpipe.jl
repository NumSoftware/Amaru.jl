# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DrainPipe

struct DrainPipeProps<:ElemProperties
    A::Float64

    function DrainPipeProps(; props...)
        names = (A="Section area",)
        required = (:A,)
        @checkmissing props required names

        props  = (; props...)
        A       = props.A

        @check A>0.0
        return new(A)
    end    
end



mutable struct DrainPipe<:Hydromech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    props ::DrainPipeProps
    active::Bool
    linked_elems::Array{Element,1}
    ctx::Context

    function DrainPipe()
        return new()
    end
end

compat_shape_family(::Type{DrainPipe}) = LINECELL
compat_elem_props(::Type{DrainPipe}) = DrainPipeProps



function elem_config_dofs(elem::DrainPipe)
    for node in elem.nodes
        add_dof(node, :uw, :fw)
    end
end


function elem_conductivity_matrix(elem::DrainPipe)
local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    H  = zeros(nnodes, nnodes)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # mount Bw
        Bw .= 0.0
        Bw = dNdR'/detJ

        # compute H
        coef = detJ*ip.w*(elem.mat.k/elem.ctx.γw)*A
        @mul H -= coef*Bw'*Bw
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map, elem.nodes
end


function elem_RHS_vector(elem::DrainPipe)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    Q  = zeros(nnodes)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        Jvert = J[end]/detJ

        # mount Bw
        Bw = dNdR'/detJ

        # compute Q 
        coef = detJ*ip.w*elem.mat.k*A*Jvert
        Q .+= coef*Bw'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end


function elem_internal_forces(elem::DrainPipe, F::Array{Float64,1})
    local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)

    A  = elem.props.A
    C  = getcoords(elem)
    Bw = zeros(1, nnodes)
    J  = Array{Float64}(undef, ndim, 1)
    dFw = zeros(nnodes)

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
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


function update_elem!(elem::DrainPipe, DU::Array{Float64,1}, Δt::Float64)
    local A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.ctx.ndim
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
        @mul J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        Jvert = J[end]/detJ

        # mount Bw
        Bw = dNdR'/detJ

        # compute Nw vector
        N  = elem.shape.func(ip.R)
        Nw = N'

        # flow gradient
        G  = dot(Bw,Uw)/(elem.ctx.γw) # flow gradient
        G += Jvert; # gradient due to gravity
        Δuw = dot(Nw,dUw) # interpolation to the integ. point

        V = update_state!(elem.mat, ip.state, Δuw, G, Δt)

        coef = Δt*detJ*ip.w*A
        dFw += coef*Bw'*V
    end

    return dFw, map_w, success()
end


function elem_vals(elem::DrainPipe)
    # get area and average fluid velocity and flow
    vals = OrderedDict(:A => elem.props.A )
    mean_va = mean( ip_state_vals(elem.mat, ip.state)[:va] for ip in elem.ips )
    vals[:va] = mean_va
    vals[:Qa] = elem.props.A*mean_va
    return vals
end
