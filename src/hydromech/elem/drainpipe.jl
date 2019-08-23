# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct DrainPipe<:Hydromechanical
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

    function DrainPipe()
        return new() 
    end
end

matching_shape_family(::Type{DrainPipe}) = LINE_SHAPE

function elem_config_dofs(elem::DrainPipe)
    nnodes = length(elem.nodes)
    for (i, node) in enumerate(elem.nodes)
        add_dof(node, :uw, :fw)
    end
end

function elem_conductivity_matrix(elem::DrainPipe)
local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.mat.A
    C  = elem_coords(elem)
    H  = zeros(nnodes, nnodes)
    Bp = zeros(1, nnodes)
    J  = Array{Float64}(undef, 1, ndim)

    for ip in elem.ips
  
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # mount Bp
        Bp .= 0.0
        for i in 1:nnodes
            Bp[1,i] = dNdR[1,i]*(1/detJ)
        end

        # compute H
        coef = detJ*ip.w*(elem.mat.k/elem.mat.γw)*A
        @gemm H -= coef*Bp'*Bp
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map
end

function elem_compressibility_matrix(elem::DrainPipe)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.mat.A
    C  = elem_coords(elem)
    Cpp = zeros(nnodes, nnodes) 
    J  = Array{Float64}(undef, 1, ndim)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np vector
        N    = elem.shape.func(ip.R)
        Np   = N'   

        # compute Cuu
        coef = detJ*ip.w*elem.mat.β*A
        Cpp  -= coef*Np'*Np
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Cpp, map, map
end

function elem_RHS_vector(elem::DrainPipe)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.mat.A
    C  = elem_coords(elem)
    Q  = zeros(nnodes)
    Bp = zeros(1, nnodes)
    dNdX = Array{Float64}(undef, 1, nnodes)
    J  = Array{Float64}(undef, 1, ndim)
    Juni  = Array{Float64}(undef, 1, ndim)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)
        Juni = J/norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        Jvert = Juni[1,end]

        # mount Bp
        Bp .= 0.0
        for i in 1:nnodes
            Bp[1,i] = dNdR[1,i]*(1/detJ)
        end

        # compute Q
        coef = detJ*ip.w*elem.mat.k*Jvert*A
        Q += coef*Bp'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end

function elem_internal_forces(elem::DrainPipe, F::Array{Float64,1})
local k::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A  = elem.mat.A
    C  = elem_coords(elem)
    H  = zeros(nnodes, nnodes)
    Bp = zeros(1, nnodes)
    dNdX = Array{Float64}(undef, 1, nnodes)
    J  = Array{Float64}(undef, 1, ndim)
    dFw = zeros(nnodes)

    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        
        # mount Bp
        Bp .= 0.0
        for i in 1:nnodes
            Bp[1,i] = dNdR[1,i]*(1/detJ)
        end

        # compute Np vector
        N    = elem.shape.func(ip.R)
        Np   = N'   

        # internal volumes dFw
        uw   = ip.data.uw
        coef = detJ*ip.w*elem.mat.β*A
        dFw -= coef*Np'*uw
        
        V    = ip.data.V
        coef = detJ*ip.w*A
        dFw += coef*Bp'*V

    end

    F[map_p] += dFw
end

function elem_update!(elem::DrainPipe, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    local A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C  = elem_coords(elem)
    A  = elem.mat.A

    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dUw = DU[map_p] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ]
    Uw += dUw # nodal pore-pressure at step n+1

    Bp = zeros(1, nnodes)
    dNdX = Array{Float64}(undef, 1, nnodes)
    J  = Array{Float64}(undef, 1, ndim)
    Juni  = Array{Float64}(undef, 1, ndim)
    dFw = zeros(nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)
        Juni = J/norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        Jvert = Juni[1,end]
        
        # mount Bp
        Bp .= 0.0
        for i in 1:nnodes
            Bp[1,i] = dNdR[1,i]*(1/detJ)
        end

        # compute Np vector
        N    = elem.shape.func(ip.R)
        Np   = N'   

        # flow gradient
        G  = dot(Bp,Uw)/elem.mat.γw # flow gradient
        G += Jvert; # gradient due to gravity
        @show Bp
        @show Uw
        @show G 
        @show Jvert

        Δuw = dot(Np,dUw) # interpolation to the integ. point

        V = update_state!(elem.mat, ip.data, Δuw, G)
        @show Δuw 
        @show V

        coef = A*detJ*ip.w*elem.mat.β
        dFw -= coef*Np'*Δuw  

        coef = Δt*A*detJ*ip.w
        dFw += coef*Bp'*V
    end

    DF[map_p] += dFw
end


function elem_vals(elem::DrainPipe)
    # get area and average fluid velocity and flow
    vals = OrderedDict(:A => elem.mat.A )
    mean_va = mean( ip_state_vals(elem.mat, ip.data)[:va] for ip in elem.ips )
    vals[:va] = mean_va
    vals[:Qa] = elem.mat.A*mean_va
    return vals
end
