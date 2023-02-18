# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct AcusticFluid<:Element
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    env::ModelEnv

    function AcusticFluid();
        return new()
    end
end

matching_shape_family(::Type{AcusticFluid}) = BULKCELL

function elem_config_dofs(elem::AcusticFluid)
    for node in elem.nodes
        add_dof(node, :up, :fp)
    end
end

function elem_init(elem::AcusticFluid)
end


function distributed_bc(elem::AcusticFluid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.thickness
    suitable_keys = (:tq,:ax,:ay,:az) # tq: mass flow acceleration?, ax: x acceleration

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a AcusticFluid element")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = getcoords(nodes, ndim)

    # Calculate the nodal values
    F     = zeros(nnodes)
    J     = Array{Float64}(undef, ndim, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        @gemm J = C'*D
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if elem.env.modeltype=="axisymmetric"
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        end

        if key in (:ax, :ay, :az, :an)
            vip *= -elem.mat.rho
        end

        if  ndim==2
            n = [J[2], -J[1]]
        else
            n = cross(J[:,1], J[:,2])
        end
        normalize!(n)

        
        if key in (:ax)
            Q[1] = vip
        elseif key in (:ay)
            Q[2] = vip
        elseif key in (:az)
            Q[3] = vip
        elseif key in (:an)
            Q = vip*n
        end
        
        vip = dot(n, Q)

        coef = vip*norm(J)*w*th
        F .+= coef*N # F is a vector
    end

    # generate a map
    map  = [ node.dofdict[:up].eq_id for node in target.nodes ]

    return F, map
end


# acustic fluid stiffness
function elem_fluid_stiffness(elem::AcusticFluid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    K      = zeros(nnodes, nnodes)
    dNdX   = zeros(nnodes, ndim) # cartesian derivatives
    J      = Array{Float64}(undef, ndim, ndim) # Jacobian

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = C'*dNdR
        detJ = det(J)
        invJ = inv(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @gemm dNdX = dNdR*invJ # Bp = dNdX'

        coef = detJ*ip.w*th
        @gemm K += coef*dNdX'*dNdX
    end

    # map
    map = [ node.dofdict[:up].eq_id for node in elem.nodes ]

    return K, map, map
end

function elem_mass_matrix(elem::AcusticFluid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)
    J      = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = detJ*ip.w*th
        M    = coef*N*N'
    end

    # map
    map = [ node.dofdict[:up].eq_id for node in elem.nodes  ]

    return M, map, map
end

#TODO 
function elem_RHS_vector(elem::AcusticFluid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    Q      = zeros(nnodes)
    Bw     = zeros(ndim, nnodes)
    KZ     = zeros(ndim)

    J      = Array{Float64}(undef, ndim, ndim)
    dNdX   = Array{Float64}(undef, nnodes, ndim)
    Z      = zeros(ndim) # hydrostatic gradient
    Z[end] = 1.0

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = C'*dNdR
        @gemm dNdX = dNdR*inv(J) # Bw = dNdX'
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Q
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @gemv KZ = K*Z
        @gemm Q += coef*dNdX*KZ
    end

    # map
    map = [ node.dofdict[:up].eq_id for node in elem.nodes  ]

    return Q, map
end

# TODO
function elem_update!(elem::AcusticFluid, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th     = elem.env.thickness

    map_w  = [ node.dofdict[:up].eq_id for node in elem.nodes ]

    C   = getcoords(elem)

    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:up].vals[:up] for node in elem.nodes ]
    Uw += dUw # nodal pore-pressure at step n+1

    dFw = zeros(nnodes)
    Bw  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @gemm dNdX = dNdR*inv(J) # Bw = dNdX'

        G  = dNdX'*Uw/elem.mat.γw # flow gradient
        G[end] += 1.0; # gradient due to gravity

        Δuw = N'*dUw # interpolation to the integ. point

        V = update_state!(elem.mat, ip.state, Δuw, G, Δt)

        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        dFw  -= coef*N*Δuw

        coef = Δt*detJ*ip.w*th
        @gemv dFw += coef*dNdX*V
    end

    DF[map_w] += dFw
    return success()
end

