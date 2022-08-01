# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct SeepSolid<:Hydromechanical
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function SeepSolid();
        return new()
    end
end

matching_shape_family(::Type{SeepSolid}) = SOLID_CELL

function elem_config_dofs(elem::SeepSolid)
    for node in elem.nodes
        add_dof(node, :uw, :fw)
    end
end

function elem_init(elem::SeepSolid)
    nothing
end


function distributed_bc(elem::SeepSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.thickness
    suitable_keys = (:tq,) # tq: fluid volumes per area

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a SeepSolid element")

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
        nJ = norm2(J)
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
        coef = vip*norm(J)*w*th
        F .+= coef*N # F is a vector
    end

    # generate a map
    map  = [ node.dofdict[:uw].eq_id for node in target.nodes ]

    return F, map
end


# hydraulic conductivity
function elem_conductivity_matrix(elem::SeepSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    dNdX   = zeros(nnodes, ndim)
    KBw    = zeros(ndim, nnodes)

    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        @gemm dNdX = dNdR*inv(J) # Bw = dNdX'

        # compute H
        K = calcK(elem.mat, ip.state)
        coef  = 1/elem.mat.γw
        coef *= detJ*ip.w*th
        @gemm KBw = K*dNdX'
        @gemm H -= coef*dNdX*KBw
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map
end

function elem_compressibility_matrix(elem::SeepSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    Cpp    = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Cpp
        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        Cpp  -= coef*N*N'
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Cpp, map, map
end

function elem_RHS_vector(elem::SeepSolid)
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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Q
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @gemv KZ = K*Z
        @gemm Q += coef*dNdX*KZ
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end

function elem_internal_forces(elem::SeepSolid, F::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C   = getcoords(elem)

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dFw = zeros(nnodes)
    Bw  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bw matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = dNdR*inv(J)

        # Bw = copy(dNdX')

        # compute N
        N    = elem.shape.func(ip.R)

        # internal volumes dFw
        uw   = ip.state.uw
        coef = detJ*ip.w*elem.mat.S*th
        dFw -= coef*N*uw

        D    = ip.state.D
        coef = detJ*ip.w*th
        @gemv dFw += coef*dNdX*D
    end

    F[map_w] += dFw
end

function elem_update!(elem::SeepSolid, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th     = elem.env.thickness

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    C   = getcoords(elem)

    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ]
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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
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

