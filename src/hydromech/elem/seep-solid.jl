# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


mutable struct SeepSolid<:Hydromechanical
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

    function SeepSolid();
        return new()
    end
end

matching_shape_family(::Type{SeepSolid}) = SOLID_SHAPE

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
    C = nodes_coords(nodes, ndim)

    # Calculate the nodal values
    F     = zeros(nnodes)
    J     = Array{Float64}(undef, ndim, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        @gemm J = D*C
        nJ = norm2(J)
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        end
        coef = vip*nJ*w
        F .+= N*coef # F is a vector
    end

    # generate a map
    map  = [ node.dofdict[:uw].eq_id for node in target.nodes ]

    return F, map
end


# conductivity
function elem_conductivity_matrix(elem::SeepSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = elem_coords(elem)
    H      = zeros(nnodes, nnodes)
    Bw     = zeros(ndim, nnodes)
    KBw    = zeros(ndim, nnodes)

    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        @gemm Bw = inv(J)*dNdR

        # compute H
        K = calcK(elem.mat, ip.data)
        coef  = 1/elem.mat.γw
        coef *= detJ*ip.w*th
        @gemm KBw = K*Bw
        @gemm H -= coef*Bw'*KBw
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map
end

function elem_compressibility_matrix(elem::SeepSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C      = elem_coords(elem)
    Cpp    = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
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
    nnodes = length(elem.nodes)
    C      = elem_coords(elem)
    Q      = zeros(nnodes)
    Bw     = zeros(ndim, nnodes)
    KZ     = zeros(ndim)

    J      = Array{Float64}(undef, ndim, ndim)
    dNdX   = Array{Float64}(undef, ndim, nnodes)
    Z      = zeros(ndim) # hydrostatic gradient
    Z[end] = 1.0

    for ip in elem.ips

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = dNdR*C
        @gemm Bw = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Q
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w
        @gemv KZ = K*Z
        @gemm Q += coef*Bw'*KZ
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end

function elem_internal_forces(elem::SeepSolid, F::Array{Float64,1})
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C   = elem_coords(elem)

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dFw = zeros(nnodes)
    Bw  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute Bw matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR

        Bw = dNdX

        # compute N
        N    = elem.shape.func(ip.R)

        # internal volumes dFw
        uw   = ip.data.uw
        coef = detJ*ip.w*elem.mat.S
        dFw -= coef*N*uw

        D    = ip.data.D
        coef = detJ*ip.w
        @gemv dFw += coef*Bw'*D
    end

    F[map_w] += dFw
end

function elem_update!(elem::SeepSolid, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th     = elem.env.thickness

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    C   = elem_coords(elem)

    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ]
    Uw += dUw # nodal pore-pressure at step n+1

    dF  = zeros(nnodes*ndim)
    dFw = zeros(nnodes)
    Bw  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute Bu matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR

        Bw = dNdX
        G  = Bw*Uw/elem.mat.γw # flow gradient
        G[end] += 1.0; # gradient due to gravity

        Δuw = N'*dUw # interpolation to the integ. point

        V = update_state!(elem.mat, ip.data, Δuw, G, Δt)

        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        dFw  -= coef*N*Δuw

        coef = Δt*detJ*ip.w*th
        @gemv dFw += coef*Bw'*V
    end

    DF[map_w] += dFw
end

