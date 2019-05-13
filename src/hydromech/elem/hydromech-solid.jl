# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


mutable struct HMSolid<:Hydromechanical
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

    function HMSolid(); 
        return new() 
    end
end

matching_shape_family(::Type{HMSolid}) = SOLID_SHAPE

function elem_init(elem::HMSolid)
    nothing
end


function distributed_bc(elem::HMSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim

    # Check bcs
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")
    !(key in (:tx, :ty, :tz, :tn)) && error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    # TODO: add tq boundary condition (fluid volume per area)

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = nodes_coords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
        nJ = norm2(J)
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*n/norm(n)
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x)
            if key == :tx
                Q = [vip, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, vip, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, vip]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = vip*n/norm(n)
            end
        end
        F += N*Q'*(nJ*w) # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


function setBu(env::ModelEnv, dNdX::Matx, detJ::Float64, B::Matx)
    ndim, nnodes = size(dNdX)
    B .= 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[6,1+j*ndim] = dNdX[2,i]/SR2; B[6,2+j*ndim] = dNdX[1,i]/SR2
        end
        if env.modeltype==:axisymmetric
            for i in 1:nnodes
                N =elem.shape.func(R)
                j = i-1
                r = R[0]
                B[1,1+j*ndim] = dNdX[1,i]
                B[2,2+j*ndim] = dNdX[2,i]
                B[3,1+j*ndim] =    N[i]/r
                B[6,1+j*ndim] = dNdX[2,i]/SR2; B[6,2+j*ndim] = dNdX[1,i]/SR2
            end
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,2+j*ndim] = dNdz/SR2;   B[4,3+j*ndim] = dNdy/SR2
            B[5,1+j*ndim] = dNdz/SR2;   B[5,3+j*ndim] = dNdx/SR2
            B[6,1+j*ndim] = dNdy/SR2;   B[6,2+j*ndim] = dNdx/SR2
        end
    end

    return detJ
end


function elem_stiffness(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C = elem_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        setBu(elem.env, dNdX, detJ, B)

        # compute K
        coef = detJ*ip.w
        D    = calcD(elem.mat, ip.data) 
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


# matrix C
function elem_coupling_matrix(elem::HMSolid) 
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C   = elem_coords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cup = zeros(nnodes*ndim, nnodes) # u-p coupling matrix

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    for ip in elem.ips

        # compute Bu matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        setBu(elem.env, dNdX, detJ, Bu)

        # compute K
        coef = detJ*ip.w*elem.mat.α 
        mNt  = m*N'

        @gemm Cup -= coef*Bu'*mNt
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    return Cup, map_u, map_p
end


function elem_conductivity_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C      = elem_coords(elem)
    H      = zeros(nnodes, nnodes)
    Bp     = zeros(ndim, nnodes)
    KBp    = zeros(ndim, nnodes)

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = dNdR*C
        @gemm Bp = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute H
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w/elem.mat.γw
        @gemm KBp = K*Bp
        @gemm H -= coef*Bp'*KBp
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map
end

function elem_compressibility_matrix(elem::HMSolid)
    if elem.mat.S == 0.0
        return zeros(0,0), zeros(Int64,0), zeros(Int64,0)
    end

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C   = elem_coords(elem)
    Cpp = zeros(nnodes, nnodes) 

    J  = Array{Float64}(undef, ndim, ndim)


    for ip in elem.ips

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Cuu
        coef = detJ*ip.w*elem.mat.S 
        Cpp  -= coef*N*N'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Cpp, map, map
end

function elem_RHS_vector(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    C      = elem_coords(elem)
    Q      = zeros(nnodes)
    Bp     = zeros(ndim, nnodes)
    KZ     = zeros(ndim)

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    Z      = zeros(ndim) 
    Z[end] = 1.0 # hydrostatic gradient

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        @gemm J  = dNdR*C
        @gemm Bp = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Q
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w
        @gemv KZ = K*Z
        @gemm Q += coef*Bp'*KZ
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end


function elem_update!(elem::HMSolid, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    C   = elem_coords(elem)

    dU  = DU[map_u] # nodal displacement increments
    dUw = DU[map_p] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ] 
    Uw += dUw # nodal pore-pressure at step n+1
    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nnodes)
    Bp  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    Δε = zeros(6)

    for ip in elem.ips

        # compute Bu matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR
        setBu(elem.env, dNdX, detJ, Bu)

       	# Compute Δε 
        @gemv Δε = Bu*dU

        # Compute Δuw
        Δuw = N'*dUw # interpolation to the integ. point

        # Compute flow gradient G 
        Bp = dNdX
        G  = Bp*Uw/elem.mat.γw
        G[end] += 1.0; # gradient due to gravity

        # internal force dF
        Δσ, V = stress_update(elem.mat, ip.data, Δε, Δuw, G)
        Δσ -= elem.mat.α*Δuw*m # get total stress

        coef = detJ*ip.w
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFw
        Δεvol = dot(m, Δε)
        coef  = elem.mat.α*Δεvol*detJ*ip.w
        dFw  -= coef*N

        if elem.mat.S != 0.0
            coef = elem.mat.S*Δuw*detJ*ip.w 
            dFw -= coef*N       
        end

        coef = Δt*detJ*ip.w
        @gemv dFw += coef*Bp'*V
    end

    DF[map_u] += dF
    DF[map_p] += dFw
end

