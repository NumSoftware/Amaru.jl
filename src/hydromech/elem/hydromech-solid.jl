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

function elem_config_dofs(elem::HMSolid)
    nbnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)
        if  i<=(nbnodes)
            add_dof(node, :uw, :fw)
        end
    end
end

function elem_init(elem::HMSolid)
    nothing
end


function distributed_bc(elem::HMSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    suitable_keys = (:tx, :ty, :tz, :tn, :tq)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = get_coords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    if key == :tq # fluid volume per area
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
                vip = eval_arith_expr(val, t=t, x=x, y=y)
            else
                x, y, z = X
                vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            end
            coef = vip*nJ*w
            F .+= N*coef # F is a vector
        end

        # generate a map
        map  = [ node.dofdict[:ut].eq_id for node in target.nodes ]

        return F, map
    end

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
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            if key == :tx
                Q = [vip, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, vip, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, vip]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = vip*normalize(n)
            end
        end
        coef = nJ*w
        @gemm F += coef*N*Q' # F is a matrix
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
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C = get_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    Bu = zeros(6, nnodes*ndim)

    DBu = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        setBu(elem.env, dNdX, detJ, Bu)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.mat, ip.data)
        @gemm DBu = D*Bu
        @gemm K += coef*Bu'*DBu
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


# matrix C
function elem_coupling_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = get_coords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cuw = zeros(nnodes*ndim, nbnodes) # u-p coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    for ip in elem.ips

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        setBu(elem.env, dNdX, detJ, Bu)

        # compute Cuw
        Nw    = elem.shape.basic_shape.func(ip.R)
        coef  = elem.mat.α
        coef *= detJ*ip.w*th
        mNw   = m*Nw'
        @gemm Cuw -= coef*Bu'*mNw
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cuw, map_u, map_w
end


function elem_conductivity_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)
    H      = zeros(nbnodes, nbnodes)
    Bw     = zeros(ndim, nbnodes)
    KBw    = zeros(ndim, nbnodes)
    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        dNdR  = elem.shape.deriv(ip.R)
        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        @gemm Bw = inv(J)*dNwdR

        # compute H
        K = calcK(elem.mat, ip.data)
        coef  = 1/elem.mat.γw
        coef *= detJ*ip.w*th
        @gemm KBw = K*Bw
        @gemm H -= coef*Bw'*KBw
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes]  ]

    return H, map, map
end

function elem_compressibility_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)
    Cpp    = zeros(nbnodes, nbnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        Nw   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Cpp
        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        Cpp  -= coef*Nw*Nw'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes]  ]

    return Cpp, map, map
end

function elem_RHS_vector(elem::HMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)
    Q      = zeros(nbnodes)
    Bw     = zeros(ndim, nbnodes)
    KZ     = zeros(ndim)

    J      = Array{Float64}(undef, ndim, ndim)
    Z      = zeros(ndim)
    Z[end] = 1.0 # hydrostatic gradient

    for ip in elem.ips
        dNdR  = elem.shape.deriv(ip.R)
        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        @gemm Bw = inv(J)*dNwdR

        # compute Q
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w*th
        @gemv KZ = K*Z
        @gemm Q += coef*Bw'*KZ
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes]  ]

    return Q, map
end

function elem_internal_forces(elem::HMSolid, F::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = get_coords(elem)
    Cp  = get_coords(elem)[1:nbnodes,:]

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbnodes)
    Bw  = zeros(ndim, nbnodes)

    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)

    for ip in elem.ips

        # compute Bu matrix and Bw
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR
        setBu(elem.env, dNdX, detJ, Bu)

        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm Bw = inv(J)*dNwdR

        # compute N
        Nw   = elem.shape.basic_shape.func(ip.R)

        # internal force
        uw   = ip.data.uw
        σ    = ip.data.σ - elem.mat.α*uw*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ

        # internal volumes dFw
        ε    = ip.data.ε
        εvol = dot(m, ε)
        coef = elem.mat.α*detJ*ip.w*th
        dFw  -= coef*Nw*εvol

        coef = detJ*ip.w*elem.mat.S*th
        dFw -= coef*Nw*uw

        D    = ip.data.D
        coef = detJ*ip.w*th
        @gemv dFw += coef*Bw'*D
    end

    F[map_u] += dF
    F[map_w] += dFw
end


function elem_update!(elem::HMSolid, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes[1:nbnodes] ]
    Uw += dUw # nodal pore-pressure at step n+1
    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbnodes)
    Bw  = zeros(ndim, nbnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, ndim, nnodes)
    dNwdX = Array{Float64}(undef, ndim, nbnodes)
    Δε = zeros(6)

    for ip in elem.ips

        # compute Bu matrix and Bw
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @gemm dNdX = invJ*dNdR
        setBu(elem.env, dNdX, detJ, Bu)

        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm dNwdX = invJ*dNwdR

        # compute Nw
        Nw = elem.shape.basic_shape.func(ip.R)

        # compute Δε
        @gemv Δε = Bu*dU

        # compute Δuw
        Δuw = Nw'*dUw # interpolation to the integ. point

        # Compute flow gradient G
        Bw = dNwdX
        G  = Bw*Uw/elem.mat.γw
        G[end] += 1.0; # gradient due to gravity

        # internal force dF
        Δσ, V = stress_update(elem.mat, ip.data, Δε, Δuw, G, Δt)
        Δσ -= elem.mat.α*Δuw*m # get total stress

        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFw
        Δεvol = dot(m, Δε)
        coef  = elem.mat.α
        coef *= detJ*ip.w*th
        dFw  -= coef*Nw*Δεvol

        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        dFw  -= coef*Nw*Δuw

        coef  = Δt
        coef *= detJ*ip.w*th
        @gemv dFw += coef*Bw'*V
    end

    DF[map_u] += dF
    DF[map_w] += dFw
end
