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
    nbsnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)     
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)    
        if  i<=(nbsnodes)
            add_dof(node, :uw, :fw)    
        end
    end
end

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
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C = elem_coords(elem)
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
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    C   = elem_coords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cup = zeros(nnodes*ndim, nbsnodes) # u-p coupling matrix

    J  = Array{Float64}(undef, ndim, ndim)
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

        # compute Cup
        Np   = elem.shape.basic_shape.func(ip.R)
        coef = detJ*ip.w*elem.mat.α 
        mNp  = m*Np'
        @gemm Cup -= coef*Bu'*mNp
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes] ]

    return Cup, map_u, map_p
end


function elem_conductivity_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    Cp     = elem_coords(elem)[1:nbsnodes,:]
    H      = zeros(nbsnodes, nbsnodes)
    Bp     = zeros(ndim, nbsnodes)
    KBp    = zeros(ndim, nbsnodes)

    J    = Array{Float64}(undef, ndim, ndim)
    dNpdX = Array{Float64}(undef, ndim, nbsnodes)
    nodes_p = elem.nodes[1:nbsnodes]

    for ip in elem.ips

        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = dNpdR*Cp
        @gemm Bp = inv(J)*dNpdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute H
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w/elem.mat.γw
        @gemm KBp = K*Bp
        @gemm H -= coef*Bp'*KBp
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes]  ]

    return H, map, map, nodes_p
end

function elem_compressibility_matrix(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    Cp  = elem_coords(elem)[1:nbsnodes,:]
    Cpp = zeros(nbsnodes, nbsnodes) 

    J  = Array{Float64}(undef, ndim, ndim)


    for ip in elem.ips

        Np    = elem.shape.basic_shape.func(ip.R)
        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J = dNpdR*Cp
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Cuu
        coef = detJ*ip.w*elem.mat.S 
        Cpp  -= coef*Np*Np'
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes]  ]

    return Cpp, map, map
end

function elem_RHS_vector(elem::HMSolid)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    Cp     = elem_coords(elem)[1:nbsnodes,:]
    Q      = zeros(nbsnodes)
    Bp     = zeros(ndim, nbsnodes)
    KZ     = zeros(ndim)

    J    = Array{Float64}(undef, ndim, ndim)
    dNpdX = Array{Float64}(undef, ndim, nbsnodes)
    Z      = zeros(ndim) 
    Z[end] = 1.0 # hydrostatic gradient

    for ip in elem.ips

        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = dNpdR*Cp
        @gemm Bp = inv(J)*dNpdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Q
        K = calcK(elem.mat, ip.data)
        coef = detJ*ip.w
        @gemv KZ = K*Z
        @gemm Q += coef*Bp'*KZ
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes]  ]

    return Q, map
end

function elem_internal_forces(elem::HMSolid, F::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    C   = elem_coords(elem)
    Cp  = elem_coords(elem)[1:nbsnodes,:]

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes] ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbsnodes)
    Bp  = zeros(ndim, nbsnodes)

    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    Jp  = Array{Float64}(undef, ndim, nbsnodes)
    dNpdX = Array{Float64}(undef, ndim, nbsnodes)

    for ip in elem.ips

        # compute Bu matrix and Bp
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR
        setBu(elem.env, dNdX, detJ, Bu)

        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        Jp = dNpdR*Cp
        @gemm dNpdX = inv(Jp)*dNpdR
        Bp = dNpdX

        # compute N
        Np   = elem.shape.basic_shape.func(ip.R)

        # internal force 
        uw   = ip.data.uw
        σ    = ip.data.σ - elem.mat.α*uw*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ

        # internal volumes dFw
        ε    = ip.data.ε
        εvol = dot(m, ε)
        coef = elem.mat.α*detJ*ip.w
        dFw  -= coef*Np*εvol
        
        coef = detJ*ip.w*elem.mat.S 
        dFw -= coef*Np*uw  

        D    = ip.data.D
        coef = detJ*ip.w
        @gemv dFw += coef*Bp'*D
    end

    F[map_u] += dF
    F[map_p] += dFw
end

function elem_update!(elem::HMSolid, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    C   = elem_coords(elem)
    Cp   = elem_coords(elem)[1:nbsnodes,:]

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbsnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUw = DU[map_p] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes[1:nbsnodes] ] 
    Uw += dUw # nodal pore-pressure at step n+1
    m = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbsnodes)
    Bp  = zeros(ndim, nbsnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, ndim, nnodes)
    Jp = Array{Float64}(undef, ndim, ndim)
    dNpdX = Array{Float64}(undef, ndim, nbsnodes)
    Δε = zeros(6)

    for ip in elem.ips

        # compute Bu matrix and Bp
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR
        setBu(elem.env, dNdX, detJ, Bu)

        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm Jp = dNpdR*Cp
        @gemm dNpdX = inv(Jp)*dNpdR
        Bp = dNpdX

        # compute Np
        Np   = elem.shape.basic_shape.func(ip.R)

       	# Compute Δε 
        @gemv Δε = Bu*dU

        # Compute Δuw
        Δuw = Np'*dUw # interpolation to the integ. point

        # Compute flow gradient G 
        G  = Bp*Uw/elem.mat.γw
        G[end] += 1.0; # gradient due to gravity

        # internal force dF
        Δσ, V = stress_update(elem.mat, ip.data, Δε, Δuw, G, Δt)
        Δσ -= elem.mat.α*Δuw*m # get total stress

        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFw
        Δεvol = dot(m, Δε)
        coef  = elem.mat.α*detJ*ip.w
        dFw  -= coef*Np*Δεvol

        coef = detJ*ip.w*elem.mat.S 
        dFw -= coef*Np*Δuw     

        coef = Δt*detJ*ip.w
        @gemv dFw += coef*Bp'*V
    end

    DF[map_u] += dF
    DF[map_p] += dFw
end

