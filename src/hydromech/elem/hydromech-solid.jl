# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export HMSolid

struct HMSolidProps<:ElemProperties
    ρ::Float64
    γ::Float64
    α::Float64 # Biot's coefficient

    function HMSolidProps(; props...)
        names = (rho="Density", gamma="Specific weight", alpha="Biot coefficient")

        default = (rho=0.0, gamma=0.0, alpha=1.0)
        props   = merge(default, props)
        rho     = props.rho
        gamma   = props.gamma
        alpha   = props.alpha

        @check rho>=0
        @check gamma>=0
        @check 0<alpha<=1

        return new(rho, gamma, alpha)
    end    
end


mutable struct HMSolid<:Hydromech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    props ::HMSolidProps
    active::Bool
    linked_elems::Array{Element,1}
    ctx::Context

    function HMSolid()
        return new()
    end
end

compat_shape_family(::Type{HMSolid}) = BULKCELL
compat_elem_props(::Type{HMSolid}) = HMSolidProps



function elem_config_dofs(elem::HMSolid)
    nbnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.ctx.ndim==3 && add_dof(node, :uz, :fz)
        if  i<=(nbnodes)
            add_dof(node, :uw, :fw)
        end
    end
end


function elem_init(::HMSolid)
end


function distributed_bc(elem::HMSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.ctx.ndim
    th    = elem.ctx.thickness
    suitable_keys = (:tx, :ty, :tz, :tn, :tq)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet !== nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.ctx.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    shape = target.shape
    ips   = get_ip_coords(shape)

    if key == :tq # fluid volume per area
        F = zeros(nnodes)
        for i in 1:size(ips,1)
            R = ips[i].coord
            w = ips[i].w
            N = shape.func(R)
            D = shape.deriv(R)
            J = C'*D
            nJ = norm2(J)
            X = C'*N
            if ndim==2
                x, y = X
                vip = evaluate(val, t=t, x=x, y=y)
            else
                x, y, z = X
                vip = evaluate(val, t=t, x=x, y=y, z=z)
            end
            coef = vip*nJ*w
            F .+= N*coef # F is a vector
        end

        # generate a map
        map  = [ node.dofdict[:uw].eq_id for node in target.nodes ]

        return F, map
    end

    F = zeros(nnodes, ndim)
    for i in 1:size(ips,1)

        R = ips[i].coord
        w = ips[i].w
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
        X = C'*N
        if ndim==2
            x, y = X
            vip = evaluate(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
            if elem.ctx.stressmodel==:axisymmetric
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = evaluate(val, t=t, x=x, y=y, z=z)
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
        coef = norm2(J)*w*th
        @mul F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


@inline function setBu(elem::HMSolid, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end


function elem_stiffness(elem::HMSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    Bu = zeros(6, nnodes*ndim)

    DBu = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        setBu(elem, ip, dNdX, Bu)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.mat, ip.state)
        @mul DBu = D*Bu
        @mul K += coef*Bu'*DBu
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


# matrix C
function elem_coupling_matrix(elem::HMSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cuw = zeros(nnodes*ndim, nbnodes) # u-p coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    m = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        setBu(elem, ip, dNdX, Bu)

        # compute Cuw
        Nw    = elem.shape.basic_shape.func(ip.R)
        coef  = elem.props.α
        coef *= detJ*ip.w*th
        mNw   = m*Nw'
        @mul Cuw -= coef*Bu'*mNw
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cuw, map_u, map_w
end


function elem_conductivity_matrix(elem::HMSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    H      = zeros(nbnodes, nbnodes)
    Bw     = zeros(ndim, nbnodes)
    KBw    = zeros(ndim, nbnodes)
    J      = Array{Float64}(undef, ndim, ndim)
    dNwdX  = Array{Float64}(undef, nbnodes, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @mul J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNwdX = dNwdR*inv(J)
	    Bw .= dNwdX'

        # compute H
        K = calcK(elem.mat, ip.state)
        coef  = 1/elem.ctx.γw
        coef *= detJ*ip.w*th
        @mul KBw = K*Bw
        @mul H -= coef*Bw'*KBw
    end

    # map
    map = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes]  ]

    return H, map, map
end


function elem_compressibility_matrix(elem::HMSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    Cpp    = zeros(nbnodes, nbnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        Nw   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

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
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    Q      = zeros(nbnodes)
    Bw     = zeros(ndim, nbnodes)
    KZ     = zeros(ndim)

    J      = Array{Float64}(undef, ndim, ndim)
    dNwdX  = Array{Float64}(undef, nbnodes, ndim)
    Z      = zeros(ndim)
    Z[end] = 1.0 # hydrostatic gradient

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @mul J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNwdX = dNwdR*inv(J)
        Bw .= dNwdX'

        # compute Q
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @mul KZ = K*Z
        @mul Q += coef*Bw'*KZ
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes]  ]

    return Q, map
end


function elem_internal_forces(elem::HMSolid, F::Array{Float64,1})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    Cp  = getcoords(elem)[1:nbnodes,:]

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbnodes)
    Bw  = zeros(ndim, nbnodes)

    m = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, nnodes, ndim)
    dNwdX = Array{Float64}(undef, nbnodes, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu matrix and Bw
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR*inv(J)
        setBu(elem, ip, dNdX, Bu)

        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @mul dNwdX = dNwdR*inv(J)
        Bw .= dNwdX'

        # compute N
        Nw   = elem.shape.basic_shape.func(ip.R)

        # internal force
        uw   = ip.state.uw
        σ    = ip.state.σ - elem.props.α*uw*m # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*σ

        # internal volumes dFw
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = elem.props.α*detJ*ip.w*th
        dFw  .-= coef*Nw*εvol

        coef = detJ*ip.w*elem.mat.S*th
        dFw .-= coef*Nw*uw

        D    = ip.state.D
        coef = detJ*ip.w*th
        @mul dFw += coef*Bw'*D
    end

    F[map_u] += dF
    F[map_w] += dFw
end


function update_elem!(elem::HMSolid, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUw = DU[map_w] # nodal pore-pressure increments
    Uw  = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes[1:nbnodes] ]
    Uw += dUw # nodal pore-pressure at step n+1
    m = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFw = zeros(nbnodes)
    Bw  = zeros(ndim, nbnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, nnodes, ndim)
    dNwdX = Array{Float64}(undef, nbnodes, ndim)
    Δε = zeros(6)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu matrix and Bw
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @mul dNdX = dNdR*invJ
        setBu(elem, ip, dNdX, Bu)

        dNwdR = elem.shape.basic_shape.deriv(ip.R)
        @mul dNwdX = dNwdR*invJ
        Bw .= dNwdX'

        # compute Nw
        Nw = elem.shape.basic_shape.func(ip.R)

        # compute Δε
        @mul Δε = Bu*dU

        # compute Δuw
        Δuw = Nw'*dUw # interpolation to the integ. point

        # Compute flow gradient G
        # Bw = dNwdX
        G  = Bw*Uw/elem.ctx.γw
        G[end] += 1.0; # gradient due to gravity

        # internal force dF
        Δσ, V = update_state!(elem.mat, ip.state, Δε, Δuw, G, Δt)
        Δσ -= elem.props.α*Δuw*m # get total stress

        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*Δσ

        # internal volumes dFw
        Δεvol = dot(m, Δε)
        coef  = elem.props.α
        coef *= detJ*ip.w*th
        dFw  .-= coef*Nw*Δεvol

        coef  = elem.mat.S
        coef *= detJ*ip.w*th
        dFw  .-= coef*Nw*Δuw

        coef  = Δt
        coef *= detJ*ip.w*th
        @mul dFw += coef*Bw'*V
    end

    return [dF; dFw], [map_u; map_w], success()
end
