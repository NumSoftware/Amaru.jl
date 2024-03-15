# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TMSolid

struct TMSolidProps<:ElemProperties
    ρ::Float64
    γ::Float64
    cv::Float64
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function TMSolidProps(; props...)
        names = (rho="Density", gamma="Specific weight", cv="Heat capacity", alpha="Thermal expansion coefficient")
        required = (:cv, :alpha)
        @checkmissing props required names

        default = (rho=0.0, gamma=0.0)
        props  = merge(default, props)

        cv  = props.cv
        rho = props.rho
        gamma = props.gamma
        alpha = props.alpha

        @check cv>0.0
        @check rho>=0.0
        @check gamma>=0.0
        @check 0<=alpha<=1

        return new(rho, gamma, cv, alpha)
    end    
end


mutable struct TMSolid<:Thermomech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    props::TMSolidProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function TMSolid()
        return new()
    end
end

compat_shape_family(::Type{TMSolid}) = BULKCELL
compat_elem_props(::Type{TMSolid}) = TMSolidProps


function elem_init(elem::TMSolid)
end

function elem_config_dofs(elem::TMSolid)
    nbnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)
        if  i<=(nbnodes)
            add_dof(node, :ut, :ft)
        end
    end
end

function distributed_bc(elem::TMSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.ana.thickness
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
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    if key == :tq # energy per area
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
        map  = [ node.dofdict[:ut].eq_id for node in target.nodes ]

        return F, map
    end

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
            if elem.env.ana.stressmodel=="axisymmetric"
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


@inline function set_Bu(elem::TMSolid, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end


function elem_stiffness(elem::TMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    C  = getcoords(elem)
    K  = zeros(nnodes*ndim, nnodes*ndim)
    Bu = zeros(6, nnodes*ndim)

    DBu = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

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
function elem_coupling_matrix(elem::TMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    E = elem.mat.E
    ν = elem.mat.ν
    α = elem.props.α

    C   = getcoords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    m    = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = E*α/(1-2*ν) # thermal stress modulus
    if elem.env.ana.stressmodel=="plane-stress"
        β = E*α/(1-ν)
        m = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ]
    end

    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        # compute Cut
        Nt    = elem.shape.basic_shape.func(ip.R)
        coef  = β
        coef *= detJ*ip.w*th
        mNt   = m*Nt'
        @mul Cut -= coef*Bu'*mNt
    end
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_t = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cut, map_u, map_t
end

# thermal conductivity
function elem_conductivity_matrix(elem::TMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    H      = zeros(nbnodes, nbnodes)
    dNtdX  = zeros(nbnodes, ndim)
    Bt     = zeros(ndim, nbnodes)
    KBt    = zeros(ndim, nbnodes)
    J      = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @mul J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNtdX = dNtdR*inv(J)
        Bt .= dNtdX'

        # compute H
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @mul KBt = K*Bt
        @mul H  -= coef*Bt'*KBt
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return H, map, map
end


function elem_mass_matrix(elem::TMSolid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    M      = zeros(nbnodes, nbnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        Nt   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Cut
        coef  = elem.props.ρ*elem.props.cv
        coef *= detJ*ip.w*th
        M    -= coef*Nt*Nt'
    end

    # map
    map = [  node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes]  ]

    return M, map, map
end


function update_elem!(elem::TMSolid, DU::Array{Float64,1}, Δt::Float64)
    ndim    = elem.env.ndim
    th      = elem.env.ana.thickness
    T0k     = elem.env.ana.T0 + 273.15
    nnodes  = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C       = getcoords(elem)

    E = elem.mat.E
    α = elem.props.α
    ρ = elem.props.ρ
    ν = elem.mat.ν
    cv = elem.props.cv
    β = E*α/(1-2*ν)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes[1:nbnodes]]
    Ut += dUt # nodal tempeture at step n+1
    m   = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    if elem.env.ana.stressmodel=="plane-stress"
        β = E*α/(1-ν)
        m = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ]
    end

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nbnodes)
    Bt  = zeros(ndim, nbnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, nnodes, ndim)
    dNtdX = Array{Float64}(undef, nbnodes, ndim)
    Δε = zeros(6)

    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @mul dNdX = dNdR*invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @mul dNtdX = dNtdR*invJ

        # compute Nt
        Nt = elem.shape.basic_shape.func(ip.R)

        # compute Δε
        @mul Δε = Bu*dU

        # compute Δut
        Δut = Nt'*dUt # interpolation to the integ. point

        # compute thermal gradient G
        Bt .= dNtdX'
        G  = Bt*Ut

        # internal force dF
        Δσ, q, status = update_state!(elem.mat, ip.state, Δε, Δut, G, Δt)
        failed(status) && return [dF; dFt], [map_u; map_t], status


        Δσ -= β*Δut*m # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0k
        coef *= detJ*ip.w*th
        dFt  .-= coef*Nt

        coef  = ρ*cv
        coef *= detJ*ip.w*th
        dFt  .-= coef*Nt*Δut

        coef  = Δt
        coef *= detJ*ip.w*th
        @mul dFt += coef*Bt'*q
    end

    return [dF; dFt], [map_u; map_t], success()
end