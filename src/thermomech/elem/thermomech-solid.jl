# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TMSolid

struct TMSolidProps<:ElemProperties
    ρ::Float64
    γ::Float64
    cv::Float64

    function TMSolidProps(;rho=0.0, gamma=0.0, cv=0.0)
        @check rho>=0
        @check gamma>=0
        @check cv>=0

        return new(rho, gamma, cv)
    end    
end

TMSolid = TMSolidProps


mutable struct TMSolidElem<:ThermomechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    props::TMSolidProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function TMSolidElem(props=TMSolidProps())
        this = new()
        this.props = props
        return this
    end
end

matching_shape_family(::Type{TMSolidElem}) = BULKCELL
matching_elem_type(::Type{TMSolidProps}) = TMSolidElem
matching_props_type(::Type{TMSolidElem}) = TMSolidProps


function elem_config_dofs(elem::TMSolidElem)
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


function elem_init(elem::TMSolidElem)
end


function distributed_bc(elem::TMSolidElem, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.anaprops.thickness
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
            R = vec(ips[i,:])
            w = R[end]
            N = shape.func(R)
            D = shape.deriv(R)
            J = C'*D
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

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
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
            if elem.env.anaprops.stressmodel=="axisymmetric"
                th = 2*pi*X[1]
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
        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


@inline function set_Bu(elem::TMSolidElem, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end


function elem_stiffness(elem::TMSolidElem)
    ndim   = elem.env.ndim
    th     = elem.env.anaprops.thickness
    nnodes = length(elem.nodes)
    C  = getcoords(elem)
    K  = zeros(nnodes*ndim, nnodes*ndim)
    Bu = zeros(6, nnodes*ndim)

    DBu = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.matparams, ip.state)
        @gemm DBu = D*Bu
        @gemm K += coef*Bu'*DBu
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


# matrix C
function elem_coupling_matrix(elem::TMSolidElem)
    ndim   = elem.env.ndim
    th     = elem.env.anaprops.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = elem.matparams.E*elem.matparams.α/(1-2*elem.matparams.nu) # thermal stress modulus

    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        # compute Cut
        Nt    = elem.shape.basic_shape.func(ip.R)
        coef  = β
        coef *= detJ*ip.w*th
        mNt   = m*Nt'
        @gemm Cut -= coef*Bu'*mNt
    end
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_t = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cut, map_u, map_t
end

# thermal conductivity
function elem_conductivity_matrix(elem::TMSolidElem)
    ndim   = elem.env.ndim
    th     = elem.env.anaprops.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    dNtdX  = zeros(nnodes, ndim)
    Bt     = zeros(ndim, nnodes)
    KBt    = zeros(ndim, nnodes)
    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @gemm dNtdX = dNtdR*inv(J)
        Bt .= dNtdX'

        # compute H
        K = calcK(elem.matparams, ip.state)
        coef = detJ*ip.w*th
        @gemm KBt = K*Bt
        @gemm H  -= coef*Bt'*KBt
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return H, map, map
end

function elem_mass_matrix(elem::TMSolidElem)
    ndim   = elem.env.ndim
    th     = elem.env.anaprops.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        Nt   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
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

#=
function elem_internal_forces(elem::TMSolidElem, F::Array{Float64,1}, DU::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = elem.env.anaprops.thickness # VERIFICAR ESPESSURA
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    T0k     = elem.env.T0k + 273.15

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nbnodes)
    Bt  = zeros(ndim, nbnodes)

    m = [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ] # = tI

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Jp  = Array{Float64}(undef, ndim, nbnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    dUt = DU[map_t] # nodal temperature increments
    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix and Bt
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @gemm dNdX = dNdR*inv(J)
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        Jp = dNtdR*Ct
        @gemm dNtdX = inv(Jp)*dNtdR
        Bt = dNtdX
        # compute N

        # internal force
        ut   = ip.state.ut + 273
        β   = elem.matparams.E*elem.matparams.α/(1-2*elem.matparams.nu)
        σ    = ip.state.σ - β*ut*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ

        # internal volumes dFt
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = β*detJ*ip.w*th
        dFt  -= coef*Nt*εvol

        coef = detJ*ip.w*elem.props.ρ*elem.props.cv*th/T0k
        dFt -= coef*Nt*ut

        QQ   = ip.state.QQ
        coef = detJ*ip.w*th/T0k
        @gemv dFt += coef*Bt'*QQ
    end

    F[map_u] += dF
    F[map_t] += dFt
end
=#


function update_elem!(elem::TMSolidElem, DU::Array{Float64,1}, Δt::Float64)
    ndim    = elem.env.ndim
    th      = elem.env.anaprops.thickness
    T0k     = elem.env.anaprops.T0 + 273.15
    nnodes  = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C       = getcoords(elem)

    E = elem.matparams.E
    α = elem.matparams.α
    ρ = elem.props.ρ
    nu = elem.matparams.nu
    cv = elem.props.cv
    β = E*α/(1-2*nu)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1
    m   = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, nnodes, ndim)
    dNtdX = Array{Float64}(undef, nbnodes, ndim)
    Δε = zeros(6)

    for ip in elem.ips
        elem.env.anaprops.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @gemm dNdX = dNdR*invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm dNtdX = dNtdR*invJ

        # compute Nt
        Nt = elem.shape.basic_shape.func(ip.R)

        # compute Δε
        @gemv Δε = Bu*dU

        # compute Δut
        Δut = Nt'*dUt # interpolation to the integ. point

        # compute thermal gradient G
        Bt .= dNtdX'
        G  = Bt*Ut

        # internal force dF
        Δσ, q = update_state(elem.matparams, ip.state, Δε, Δut, G, Δt)
        #@show Δσ

        #@show "HIIIIIIIIIIIIIIIIIII"
        Δσ -= β*Δut*m # get total stress
        #@show Δσ

        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0k
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt

        coef  = ρ*cv
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt*Δut

        coef  = Δt
        coef *= detJ*ip.w*th
        @gemv dFt += coef*Bt'*q
    end

    return [dF; dFt], [map_u; map_t], success()
end
