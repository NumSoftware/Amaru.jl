# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TMSolid

TMSolid_params = [
    FunInfo(:TMSolid, "Thermo-mechanical solid element."),
    KwArgInfo(:rho, "Density", cond=:(rho>=0)),
    KwArgInfo(:gamma, "Specific weight", 0.0, cond=:(gamma>=0)),
    KwArgInfo(:cv, "Heat capacity", cond=:(cv>0)),
    KwArgInfo(:alpha, "Thermal expansion coefficient", cond=:(0<=alpha<=1)),
]
@doc docstring(TMSolid_params) TMSolid(; kwargs...)

struct TMSolidProps<:ElemProperties
    ρ::Float64
    γ::Float64
    cv::Float64
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function TMSolidProps(; kwargs...)
        args = checkargs(kwargs, TMSolid_params)
        this = new(args.rho, args.gamma, args.cv, args.alpha)
        return this
    end    
end


mutable struct TMSolid<:Thermomech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::TMSolidProps
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
    return mech_boundary_forces(elem, facet, key, val)
end


function body_c(elem::TMSolid, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_solid_body_forces(elem, key, val)
end


@inline function elem_map_u(elem::TMSolid)
    keys =(:ux, :uy, :uz)[1:elem.env.ndim]
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


@inline function elem_map_t(elem::TMSolid)
    nbnodes = elem.shape.basic_shape.npoints
    return [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]
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
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

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

    map = elem_map_u(elem)
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
    if elem.env.ana.stressmodel==:planestress
        β = E*α/(1-ν)
        m = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ]
    end

    for ip in elem.ips
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

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
    
    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

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
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

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

    map = elem_map_t(elem)
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
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

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

    map = elem_map_t(elem)

    return M, map, map
end


function elem_internal_forces(elem::TMSolid, F::Array{Float64,1})
    ndim    = elem.env.ndim
    th      = elem.env.ana.thickness
    T0k     = elem.env.ana.T0 + 273.15
    nnodes  = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C       = getcoords(elem)

    E = elem.mat.E
    ν = elem.mat.ν
    ρ = elem.props.ρ
    α = elem.props.α
    cv = elem.props.cv
    β = E*α/(1-2*ν)

    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    m   = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    if elem.env.ana.stressmodel==:planestress
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

    for ip in elem.ips
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in element $(elem.id)")
        invJ = inv(J)
        @mul dNdX = dNdR*invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @mul dNtdX = dNtdR*invJ

        # compute Nt
        Nt = elem.shape.basic_shape.func(ip.R)

        # compute thermal gradient G
        Bt .= dNtdX'

        # internal force dF
        σ  = ip.state.σ
        ut = ip.state.ut

        σ -= β*ut*m # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*σ

        # internal volumes dFt
        ε     = ip.state.ε
        εvol  = dot(m, ε)
        coef  = β*εvol*T0k
        coef *= detJ*ip.w*th
        dFt .-= coef*Nt

        coef  = ρ*cv
        coef *= detJ*ip.w*th
        dFt .-= coef*Nt*ut

        coef  = Δt
        coef *= detJ*ip.w*th
        q = ip.state.q
        @mul dFt += coef*Bt'*q
    end

    F[map_u] = dF
    F[map_t] = dFt

end


function update_elem!(elem::TMSolid, DU::Array{Float64,1}, Δt::Float64)
    ndim    = elem.env.ndim
    th      = elem.env.ana.thickness
    T0k     = elem.env.ana.T0 + 273.15
    nnodes  = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C       = getcoords(elem)

    E = elem.mat.E
    ν = elem.mat.ν
    ρ = elem.props.ρ
    α = elem.props.α
    cv = elem.props.cv
    β = E*α/(1-2*ν)

    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    dU  = DU[map_u] # nodal displacement increments
    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes[1:nbnodes]]
    Ut += dUt # nodal tempeture at step n+1
    m   = I2  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    if elem.env.ana.stressmodel==:planestress
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
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in element $(elem.id)")
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


function post_process(elem::TMSolid)
    # update temperature at non-corner nodes
    elem.shape==elem.shape.basic_shape && return
    npoints  = elem.shape.npoints
    nbpoints = elem.shape.basic_shape.npoints
    Ut = [ elem.nodes[i].vals[:ut] for i in 1:nbpoints ]
    C = elem.shape.nat_coords
    for i in nbpoints+1:npoints
        R = C[i,:]
        N = elem.shape.basic_shape.func(R)
        elem.nodes[i].vals[:ut] = dot(N, Ut)
    end
end