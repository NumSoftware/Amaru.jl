# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


export ThermoSolid

ThermoSolid_params = [
    FunInfo(:ThermoSolid, "Finite element for thermal analyses"),
    KwArgInfo(:rho, "Density", 0.0, cond=:(rho>=0)),
    KwArgInfo(:cv, "Heat capacity", cond=:(cv>0)),
]
@doc docstring(ThermoSolid_params) ThermoSolid(; kwargs...)

struct ThermoSolidProps<:ElemProperties
    ρ::Float64  # material specific weight Ton/m3
    cv::Float64 # Specific heat J/Ton/k

    function ThermoSolidProps(; kwargs...)
        args = checkargs(kwargs, ThermoSolid_params)
        this = new(args.rho, args.cv)

        return this
    end    
end


mutable struct ThermoSolid<:ThermoMech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    props ::ThermoSolidProps
    active::Bool
    linked_elems::Array{Element,1}
    ctx::Context

    function ThermoSolid()
        return new()
    end
end


compat_shape_family(::Type{ThermoSolid}) = BULKCELL
compat_elem_props(::Type{ThermoSolid}) = ThermoSolidProps


function elem_config_dofs(elem::ThermoSolid)
    for node in elem.nodes
        add_dof(node, :ut, :ft)
    end
end


function elem_init(elem::ThermoSolid)
    nothing
end


function distributed_bc(elem::ThermoSolid, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.ctx.ndim
    th    = elem.ctx.thickness
    suitable_keys = (:tq,)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a ThermoSolid element")

    target = facet !== nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.ctx.t

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
        R = ips[i].coord
        w = ips[i].w
        N = shape.func(R)
        D = shape.deriv(R)
        @mul J = C'*D
        nJ = norm2(J)
        X = C'*N
        if ndim==2
            x, y = X
            vip = evaluate(val, t=t, x=x, y=y)
            if elem.ctx.stressmodel==:axisymmetric
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = evaluate(val, t=t, x=x, y=y, z=z)
        end
        coef = vip*norm(J)*w*th
        F .+= coef*N # F is a vector
    end

    # generate a map
    map  = [ node.dofdict[:ut].eq_id for node in target.nodes ]

    return F, map
end


# thermal conductivity
function elem_conductivity_matrix(elem::ThermoSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    dNdX   = zeros(nnodes, ndim)
    Bt     = zeros(ndim, nnodes)
    KBt    = zeros(ndim, nnodes)

    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)
        dNdR = elem.shape.deriv(ip.R)
        @mul J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        dNdX = dNdR*inv(J)
        Bt .= dNdX'

        # compute H
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @mul KBt = K*Bt
        @mul H -= coef*Bt'*KBt
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes ]

    return H, map, map
end


function elem_mass_matrix(elem::ThermoSolid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Cut
        coef  = elem.props.ρ*elem.props.cv
        coef *= detJ*ip.w*th
        M    -= coef*N*N'
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes ]

    return M, map, map
end

#=
function elem_internal_forces(elem::ThermoSolid, F::Array{Float64,1})
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    C   = getcoords(elem)
    th     = elem.ctx.thickness
    θ0     = elem.ctx.T0 + 273.15
    map_p  = [ node.dofdict[:ut].eq_id for node in elem.nodes ]

    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)


        # compute Bt matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @mul dNdX = dNdR*inv(J)

        Bt = dNdX

        # compute N
        N    = elem.shape.func(ip.R)

        # internal volumes dFw
        ut   = ip.state.ut
        coef = detJ*ip.w*elem.props.cv*elem.props.ρ*th # VERIFICAR
        dFt -= coef*N*ut

        D    = ip.state.D
        coef = detJ*ip.w*th/θ0
        @mul dFt -= coef*Bt'*D
    end

    F[map_p] += dFt
end
=#


function update_elem!(elem::ThermoSolid, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    th     = elem.ctx.thickness

    map_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes ]

    C   = getcoords(elem)

    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes ]
    Ut += dUt # nodal temperature at step n+1

    Bt  = zeros(ndim, nnodes)
    dFt = zeros(nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @mul dNdX = dNdR*inv(J)

        Bt .= dNdX'
        G  = Bt*Ut # temperature gradient

        Δut = N'*dUt # interpolation to the integ. point

        q = update_state!(elem.mat, ip.state, Δut, G, Δt)

        #@showm q
        #error()

        coef  = elem.props.ρ*elem.props.cv
        coef *= detJ*ip.w*th
        dFt  .-= coef*N*Δut

        coef = Δt*detJ*ip.w*th
        @mul dFt += coef*Bt'*q
    end

    return dFt, map_t, success()
end
