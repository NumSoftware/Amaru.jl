# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export AcousticFluid

AcousticFluid_params = [
    FunInfo(:AcousticFluid, "A finite element for wave analyses"),
    KwArgInfo(:rho, "Density", 0.0, cond=:(rho>=0.0)),
    KwArgInfo(:c, "Speed of sound", 0.0, cond=:(c>0.0)),
]
@doc docstring(AcousticFluid_params) AcousticFluid


struct AcousticFluidProps<:ElemProperties
    ρ::Float64
    c::Float64

    function AcousticFluidProps(; kwargs...)
        args = checkargs(kwargs, AcousticFluid_params)
        return new(args.rho, args.c)
    end    
end


mutable struct AcousticFluid<:AcousticMech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::AcousticFluidProps
    active::Bool
    linked_elems::Array{Element,1}
    ctx::Context

    function AcousticFluid()
        return new()
    end
end

compat_shape_family(::Type{AcousticFluid}) = BULKCELL
compat_elem_props(::Type{AcousticFluid}) = AcousticFluidProps


function elem_config_dofs(elem::AcousticFluid)
    for node in elem.nodes
        add_dof(node, :up, :fq) # up: pressure; fq: mass acceleration (e.g. kg/s2)
    end
end


function elem_init(elem::AcousticFluid)
    
end


function distributed_bc(elem::AcousticFluid, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return acoustic_mech_bc(elem, facet, t, key, val)
end


# acoustic fluid stiffness
function elem_acoustic_stiffness(elem::AcousticFluid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    K      = zeros(nnodes, nnodes)
    dNdX   = zeros(nnodes, ndim) # cartesian derivatives
    J      = Array{Float64}(undef, ndim, ndim) # Jacobian

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR*inv(J)
        Bp = dNdX'

        coef = detJ*ip.w*th

        # @mul 
        K += coef*Bp'*Bp
    end

    # map
    map = [ node.dofdict[:up].eq_id for node in elem.nodes ]

    return K, map, map
end

         
function elem_acoustic_mass(elem::AcousticFluid)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)
    J      = Array{Float64}(undef, ndim, ndim)
    c      = elem.props.c # sound speed

    for ip in elem.ips
        elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = detJ*ip.w*th/c^2
        M   += coef*N*N'
    end

    # map
    map = [ node.dofdict[:up].eq_id for node in elem.nodes  ]

    return M, map, map
end


# TODO
function update_elem!(elem::AcousticFluid, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    th     = elem.ctx.thickness

    map_p  = [ node.dofdict[:up].eq_id for node in elem.nodes ]

    C   = getcoords(elem)

    dP = DU[map_p] # nodal pore-pressure increments
    P  = [ node.dofdict[:up].vals[:up] for node in elem.nodes ]
    P += dP # nodal pore-pressure at step n+1
    A  = [ node.dofdict[:up].vals[:ap] for node in elem.nodes ]

    K, m, m = elem_acoustic_stiffness(elem)
    M, m, m = elem_acoustic_mass(elem)

    # dF = K*P + M*A
    # dF = K*dP + M*A
    dF = K*dP #+ M*A

    # dFw = zeros(nnodes)
    # Bw  = zeros(ndim, nnodes)

    # J  = Array{Float64}(undef, ndim, ndim)
    # dNdX = Array{Float64}(undef, nnodes, ndim)

    # for ip in elem.ips
    #     elem.ctx.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

    #     # compute Bu matrix
    #     N    = elem.shape.func(ip.R)
    #     dNdR = elem.shape.deriv(ip.R)
    #     @mul J = C'*dNdR
    #     detJ = det(J)
    #     detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
    #     @mul dNdX = dNdR*inv(J) # Bw = dNdX'

    #     Δuw = N'*dUp # interpolation to the integ. point

    #     V = update_state!(elem.mat, ip.state, Δuw, G, Δt)

    #     coef  = elem.mat.S
    #     coef *= detJ*ip.w*th
    #     dFw  -= coef*N*Δuw

    #     coef = Δt*detJ*ip.w*th
    #     @mul dFw += coef*dNdX*V
    # end

    return dF, map_p, success()
    # DF[map_p] += dF
    # return success()
end

