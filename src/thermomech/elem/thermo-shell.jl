# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ThermoShell

struct ThermoShellProps<:ElemProperties
    ρ::Float64
    th::Float64

    function ThermoShellProps(; args...)
        args = checkargs(args, arg_rules(ThermoShellProps))

        return new(args.rho, args.thickness)
    end
end

arg_rules(::Type{ThermoShellProps}) = 
[
    @arginfo rho rho>=0.0 "Density"
    @arginfo thickness thickness>0.0 "Thickness"
]

mutable struct ThermoShell<:Thermomech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::ThermoShellProps
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv

    function ThermoShell();
        return new()
    end
end

compat_shape_family(::Type{ThermoShell}) = BULKCELL
compat_elem_props(::Type{ThermoShell}) = ThermoShellProps


function elem_config_dofs(elem::ThermoShell)
    for node in elem.nodes
        add_dof(node, :ut, :ft)
    end
end


function elem_init(elem::ThermoShell)
    nothing
end

@inline function elem_map_t(elem::ThermoShell)
    return [ node.dofdict[:ut].eq_id for node in elem.nodes ]
end


function distributed_bc(elem::ThermoShell, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.ana.thickness
    suitable_keys = (:tq,)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a ThermoSolid element")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

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
            if elem.env.ana.stressmodel=="axisymmetric"
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

# Rotation Matrix
function set_rot_x_xp(elem::ThermoShell, J::Matx, R::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    V2 = cross(V3, V1)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    R[1,:] .= V1
    R[2,:] .= V2
    R[3,:] .= V3
end

# thermal conductivity
function elem_conductivity_matrix(elem::ThermoShell)
    ndim   = elem.env.ndim
    th     = elem.props.th   # elem.env.ana.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    Bt     = zeros(ndim, nnodes)
    L      = zeros(3,3)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR

        set_rot_x_xp(elem, J2D, L)
        J′    = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        Bt   .= dNdX′'

        K = calcK(elem.mat, ip.state)
        detJ′ = det(J′)
        @assert detJ′>0
        
        coef = detJ′*ip.w
        H -= coef*Bt'*K*Bt

    end

    map = elem_map_t(elem)
    return H, map, map
end


function elem_mass_matrix(elem::ThermoShell)
    ndim   = elem.env.ndim
    th     = elem.props.th
    nnodes = length(elem.nodes)
    
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    L      = zeros(3,3)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR

        set_rot_x_xp(elem, J2D, L)
        J′    = [ L*J2D [ 0,0,th/2]  ]
        detJ′ = det(J′)
        @assert detJ′>0

        # compute Cut
        coef  = elem.props.ρ*elem.mat.cv
        coef *= detJ′*ip.w
        M    -= coef*N*N'

    end

    # map
    map = elem_map_t(elem)

    return M, map, map
end

#=
function elem_internal_forces(elem::ThermoShell, F::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = thickness   # elem.env.ana.thickness
    nnodes = length(elem.nodes)
    
    C   = getcoords(elem)
    T0     = elem.env.T0 + 273.15
    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]
    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nbnodes)
    Bt  = zeros(ndim, nbnodes)
    m = [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]  I2
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Jp  = Array{Float64}(undef, ndim, nbnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    dUt = DU[mat_t] # nodal temperature increments
    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)
        # compute Bu matrix and Bt
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @mul dNdX = dNdR*inv(J)
        set_Bu(elem, ip, dNdX, Bu)
        dNtdR = elem.shape.deriv(ip.R)
        Jp = dNtdR*Ct
        @mul dNtdX = inv(Jp)*dNtdR
        Bt = dNtdX
        # compute N
        # internal force
        ut   = ip.state.ut + 273
        β   = elem.mat.E*elem.mat.α/(1-2*elem.mat.ν)
        σ    = ip.state.σ - β*ut*m # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*σ
        # internal volumes dFt
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = β*detJ*ip.w*th
        dFt  -= coef*N*εvol
        coef = detJ*ip.w*elem.props.ρ*elem.mat.cv*th/T0
        dFt -= coef*N*ut
        QQ   = ip.state.QQ
        coef = detJ*ip.w*th/T0
        @mul dFt += coef*Bt'*QQ
    end
    F[map_u] += dF
    F[mat_t] += dFt
end
=#


function update_elem!(elem::ThermoShell, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th     = elem.props.th
    ndof = 6

    map_t = elem_map_t(elem)

    C      = getcoords(elem)

    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1

    Bt  = zeros(ndim, nnodes)
    dFt = zeros(nnodes)

    L = zeros(3,3)
    #Rrot = zeros(5,ndof)

    for ip in elem.ips
        #elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        detJ′ = det(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        Bt .= dNdX′'
        G  = Bt*Ut # temperature gradient

        # compute Δut
        Δut = N'*dUt # interpolation to the integ. point

        q = update_state!(elem.mat, ip.state, Δut, G, Δt)

        coef  = elem.props.ρ*elem.mat.cv
        coef *= detJ′*ip.w
        dFt  -= coef*N*Δut

        coef = Δt*detJ′*ip.w
        dFt += coef*Bt'*q

    end
    return dFt, map_t, success()
end
