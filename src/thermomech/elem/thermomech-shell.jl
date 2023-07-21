# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


export TMShell

struct TMShellProps<:ElemProperties
    ρ::Float64
    γ::Float64
    cv::Float64
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C
    th::Float64

    function TMShellProps(; props...)
        names = (rho="Density", gamma="Specific weight", thickness="Thickness", cv="Heat capacity", alpha="Thermal expansion coefficient")
        required = (:thickness, :cv, :alpha)
        @checkmissing props required names

        default = (rho=0.0, gamma=0.0)
        props  = merge(default, props)

        rho       = props.rho
        gamma     = props.gamma
        cv        = props.cv
        thickness = props.thickness
        alpha = props.alpha

        @check rho>=0.0
        @check gamma>=0.0
        @check cv>0.0
        @check thickness>0.0
        @check 0<=alpha<=1

        return new(rho, gamma, cv, alpha, thickness)
    end    
end


mutable struct TMShell<:Thermomech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::TMShellProps
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv
    Dlmn::Array{ Array{Float64,2}, 1}

    function TMShell();
        return new()
    end
end


matching_shape_family(::Type{TMShell}) = BULKCELL
matching_elem_props(::Type{TMShell}) = TMShellProps


function elem_init(elem::TMShell)
    
    nnodes = length(elem.nodes)
    Dlmn = Array{Float64,2}[]
    C = getcoords(elem)

    for i in 1:nnodes
        Ri = elem.shape.nat_coords[i,:]
        dNdR = elem.shape.deriv(Ri)
        J = C'*dNdR

        V1 = J[:,1]
        V2 = J[:,2]
        V3 = cross(V1, V2)
        V2 = cross(V3, V1)
        normalize!(V1)
        normalize!(V2)
        normalize!(V3)

        push!(Dlmn, [V1'; V2'; V3' ])
    end
    elem.Dlmn = Dlmn

    return nothing
end


function setquadrature!(elem::TMShell, n::Int=0)

    if n in (8, 18)
        n = div(n,2)
    end
    ip2d = get_ip_coords(elem.shape, n)
    ip1d = get_ip_coords(LIN2, 2)
    n = size(ip2d,1)

    resize!(elem.ips, 2*n)
    for k in 1:2
        for i in 1:n
            R = [ ip2d[i,1:2]; ip1d[k,1] ]
            w = ip2d[i,4]*ip1d[k,4]
            j = (k-1)*n + i
            elem.ips[j] = Ip(R, w)
            elem.ips[j].id = j
            elem.ips[j].state = ip_state_type(typeof(elem.mat))(elem.env)
            elem.ips[j].owner = elem
        end
    end

    # finding ips global coordinates
    C     = getcoords(elem)
    shape = elem.shape

    for ip in elem.ips
        R = [ ip.R[1:2]; 0.0 ]
        N = shape.func(R)
        ip.coord = C'*N
    end

end

# DUVIDA!!!!!!
function distributed_bc(elem::TMShell, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
        return mech_shell_boundary_forces(elem, facet, key, val)
end


function body_c(elem::TMShell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, key, val)
end


function elem_config_dofs(elem::TMShell)
    ndim = elem.env.ndim
    ndim in (1,2) && error("TMShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
        add_dof(node, :ut, :ft) # VERIFICAR
    end
end

@inline function elem_map_u(elem::TMShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end

@inline function elem_map_t(elem::TMShell)
    return [ node.dofdict[:ut].eq_id for node in elem.nodes ]
end


# Rotation Matrix
function set_rot_x_xp(elem::TMShell, J::Matx, R::Matx)
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


function setB(elem::TMShell, ip::Ip, N::Vect, L::Matx, dNdX::Matx, Rrot::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.props.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        ζ = ip.R[3]
        Rrot[1:3,1:3] .= L
        #Rrot[4:5,4:6] .= L[1:2,:]
        #  Rrot[1:3,1:3] .= elem.Dlmn[i]
        Rrot[4:5,4:6] .= elem.Dlmn[i][1:2,:]
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        Ni = N[i]

        
        Bil[1,1] = dNdx;                                                                                Bil[1,5] = dNdx*ζ*th/2
                             Bil[2,2] = dNdy;                           Bil[2,4] = -dNdy*ζ*th/2
                                                  
                                                  Bil[4,3] = dNdy/SR2;  Bil[4,4] = -1/SR2*Ni
                                                  Bil[5,3] = dNdx/SR2;                                  Bil[5,5] = 1/SR2*Ni
        Bil[6,1] = dNdy/SR2; Bil[6,2] = dNdx/SR2;                       Bil[6,4] = -1/SR2*dNdx*ζ*th/2;  Bil[6,5] = 1/SR2*dNdy*ζ*th/2

        c = (i-1)*ndof
        @gemm Bi = Bil*Rrot
        B[:, c+1:c+6] .= Bi
    end 
end


function elem_stiffness(elem::TMShell)
    nnodes = length(elem.nodes)
    th     = elem.props.th
    ndof   = 6
    nstr   = 6
    C      = getcoords(elem)
    K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rrot   = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′   = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        D = calcD(elem.mat, ip.state, "shell")
        #@show D
        #error()

        detJ′ = det(J′)
        @assert detJ′>0
        
        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)
        coef = detJ′*ip.w
        K += coef*B'*D*B
    end

    δ = 1e-5
    for i in 1:ndof*nnodes
        K[i,i] += δ
    end

    map = elem_map_u(elem)
    #@show K
    #error()
    return K, map, map
end

#=
# DUVIDA !!!!!!!!!!!!!!
@inline function set_Bu1(elem::Element, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end
=#

# matrix C
function elem_coupling_matrix(elem::TMShell)
    #ndim   = elem.env.ndim
    #th     = elem.env.ana.thickness
    #nnodes = length(elem.nodes)
    #
    #C   = getcoords(elem)
    #Bu  = zeros(6, nnodes*ndim)
    #Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    #J    = Array{Float64}(undef, ndim, ndim)
    #dNdX = Array{Float64}(undef, nnodes, ndim)
    #m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    #β    = elem.mat.E*elem.props.α/(1-2*elem.mat.ν) # thermal stress modulus

    nnodes = length(elem.nodes)
    # 
    ndim   = elem.env.ndim

    th     = elem.props.th
    ndof   = 6
    nstr   = 6
    C      = getcoords(elem)
    #K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rrot   = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)

    C   = getcoords(elem)
    Cut = zeros(ndof*nnodes, nnodes) # u-t coupling matrix
    m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = elem.mat.E*elem.props.α/(1-2*elem.mat.ν) 

    for ip in elem.ips
        #elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        #dNdR = elem.shape.deriv(ip.R)
        #@gemm J = C'*dNdR
        #@gemm dNdX = dNdR*inv(J)
        #detJ = det(J)
        #detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′   = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′      
        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)

        detJ′ = det(J′)
        @assert detJ′>0
        # compute Cut
        coef  = β
        coef *= detJ′*ip.w
        mN   = m*N'
        @gemm Cut -= coef*B'*mN
    end
    # map
    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    # keys = (:ux, :uy, :uz)[1:ndim]
    # map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    # mat_t = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nnodes] ]

    #@show Cut

    return Cut, map_u, map_t
end

# thermal conductivity
function elem_conductivity_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.props.th   # elem.env.ana.thickness
    nnodes = length(elem.nodes)
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    Bt     = zeros(ndim, nnodes)
    KBt    = zeros(ndim, nnodes)
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
        #coef = detJ′*ip.w*th
        @gemm KBt = K*Bt
        @gemm H -= coef*Bt'*KBt
    end

    map = elem_map_t(elem)
    #@show H
    return H, map, map
end


function elem_mass_matrix(elem::TMShell)
    th     = elem.props.th
    nnodes = length(elem.nodes)
    
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)

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
        coef  = elem.props.ρ*elem.props.cv
        coef *= detJ′*ip.w
        M    -= coef*N*N'
    end

    # map
    map = elem_map_t(elem)

    #@show "HIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    #@show M
    return M, map, map
    
end

#=
function elem_internal_forces(elem::TMShell, F::Array{Float64,1}, DU::Array{Float64,1})
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
    m = [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]  tI
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Jp  = Array{Float64}(undef, ndim, nbnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    dUt = DU[mat_t] # nodal temperature increments
    for ip in elem.ips
        elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)
        # compute Bu matrix and Bt
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @gemm dNdX = dNdR*inv(J)
        set_Bu(elem, ip, dNdX, Bu)
        dNtdR = elem.shape.deriv(ip.R)
        Jp = dNtdR*Ct
        @gemm dNtdX = inv(Jp)*dNtdR
        Bt = dNtdX
        # compute N
        # internal force
        ut   = ip.state.ut + 273
        β   = elem.mat.E*elem.props.α/(1-2*elem.mat.ν)
        σ    = ip.state.σ - β*ut*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ
        # internal volumes dFt
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = β*detJ*ip.w*th
        dFt  -= coef*N*εvol
        coef = detJ*ip.w*elem.props.ρ*elem.props.cv*th/T0
        dFt -= coef*N*ut
        QQ   = ip.state.QQ
        coef = detJ*ip.w*th/T0
        @gemv dFt += coef*Bt'*QQ
    end
    F[map_u] += dF
    F[mat_t] += dFt
end
=#


function update_elem!(elem::TMShell, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.props.th
    T0k     = elem.env.ana.T0 + 273.15
    nnodes = length(elem.nodes)
    ndof = 6
    C      = getcoords(elem)

    E = elem.mat.E
    α = elem.props.α
    ρ = elem.props.ρ
    nu = elem.mat.ν
    cv = elem.props.cv
    β = E*α/(1-2*nu)
    
    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    dU = DU[map_u]
    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1
    m   = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]  #

    dF = zeros(length(dU))
    B = zeros(6, ndof*nnodes)
    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    L = zeros(3,3)
    Rrot = zeros(5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)
    Δε = zeros(6)

    for ip in elem.ips
        #elem.env.ana.stressmodel=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)

        # compute Δε
        @gemv Δε = B*dU

        # compute Δut
        Δut = N'*dUt # interpolation to the integ. point

        #@show size(Ut)
        #@show size(Bt)
        # compute thermal gradient G
        Bt .= dNdX′'
        G  = Bt*Ut

        # internal force dF
        Δσ, q, status = update_state!(elem.mat, ip.state, Δε, Δut, G, Δt, "shell")
        failed(status) && return [dF; dFt], [map_u; map_t], status
        #error()
        Δσ -= β*Δut*m # get total stress
        #@showm Δσ
        #error()
    
        detJ′ = det(J′)
        #coef = detJ′*ip.w*th
        coef = detJ′*ip.w
        @gemv dF += coef*B'*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0k
        coef *= detJ′*ip.w
        dFt  -= coef*N

        coef  = ρ*cv
        coef *= detJ′*ip.w
        dFt  -= coef*N*Δut

        coef  = Δt
        coef *= detJ′*ip.w
        @gemv dFt += coef*Bt'*q

    end
    #@show "HIIIIIIIIIIIIIIIIIIIIII UPDATED"
    return [dF; dFt], [map_u; map_t], success()
end