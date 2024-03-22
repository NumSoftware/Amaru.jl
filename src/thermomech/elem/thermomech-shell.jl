# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TMShell

TMShellparams = [
    FunInfo(:TMShell, "A shell finite element for thermomechanical analyses"),
    KwArgInfo(:rho, "Density", 0.0, cond=:(rho>0)),
    KwArgInfo(:gamma, "Specific weight", 0.0, cond=:(gamma>=0)),
    KwArgInfo(:thickness, "Thickness", cond=:(thickness>0)),
    KwArgInfo(:alpha_s, "Shear correction coef.", 5/6, cond=:(alpha_s>0))
]
@doc docstring(TMShellparams) TMShell

struct TMShellProps<:ElemProperties
    ρ::Float64
    γ::Float64
    αs::Float64
    th::Float64

    function TMShellProps(; args...)
        args = checkargs(args, TMShellparams)

        return new(args.rho, args.gamma, args.alpha_s, args.thickness)
    end
end


mutable struct TMShell<:Thermomech
    env   ::ModelEnv
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::TMShellProps
    active::Bool
    linked_elems::Array{Element,1}
    Dlmn::Array{ SMatrix{3,3,Float64}, 1}

    function TMShell();
        return new()
    end
end


compat_shape_family(::Type{TMShell}) = BULKCELL
compat_elem_props(::Type{TMShell}) = TMShellProps


function elem_init(elem::TMShell)
    # check element dimension
    elem.shape.ndim==2 || throw(AmaruException("TMShell: Invalid element shape. Got $(elem.shape.name)"))

    # Compute nodal rotation matrices
    nnodes = length(elem.nodes)
    elem.Dlmn = Array{SMatrix{3,3,Float64}}(undef,nnodes)
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

        elem.Dlmn[i] = [V1'; V2'; V3' ]
    end

    return nothing
end


function setquadrature!(elem::TMShell, n::Int=0)
    # Set integration points
    if n in (8, 18)
        n = div(n,2)
    end
    ip2d = get_ip_coords(elem.shape, n)
    ip1d = get_ip_coords(LIN2, 2)
    n = size(ip2d,1)

    resize!(elem.ips, 2*n)
    for k in 1:2
        for i in 1:n
            R = [ ip2d[i].coord[1:2]; ip1d[k].coord[1] ]
            w = ip2d[i].w*ip1d[k].w
            j = (k-1)*n + i
            elem.ips[j] = Ip(R, w)
            elem.ips[j].id = j
            elem.ips[j].state = compat_state_type(typeof(elem.mat), typeof(elem), elem.env)(elem.env)
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

        dNdR = elem.shape.deriv(ip.R) # 3xn
        J = C'*dNdR
        No = normalize(cross(J[:,1], J[:,2]))
        ip.coord += elem.props.th/2*ip.R[3]*No
    end

end


function distributed_bc(elem::TMShell, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_boundary_forces(elem, facet, key, val)
end


function body_c(elem::TMShell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, key, val)
end


function elem_config_dofs(elem::TMShell)
    ndim = elem.env.ndim
    ndim==3 || error("MechShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
        add_dof(node, :ut, :ft)
    end
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


function calcS(elem::TMShell, αs::Float64)
    return @SMatrix [ 
        1.  0.  0.  0.  0.  0.
        0.  1.  0.  0.  0.  0.
        0.  0.  1.  0.  0.  0.
        0.  0.  0.  αs  0.  0.
        0.  0.  0.  0.  αs  0.
        0.  0.  0.  0.  0.  1. ]
end


function elem_map_u(elem::TMShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


function elem_map_t(elem::TMShell)
    return [ node.dofdict[:ut].eq_id for node in elem.nodes ]
end


function setB(elem::TMShell, ip::Ip, N::Vect, L::Matx, dNdX::Matx, Rrot::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.props.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        ζ = ip.R[3]
        Rrot[1:3,1:3] .= L
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
        @mul Bi = Bil*Rrot
        B[:, c+1:c+6] .= Bi
    end 
end


function elem_stiffness(elem::TMShell)
    th     = elem.props.th
    nnodes = length(elem.nodes)
    ndof   = 6
    nstr   = 6
    C      = getcoords(elem)
    K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rrot   = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)
    S      = calcS(elem, elem.props.αs)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′   = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        detJ′ = det(J′)
        @assert detJ′>0
        
        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)

        coef  = detJ′*ip.w
        D     = calcD(elem.mat, ip.state)
        K    += coef*B'*S*D*B
    end

    δ = 1e-5 # todo: add drilling control
    for i in 1:ndof*nnodes
        K[i,i] += δ
    end

    map = elem_map_u(elem)
    return K, map, map
end


# matrix C
function elem_coupling_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.props.th
    nnodes = length(elem.nodes)

    E = elem.mat.E
    ν = elem.mat.ν

    ndof = 6
    nstr = 6
    C    = getcoords(elem)
    B    = zeros(nstr, ndof*nnodes)
    L    = zeros(3,3)
    Rrot = zeros(5,ndof)
    Bil  = zeros(nstr,5)
    Bi   = zeros(nstr,ndof)
    C    = getcoords(elem)
    Cut  = zeros(ndof*nnodes, nnodes) # u-t coupling matrix
    m    = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ]



    for ip in elem.ips
        #elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)
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
        α = calc_α(elem.mat, ip.state.ut)
        β = E*α/(1-ν)

        coef  = β
        coef *= detJ′*ip.w
        mN    = m*N'
        @mul Cut -= coef*B'*mN
    end
    
    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

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
        @mul KBt = K*Bt
        @mul H -= coef*Bt'*KBt
    end

    map = elem_map_t(elem)
    return H, map, map
end


function elem_mass_matrix(elem::TMShell)
    th     = elem.props.th
    nnodes = length(elem.nodes)
    
    C = getcoords(elem)
    M = zeros(nnodes, nnodes)
    L = zeros(3,3)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR

        set_rot_x_xp(elem, J2D, L)
        J′    = [ L*J2D [ 0,0,th/2]  ]
        detJ′ = det(J′)
        @assert detJ′>0

        # compute Cut
        cv = calc_cv(elem.mat, ip.state.ut)
        coef  = elem.props.ρ*cv
        coef *= detJ′*ip.w
        M    -= coef*N*N'
    end

    map = elem_map_t(elem)

    return M, map, map
end


function update_elem!(elem::TMShell, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.props.th
    T0k    = elem.env.ana.T0 + 273.15
    nnodes = length(elem.nodes)
    ndof   = 6
    C      = getcoords(elem)

    E = elem.mat.E
    ν = elem.mat.ν
    ρ = elem.props.ρ

    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    dU  = DU[map_u]
    dUt = DU[map_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1
    m   = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    dF  = zeros(length(dU))
    B   = zeros(6, ndof*nnodes)
    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    L    = zeros(3,3)
    Rrot = zeros(5,ndof) 
    Bil  = zeros(6,5)
    Bi   = zeros(6,ndof)
    Δε   = zeros(6)
    S    = calcS(elem, elem.props.αs)


    for ip in elem.ips
        # compute Bu and Bt matrices
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′
        detJ′ = det(J′)

        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)

        # compute Δε
        @mul Δε = B*dU

        # compute Δut
        Δut = N'*dUt # interpolation to the integ. point

        # compute thermal gradient G
        Bt .= dNdX′'
        G  = Bt*Ut

        # internal force dF
        Δσ, q, status = update_state!(elem.mat, ip.state, Δε, Δut, G, Δt)
        failed(status) && return [dF; dFt], [map_u; map_t], status

        α = calc_α(elem.mat, ip.state.ut)
        β = E*α/(1-ν)

        Δσ -= β*Δut*m # get total stress
    
        coef = detJ′*ip.w
        dF += coef*B'*S*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0k
        coef *= detJ′*ip.w
        dFt  -= coef*N

        cv = calc_cv(elem.mat, ip.state.ut)
        coef  = ρ*cv
        coef *= detJ′*ip.w
        dFt  -= coef*N*Δut

        coef  = Δt
        coef *= detJ′*ip.w
        dFt += coef*Bt'*q

    end

    return [dF; dFt], [map_u; map_t], success()
end
