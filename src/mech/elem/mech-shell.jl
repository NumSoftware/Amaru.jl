# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechShell


MechShell_params = [
    FunInfo(:MechShell, "A Reissner-Mindlin curved shell element with drilling control"),
    KwArgInfo(:rho, "Density", 0.0, cond=:(rho>=0.0)),
    KwArgInfo(:gamma, "Specific weight", 0.0, cond=:(gamma>=0.0)),
    KwArgInfo(:thickness, "Thickness", cond=:(thickness>0.0)),
    KwArgInfo(:alpha_s, "Shear correction coef.", 5/6, cond=:(alpha_s>0)),
    KwArgInfo(:kappa, "Drilling penalty coef.", 1e-8)
]
@doc docstring(MechShell_params) MechShell

struct MechShellProps<:ElemProperties
    ρ::Float64
    γ::Float64
    αs::Float64
    th::Float64
    kappa::Float64

    function MechShellProps(; args...)
        args = checkargs(args, MechShell_params)

        return new(args.rho, args.gamma, args.alpha_s, args.thickness, args.kappa)
    end
end


mutable struct MechShell<:Mech
    ctx::Context
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::MechShellProps
    active::Bool
    linked_elems::Array{Element,1}
    Dlmn::Array{ SMatrix{3,3,Float64}, 1}

    function MechShell()
        return new()
    end
end 


compat_shape_family(::Type{MechShell}) = BULKCELL
compat_elem_props(::Type{MechShell}) = MechShellProps


function elem_init(elem::MechShell)
    # check element dimension
    elem.shape.ndim==2 || throw(AmaruException("MechShell: Invalid element shape. Got $(elem.shape.name)"))

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


function setquadrature!(elem::MechShell, n::Int=0)
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
            elem.ips[j].state = compat_state_type(typeof(elem.mat), typeof(elem), elem.ctx)(elem.ctx)
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


function distributed_bc(elem::MechShell, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_boundary_forces(elem, facet, t, key, val)
end


function body_c(elem::MechShell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, 0.0, key, val)
end


function elem_config_dofs(elem::MechShell)
    ndim = elem.ctx.ndim
    ndim==3 || error("MechShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end



# Rotation Matrix
function set_rot_x_xp(elem::MechShell, J::Matx, R::Matx)
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


function calcS(elem::MechShell, αs::Float64)
    return @SMatrix [ 
        1.  0.  0.  0.  0.  0.
        0.  1.  0.  0.  0.  0.
        0.  0.  1.  0.  0.  0.
        0.  0.  0.  αs  0.  0.
        0.  0.  0.  0.  αs  0.
        0.  0.  0.  0.  0.  1. ]
end


function elem_map(elem::MechShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


function setB(elem::MechShell, ip::Ip, N::Vect, L::Matx, dNdX::Matx, Rθ::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.props.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        ζ = ip.R[3]
        Rθ[1:3,1:3] .= L
        # Rθ[1:3,1:3] .= elem.Dlmn[i]
        Rθ[4:5,4:6] .= elem.Dlmn[i][1:2,:]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        Ni = N[i]

        Bil[1,1] = dNdx;                                                                                Bil[1,5] = dNdx*ζ*th/2
                             Bil[2,2] = dNdy;                           Bil[2,4] = -dNdy*ζ*th/2
                                                  
                                                  Bil[4,3] = dNdy/SR2;  Bil[4,4] = -1/SR2*Ni
                                                  Bil[5,3] = dNdx/SR2;                                  Bil[5,5] = 1/SR2*Ni
        Bil[6,1] = dNdy/SR2; Bil[6,2] = dNdx/SR2;                       Bil[6,4] = -1/SR2*dNdx*ζ*th/2;  Bil[6,5] = 1/SR2*dNdy*ζ*th/2

        c = (i-1)*ndof
        @mul Bi = Bil*Rθ
        B[:, c+1:c+6] .= Bi
    end 
end

function setB_dr(elem::MechShell, N::Vect, L::Matx, dNdX::Matx, Rθ_dr::Matx, Bil_dr::Matx, Bi_dr::Matx, B_dr::Matx)
    nnodes = size(dNdX,1)

    ndof = 6
    for i in 1:nnodes

        Rθ_dr[1:2,1:3] .= L[1:2,1:3]
        Rθ_dr[3,4:6] .= elem.Dlmn[i][3,:]
        
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        Ni   = N[i]
        

        Bil_dr[1,1] = +0.5*dNdy
        Bil_dr[1,2] = -0.5*dNdx
        Bil_dr[1,3] = Ni                                                          
      
        c = (i-1)*ndof
        @mul Bi_dr = Bil_dr*Rθ_dr
        
        B_dr[:, c+1:c+6] .= Bi_dr
    end 
end

function setNN(elem::MechShell, ip::Ip, N::Vect, NNil::Matx, NNi::Matx, L::Matx, Rθ::Matx, NN::Matx)
    nnodes = length(N)
    ndof = 6
    th = elem.props.th
    ζ = ip.R[3]

    for i in 1:nnodes
        
        Rθ[1:3,1:3] .= L
        Rθ[4:5,4:6] .= L[1:2,:]

        NNil[1,1] = N[i]
        NNil[2,2] = N[i]
        NNil[3,3] = N[i]
        NNil[1,5] = th/2*ζ*N[i]
        NNil[2,4] = -th/2*ζ*N[i]

        c = (i-1)*ndof
        @mul NNi = NNil*Rθ

        NN[:, c+1:c+6] .= NNi
    end 
end


function elem_stiffness(elem::MechShell)
    nnodes = length(elem.nodes)
    th     = elem.props.th
    kappa  = elem.props.kappa
    ndof   = 6
    nstr   = 6
    C      = getcoords(elem)
    K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rθ     = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)
    S      = calcS(elem, elem.props.αs)

    B_dr   = zeros(1, ndof*nnodes)
    Bil_dr = zeros(1,3)
    Bi_dr  = zeros(1,ndof)
    Rθ_dr  = zeros(3,ndof)
    m      = div(length(elem.ips), 2) # half the number of integration points

    for (i,ip) in enumerate(elem.ips)
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
        
        setB(elem, ip, N, L, dNdX′, Rθ, Bil, Bi, B)

        E  = elem.mat.E
        nu = elem.mat.ν
        G  = E/(2*(1+nu))

        coef  = detJ′*ip.w
        D     = calcD(elem.mat, ip.state)
        K    += coef*B'*S*D*B
        
        if i<=m # drilling stiffness (area integration)
            setB_dr(elem, N, L, dNdX′, Rθ_dr, Bil_dr, Bi_dr, B_dr)       
            coef = kappa*G*norm2(J2D)*th*ip.w
            @mul K += coef*B_dr'*B_dr
        end

    end

#=
    δ = 1e-7
    for i in 1:ndof*nnodes
        K[i,i] += δ
    end
=#
    map = elem_map(elem)
    return K, map, map
end


function elem_mass(elem::MechShell)
        nnodes = length(elem.nodes)
        th     = elem.props.th
        ndof   = 6 #6
        ρ      = elem.props.ρ
        C      = getcoords(elem)
        M      = zeros(nnodes*ndof, nnodes*ndof)
        L      = zeros(3,3)
        Rθ   = zeros(5,ndof)
        
        NN     = zeros(3, nnodes*ndof)
        NNil    = zeros(3,5)
        NNi     = zeros(3,ndof)
    
        for ip in elem.ips
            # compute N matrix
            N    = elem.shape.func(ip.R)
            dNdR = elem.shape.deriv(ip.R)
            J2D  = C'*dNdR
            set_rot_x_xp(elem, J2D, L)
            J′   = [ L*J2D [ 0,0,th/2]  ] 
                                  
            detJ′ = det(J′)
            @assert detJ′>0

            setNN(elem, ip, N, NNil, NNi, L, Rθ, NN)

            # compute M
            coef = ρ*detJ′*ip.w
            @mul M += coef*NN'*NN
        end

        map = elem_map(elem)
        return M, map, map
end


function update_elem!(elem::MechShell, U::Array{Float64,1}, dt::Float64)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    th   = elem.props.th
    ndof = 6

    map = elem_map(elem)
    dU  = U[map]
    dF  = zeros(length(dU))

    C = getcoords(elem)
    B = zeros(6, ndof*nnodes)

    L    = zeros(3,3)
    Rθ = zeros(5,ndof)
    Bil  = zeros(6,5)
    Bi   = zeros(6,ndof)
    Δε   = zeros(6)
    S    = calcS(elem, elem.props.αs)


    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        setB(elem, ip, N, L, dNdX′, Rθ, Bil, Bi, B)
        Δε = B*dU
        Δσ, status = update_state!(elem.mat, ip.state, Δε)
        failed(status) && return dF, map, failure("MechShell: Error at integration point $(ip.id)")

        detJ′ = det(J′)
        coef  = detJ′*ip.w
        dF   += coef*B'*S*Δσ

    end

     return dF, map, success()
end
