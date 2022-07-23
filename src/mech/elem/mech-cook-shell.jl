# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    CookShell
A bulk finite element for mechanical equilibrium analyses.
"""
mutable struct CookShell<:Mechanical
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv
    Dlmn::Array{ Array{Float64,2}, 1}

    function CookShell();
        return new()
    end
end

matching_shape_family(::Type{CookShell}) = SOLID_CELL


function elem_init(elem::CookShell)

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

function setquadrature!(elem::CookShell, n::Int=0)

    # n = 9
    ip2d = get_ip_coords(elem.shape, n)
    ip1d = get_ip_coords(LIN2, 2)
    n = size(ip2d,1)

    resize!(elem.ips, 2*n)
    for k in 1:2
        for i=1:n
            R = [ ip2d[i,1:2]; ip1d[k,1] ]
            w = ip2d[i,4]*ip1d[k,4]
            j = (k-1)*n + i
            elem.ips[j] = Ip(R, w)
            elem.ips[j].id = j
            elem.ips[j].state = ip_state_type(elem.mat)(elem.env)
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


function distributed_bc(elem::CookShell, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.mat.th
    suitable_keys = (:tx, :ty, :tz, :tn)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    # @show facet
    # @show elem
    # error()

    target = facet!==nothing ? facet : elem
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

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
        X = C'*N

        x, y, z = X
        vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        if key == :tx
            Q = [vip, 0.0, 0.0]
        elseif key == :ty
            Q = [0.0, vip, 0.0]
        elseif key == :tz
            Q = [0.0, 0.0, vip]
        elseif key == :tn
            n = cross(J[:,1], J[:,2])
            Q = vip*normalize(n)
        end

        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


# Rotation Matrix
function set_rot_x_xp(elem::CookShell, J::Matx, R::Matx)
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


function setB(elem::CookShell, ip::Ip, invJ::Matx, N::Vect, dNdX::Matx, Rrot::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.mat.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6

    for i in 1:nnodes
        ζ = ip.R[3]
        lx, ly, lz = elem.Dlmn[i][1,:] 
        mx, my, mz = elem.Dlmn[i][2,:]

        Rrot[4,4] = lx;  Rrot[4,5] = ly;  Rrot[4,6] = lz;
        Rrot[5,4] = mx;  Rrot[5,5] = my;  Rrot[5,6] = mz;

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        dNdz = dNdX[i,3]
        dζdx = invJ[3,1]
        dζdy = invJ[3,2]
        dζdz = invJ[3,3]
        dNζdx = dNdx*ζ + N[i]*dζdx
        dNζdy = dNdy*ζ + N[i]*dζdy
        dNζdz = dNdz*ζ + N[i]*dζdz

        Bil[1,1] = dNdx;                                                Bil[1,4] = -dNζdx*th/2*mx;                       Bil[1,5] = dNζdx*th/2*lx
                             Bil[2,2] = dNdy;                           Bil[2,4] = -dNζdy*th/2*my;                       Bil[2,5] = dNζdy*th/2*ly
                                                  Bil[3,3] = dNdz;      Bil[3,4] = -dNζdz*th/2*mz;                       Bil[3,5] = dNζdz*th/2*lz
                             Bil[4,2] = dNdz/SR2; Bil[4,3] = dNdy/SR2;  Bil[4,4] = -1/SR2*(dNζdz*th/2*my+dNζdy*th/2*mz);  Bil[4,5] = 1/SR2*(dNζdz*th/2*ly+dNζdy*th/2*lz)
        Bil[5,1] = dNdz/SR2;                      Bil[5,3] = dNdx/SR2;  Bil[5,4] = -1/SR2*(dNζdz*th/2*mx+dNζdx*th/2*mz);  Bil[5,5] = 1/SR2*(dNζdz*th/2*lx+dNζdx*th/2*lz)
        Bil[6,1] = dNdy/SR2; Bil[6,2] = dNdx/SR2;                       Bil[6,4] = -1/SR2*(dNζdy*th/2*mx+dNζdx*th/2*my);  Bil[6,5] = 1/SR2*(dNζdy*th/2*lx+dNζdx*th/2*ly)
        
        c = (i-1)*ndof
        @gemm Bi = Bil*Rrot
        B[:, c+1:c+6] .= Bi

    end 

    # @showm round.(B[:,1:5], digits=5)
    # error()
end


function elem_config_dofs(elem::CookShell)
    ndim = elem.env.ndim
    ndim in (1,2) && error("CookShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end

function elem_map(elem::CookShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


function elem_stiffness(elem::CookShell)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th = elem.mat.th
    ndof = 6

    C = getcoords(elem)
    K = zeros(ndof*nnodes, ndof*nnodes)
    B = zeros(6, ndof*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    L = zeros(3,3)
    T = zeros(6,6)
    Rrot = Array{Float64}(I,5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)

    Dn = [ elem.Dlmn[i][3,j] for i in 1:nnodes, j in 1:3 ] # nx3

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        J = [ C'*dNdR + th/2*ip.R[3]*Dn'*dNdR   th/2*Dn'*N ] # 3x3

        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)

        set_tensor_rot!(L, T)
        
        dNdR = [ dNdR zeros(nnodes) ]
        invJ = inv(J)
        dNdX = dNdR*invJ

        D = calcD(elem.mat, ip.state)
        setB(elem, ip, invJ, N, dNdX, Rrot, Bil, Bi, B)
        detJ = det(J)

        @assert detJ>0
        coef = detJ*ip.w

        K += coef*B'*T'*D*T*B
    end

    # δ = 1e-10
    δ = 1e-7
    # δ = 0.0
    for i in 1:ndof*nnodes
        K[i,i] += δ
    end

    map = elem_map(elem)
    return K, map, map
end

function elem_update!(elem::CookShell, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    th = elem.mat.th
    ndof = 6

    map = elem_map(elem)
    dU = U[map]
    dF = zeros(length(dU))

    C = getcoords(elem)
    B = zeros(6, ndof*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    L = zeros(3,3)
    T = zeros(6,6)
    Rrot = Array{Float64}(I,5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)
    Δε = zeros(6)

    Dn = [ elem.Dlmn[i][3,j] for i in 1:nnodes, j in 1:3 ] # nx3


    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J = [ C'*dNdR + th/2*ip.R[3]*Dn'*dNdR   th/2*Dn'*N ] # 3x3
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        set_tensor_rot!(L, T)

        dNdR = [ dNdR zeros(nnodes) ]
        invJ = inv(J)
        dNdX = dNdR*invJ
        setB(elem, ip, invJ, N, dNdX, Rrot, Bil, Bi, B)
        Δε = T*B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return failure("MechSolid: Error at integration point $(ip.id)")

        detJ = det(J)
        coef = detJ*ip.w
        dF += coef*B'*T'*Δσ

    end

    F[map] += dF
    return success()
end
