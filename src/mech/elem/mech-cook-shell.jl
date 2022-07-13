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
    # elem.shape==QUAD8 || error("elem_init: CookShell only works with shape QUAD8.")

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

        push!(Dlmn, [V1 V2 V3 ])
    end
    elem.Dlmn = Dlmn
   
    return nothing
end

function setquadrature!(elem::CookShell, n::Int=0)

    # if !(n in keys(elem.shape.quadrature))
    #     alert("setquadrature!: cannot set $n integration points for shape $(elem.shape.name)")
    #     return
    # end

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


# Rotation Matrix
function set_dir_matrix(elem::CookShell, J::Matx, dir::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    V2 = cross(V3, V1)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    dir[:,1] .= V1
    dir[:,2] .= V2
    dir[:,3] .= V3
end


function setB(elem::CookShell, ip::Ip, invJ::Matx, N::Vect, dNdX::Matx, Rrot::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    t = elem.mat.t
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6

    for i in 1:nnodes
        ζ = ip.R[3]
        lx, ly, lz = elem.Dlmn[i][:,1] 
        mx, my, mz = elem.Dlmn[i][:,2]

        Rrot[4,4] = lx;  Rrot[4,5] = ly;  Rrot[4,6] = lz;
        Rrot[5,4] = mx;  Rrot[5,5] = my;  Rrot[5,6] = mz;

        # @showm Rrot
        # error()
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        dNdz = dNdX[i,3]
        # dζdx = invJ[1,3]
        # dζdy = invJ[2,3]
        # dζdz = invJ[3,3]
        dζdx = invJ[3,1]
        dζdy = invJ[3,2]
        dζdz = invJ[3,3]
        dNζdx = dNdx*ζ + N[i]*dζdx
        dNζdy = dNdy*ζ + N[i]*dζdy
        dNζdz = dNdz*ζ + N[i]*dζdz
        # @showm N
        # @show dNζdx
        # @show dNζdy
        # @show dNζdz
        # @show mx
        # @show my
        # @show mz
        # @show lx
        # @show ly
        # @show lz

        Bil[1,1] = dNdx;                                                Bil[1,4] = -dNζdx*t/2*mx;                       Bil[1,5] = dNζdx*t/2*lx
                             Bil[2,2] = dNdy;                           Bil[2,4] = -dNζdy*t/2*my;                       Bil[2,5] = dNζdy*t/2*ly
                                                  Bil[3,3] = dNdz;      Bil[3,4] = -dNζdz*t/2*mz;                       Bil[3,5] = dNζdz*t/2*lz
                             Bil[4,2] = dNdz/SR2; Bil[4,3] = dNdy/SR2;  Bil[4,4] = -1/SR2*(dNζdz*t/2*my+dNζdy*t/2*mz);  Bil[4,5] = 1/SR2*(dNζdz*t/2*ly+dNζdy*t/2*lz)
        Bil[5,1] = dNdz/SR2;                      Bil[5,3] = dNdx/SR2;  Bil[5,4] = -1/SR2*(dNζdz*t/2*mx+dNζdx*t/2*mz);  Bil[5,5] = 1/SR2*(dNζdz*t/2*lx+dNζdx*t/2*lz)
        Bil[6,1] = dNdy/SR2; Bil[6,2] = dNdx/SR2;                       Bil[6,4] = -1/SR2*(dNζdy*t/2*mx+dNζdx*t/2*my);  Bil[6,5] = 1/SR2*(dNζdy*t/2*lx+dNζdx*t/2*ly)
        
        # B[1,1+c] = dNdx;                                                B[1,4+c] = -dNζdx*t/2*mx;               B[1,5+c] = dNζdx*t/2*lx
        #                      B[2,2+c] = dNdy;                           B[2,4+c] = -dNζdy*t/2*my;               B[2,5+c] = dNζdy*t/2*ly
        #                                           B[3,3+c] = dNdz;      B[3,4+c] = -dNζdz*t/2*mz;               B[3,5+c] = dNζdz*t/2*lz
        #                      B[4,2+c] = dNdz/SR2; B[4,3+c] = dNdy/SR2;  B[4,4+c] = -1/SR2*(dNζdz*t/2*my+dNζdy*t/2*mz);  B[4,5+c] = 1/SR2*(dNζdz*t/2*ly+dNζdy*t/2*lz)
        # B[5,1+c] = dNdz/SR2;                      B[5,3+c] = dNdx/SR2;  B[5,4+c] = -1/SR2*(dNζdz*t/2*mx+dNζdx*t/2*mz);  B[5,5+c] = 1/SR2*(dNζdz*t/2*lx+dNζdx*t/2*lz)
        # B[6,1+c] = dNdy/SR2; B[6,2+c] = dNdx/SR2;                       B[6,4+c] = -1/SR2*(dNζdy*t/2*mx+dNζdx*t/2*my);  B[6,5+c] = 1/SR2*(dNζdy*t/2*lx+dNζdx*t/2*ly)

        # @showm invJ

        # Bil[4:6,:] .*= SR2
        # # @showm Bil
        # Bil = Bil[[1,2,3,6,5,4],:]
        # @showm Bil
        # error()

        c = (i-1)*ndof
        @gemm Bi = Bil*Rrot
        B[:, c+1:c+6] .= Bi

    end 

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
    t = elem.mat.t
    ndof = 6

    C = getcoords(elem)
    K = zeros(ndof*nnodes, ndof*nnodes)
    B = zeros(ndof, ndof*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    L = zeros(3,3)
    T = zeros(6,6)
    Rrot = Array{Float64}(I,5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)

    Dn = [ elem.Dlmn[i][j,3] for i in 1:nnodes, j in 1:3 ] # nx3

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        J = [ C'*dNdR + t/2*ip.R[3]*Dn'*dNdR   t/2*Dn'*N ] # 3x3

        J2D = C'*dNdR
        set_dir_matrix(elem, J2D, L)

        set_tensor_rot!(L, T)
        
        dNdR = [ dNdR zeros(nnodes) ]
        invJ = inv(J)
        dNdX = dNdR*invJ

        D = calcD(elem.mat, ip.state)
        # @showm J
        setB(elem, ip, invJ, N, dNdX, Rrot, Bil, Bi, B)
        detJ = det(J)

        coef = detJ*ip.w

        # @showm L
        # @showm T
        # @showm D
        # @showm B[:,1:6]
        # error()

        K += coef*B'*T'*D*T*B
    end

    δ = 1e-5
    # δ = 0.0
    for i in 1:ndof*nnodes
        # if abs(K[i,i]) < δ
            K[i,i] += δ
        # end
    end

    # @show diag(K)
    # error()

    map = elem_map(elem)
    return K, map, map
end

function elem_update!(elem::CookShell, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    t = elem.mat.t
    ndof = 6

    map = elem_map(elem)
    dU = U[map]
    dF = zeros(length(dU))

    C = getcoords(elem)
    B = zeros(ndof, ndof*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    L = zeros(3,3)
    T = zeros(6,6)
    Rrot = Array{Float64}(I,5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)
    Δε = zeros(6)

    Dn = [ elem.Dlmn[i][j,3] for i in 1:nnodes, j in 1:3 ] # nx3


    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J = [ C'*dNdR + t/2*ip.R[3]*Dn'*dNdR   t/2*Dn'*N ] # 3x3
        J2D = C'*dNdR
        set_dir_matrix(elem, J2D, L)
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
