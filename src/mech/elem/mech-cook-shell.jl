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
    elem.shape==QUAD8 || error("elem_init: CookShell only works with shape QUAD8.")

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

    dir[:,1] = V1
    dir[:,2] = V2
    dir[:,3] = V3
end

# Rotation Matrix
function set_trans_matrix(elem::CookShell, dir::Matx, T::Matx)
    
    lx, ly, lz = dir[:,1]
    mx, my, mz = dir[:,2]
    nx, ny, nz = dir[:,3]

    T[1,1] = lx*lx;  T[1,2] = mx*mx;  T[1,3] = nx*nx;   T[1,4] =     2*mx*nx;  T[1,5] =     2*nx*lx;  T[1,6] =     2*lx*mx;
    T[2,1] = ly*ly;  T[2,2] = my*my;  T[2,3] = ny*ny;   T[2,4] =     2*my*ny;  T[2,5] =     2*ny*ly;  T[2,6] =     2*ly*my;
    T[3,1] = lz*lz;  T[3,2] = mz*mz;  T[3,3] = nz*nz;   T[3,4] =     2*mz*nz;  T[3,5] =     2*nz*lz;  T[3,6] =     2*lz*mz;
    T[4,1] = ly*lz;  T[4,2] = my*mz;  T[4,3] = ny*nz;   T[4,4] = my*nz+mz*ny;  T[4,5] = ly*nz+lz*ny;  T[4,6] = ly*mz+lz*my;
    T[5,1] = lz*lx;  T[5,2] = mz*mx;  T[5,3] = nz*nx;   T[5,4] = mz*nx+mx*nz;  T[5,5] = lz*nx+lx*nz;  T[5,6] = lz*mx+lx*mz;
    T[6,1] = lx*ly;  T[6,2] = mx*my;  T[6,3] = nx*ny;   T[6,4] = mx*ny+my*nx;  T[6,5] = lx*ny+ly*nx;  T[6,6] = lx*my+ly*mx;

    # @showm T
    # @showm R
    # error()
end

function setB(elem::CookShell, ip::Ip, invJ::Matx, N::Vect, dNdX::Matx, B::Matx)
    nnodes, ndim = size(dNdX)
    t = elem.mat.t
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 5

    for i in 1:nnodes
        c = (i-1)*ndof
        ζ = ip.R[3]

        lx, ly, lz = elem.Dlmn[i][:,1]
        mx, my, mz = elem.Dlmn[i][:,2]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        dNdz = dNdX[i,3]
        dζdx = invJ[1,3]
        dζdy = invJ[2,3]
        dζdz = invJ[3,3]
        dNζdx = dNdx*ζ+N[i]*dζdx
        dNζdy = dNdy*ζ+N[i]*dζdy
        dNζdz = dNdz*ζ+N[i]*dζdz

        B[1,1+c] = dNdx;                                                B[1,4+c] = -dNζdx*t/2*mx;               B[1,5+c] = dNζdx*t/2*lx
                             B[2,2+c] = dNdy;                           B[2,4+c] = -dNζdy*t/2*my;               B[2,5+c] = dNζdy*t/2*ly
                                                  B[3,3+c] = dNdz;      B[3,4+c] = -dNζdz*t/2*mz;               B[3,5+c] = dNζdz*t/2*lz
                             B[4,2+c] = dNdz/SR2; B[4,3+c] = dNdy/SR2;  B[4,4+c] = -1/SR2*(dNζdz*t/2*my+dNζdy*t/2*mz);  B[4,5+c] = 1/SR2*(dNζdz*t/2*ly+dNζdy*t/2*lz)
        B[5,1+c] = dNdz/SR2;                      B[5,3+c] = dNdx/SR2;  B[5,4+c] = -1/SR2*(dNζdz*t/2*mx+dNζdx*t/2*mz);  B[5,5+c] = 1/SR2*(dNζdz*t/2*lx+dNζdx*t/2*lz)
        B[6,1+c] = dNdy/SR2; B[6,2+c] = dNdx/SR2;                       B[6,4+c] = -1/SR2*(dNζdy*t/2*mx+dNζdx*t/2*my);  B[6,5+c] = 1/SR2*(dNζdy*t/2*lx+dNζdx*t/2*ly)

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
        #add_dof(node, :rz, :mz)
    end
end

function elem_map(elem::CookShell)
    keys =(:ux, :uy, :uz, :rx, :ry)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


function elem_stiffness(elem::CookShell)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    t = elem.mat.t

    C = getcoords(elem)
    K = zeros(5*nnodes, 5*nnodes)
    B = zeros(6, 5*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    L = zeros(3,3)
    T = zeros(6,6)

    Dn = [ elem.Dlmn[i][j,3] for i in 1:nnodes, j in 1:3 ] # nx3

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        J = [ C'*dNdR + t/2*ip.R[3]*Dn'*dNdR   t/2*Dn'*N ] # 3x3

        J2D = C'*dNdR
        set_dir_matrix(elem, J2D, L)
        set_trans_matrix(elem, L, T)

        dNdR = [ dNdR zeros(nnodes) ]
        invJ = inv(J)
        dNdX = dNdR*invJ

        D = calcD(elem.mat, ip.state)
        setB(elem, ip, invJ, N, dNdX, B)
        detJ = det(J)

        coef = detJ*ip.w
        
        K += coef*B'*T'*D*T*B
    end

    map = elem_map(elem)
    return K, map, map
end

function elem_update!(elem::CookShell, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    t = elem.mat.t

    map = elem_map(elem)
    dU = U[map]
    dF = zeros(length(dU))

    C = getcoords(elem)
    B = zeros(6, 5*nnodes)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Δε = zeros(6)
    L = zeros(3,3)
    T = zeros(6,6)

    Dn = [ elem.Dlmn[i][j,3] for i in 1:nnodes, j in 1:3 ] # nx3


    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        J = [ C'*dNdR + t/2*ip.R[3]*Dn'*dNdR   t/2*Dn'*N ] # 3x3
        J2D = C'*dNdR
        set_dir_matrix(elem, J2D, L)
        set_trans_matrix(elem, L, T)

        dNdR = [ dNdR zeros(nnodes) ]
        invJ = inv(J)
        dNdX = dNdR*invJ
        setB(elem, ip, invJ, N, dNdX, B)
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
