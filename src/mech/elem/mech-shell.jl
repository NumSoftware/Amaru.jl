# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    MechShell
A shell finite element for mechanical equilibrium analyses.
"""
mutable struct MechShell<:Mechanical
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

    function MechShell();
        return new()
    end
end

matching_shape_family(::Type{MechShell}) = BULKCELL


function elem_init(elem::MechShell)

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

function setquadrature!(elem::MechShell, n::Int=0)

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


function distributed_bc(elem::MechShell, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_boundary_forces(elem, facet, key, val)
end

function body_c(elem::MechShell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, key, val)
end


function elem_config_dofs(elem::MechShell)
    ndim = elem.env.ndim
    ndim in (1,2) && error("MechShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end

function elem_map(elem::MechShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
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


function setB(elem::MechShell, ip::Ip, N::Vect, L::Matx, dNdX::Matx, Rrot::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.mat.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        ζ = ip.R[3]
        Rrot[1:3,1:3] .= L
        # Rrot[4:5,4:6] .= L[1:2,:]
        # Rrot[1:3,1:3] .= elem.Dlmn[i]
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


function elem_stiffness(elem::MechShell)
    nnodes = length(elem.nodes)
    th     = elem.mat.th
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
        J′    = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        D = calcD(elem.mat, ip.state)
        detJ′ = det(J′)
        @assert detJ′>0
        
        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)
        coef = detJ′*ip.w
        K += coef*B'*D*B
    end

    δ = 1e-7
    for i in 1:ndof*nnodes
        K[i,i] += δ
    end

    map = elem_map(elem)
    return K, map, map
end

function elem_update!(elem::MechShell, U::Array{Float64,1}, dt::Float64)
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
    L = zeros(3,3)
    Rrot = zeros(5,ndof)
    Bil = zeros(6,5)
    Bi = zeros(6,ndof)
    Δε = zeros(6)

    # Dn = [ elem.Dlmn[i][3,j] for i in 1:nnodes, j in 1:3 ] # nx3


    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn
        
        # J = [ C'*dNdR + th/2*ip.R[3]*Dn'*dNdR   th/2*Dn'*N ] # 3x3
        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        # J′ = L*J
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        setB(elem, ip, N, L, dNdX′, Rrot, Bil, Bi, B)
        Δε = B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return failure("MechShell: Error at integration point $(ip.id)")

        # detJ = det(J)
        detJ′ = det(J′)
        coef = detJ′*ip.w
        dF += coef*B'*Δσ

    end

     return dF, map, success()
end
