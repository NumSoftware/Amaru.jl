# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    MechBeam
A beam finite element for mechanical equilibrium analyses.
"""
mutable struct MechBeam<:Mechanical
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

    function MechBeam();
        return new()
    end
end

matching_shape_family(::Type{MechBeam}) = SOLID_CELL


function elem_init(elem::MechBeam)

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

function setquadrature!(elem::MechBeam, n::Int=0)

    ipL = get_ip_coords(LIN2, n)
    ipT = get_ip_coords(LIN2, 2)

    resize!(elem.ips, 4*n)
    for i=1:n
        for j in 1:2
            for k in 1:2
                R = [ ipT[i,1], ipT[j,1], ipL[k,1] ]
                w = ipL[i,4]*ipT[j,4]*ipT[k,4]
                m = (i-1)*n + (j-1)*2 + k
                elem.ips[m] = Ip(R, w)
                elem.ips[m].id = m
                elem.ips[m].state = ip_state_type(elem.mat)(elem.env)
                elem.ips[m].owner = elem
            end
        end
    end

    # finding ips global coordinates
    C     = getcoords(elem)
    shape = elem.shape

    for ip in elem.ips
        R = [ ip.R[1], 0.0, 0.0 ]
        N = shape.func(R)
        ip.coord = C'*N
    end

end


function distributed_bc(elem::MechBeam, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.mat.th
    suitable_keys = (:tx, :ty, :tz, :tn)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")


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


function elem_config_dofs(elem::MechBeam)
    ndim = elem.env.ndim
    ndim in (1,2) && error("MechBeam: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end

function elem_map(elem::MechBeam)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


# Rotation Matrix

function mount_T(elem::MechBeam, J::Matx, R::Matx)
    V1 = normalize(J)

    if V1[1]==0.0
        V3 = [ 1.0, 0.0, 0.0 ]
    elseif V[2]==0.0
        V3 = [ 0.0, 1.0, 0.0 ]
    else
        V3 = [ 0.0, 0.0, 1.0 ]
    end

    V2 = cross(V3, V1)
    V3 = cross(V1, V2)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    R[1,:] .= V1
    R[2,:] .= V2
    R[3,:] .= V3
end


function setB(elem::MechBeam, ip::Ip, L::Matx, N::Vect, dNdX::Matx, Rθ::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.mat.th
    ts = elem.mat.ts
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        η = ip.R[2]
        ζ = ip.R[3]
        Rθ[1:3,1:3] .= L
        Rθ[4:6,4:6] .= elem.Dlmn[i]
        Ni = N[i]
        dNdx = dNdR[i]/normJ

        Bil[1,1] = dNdx;                                                                             Bil[1,5] = dNdx*ζ*th/2;  Bil[1,6] = dNdx*η*ts/2
                                               Bil[2,3] = dNdx/SR2;  Bil[2,4] = 1/SR2*dNdz*η*ts/2;   Bil[2,5] = 1/SR2*Ni
                          Bil[3,2] = dNdx/SR2;                       Bil[3,4] = -1/SR2*dNdz*ζ*th/2;                           Bil[3,6] = -1/SR2*Ni

        c = (i-1)*ndof
        @gemm Bi = Bil*Rθ
        B[:, c+1:c+6] .= Bi
    end
end


function elem_stiffness(elem::MechBeam)
    nnodes = length(elem.nodes)
    th = elem.mat.th
    ts = elem.mat.ts
    ndof = 6
    nstr = 3
    C = getcoords(elem)
    K = zeros(ndof*nnodes, ndof*nnodes)
    dNdX′ = zeros(nnodes)
    B = zeros(nstr, ndof*nnodes)
    L = zeros(3,3)
    Rθ = zeros(6,ndof)
    Bil = zeros(nstr,5)
    Bi = zeros(nstr,ndof)

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        normJ = dx′dξ*th/2*ts/2
        D = calcD(elem.mat, ip.state)
        
        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)
        coef = normJ*ip.w
        K += coef*B'*D*B
    end

    map = elem_map(elem)
    return K, map, map
end

function elem_update!(elem::MechBeam, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    nnodes = length(elem.nodes)
    th = elem.mat.th
    ndof = 6
    nstr = 3
    C = getcoords(elem)
    B = zeros(nstr, ndof*nnodes)
    L = zeros(3,3)
    Rθ = zeros(6,ndof)
    Bil = zeros(nstr,5)
    Bi = zeros(nstr,ndof)

    map = elem_map(elem)
    dU = U[map]
    dF = zeros(length(dU))
    Δε = zeros(6)

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        normJ = dx′dξ*th/2*ts/2
        
        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)
        Δε = B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return failure("MechBeam: Error at integration point $(ip.id)")
        
        coef = normJ*ip.w
        dF += coef*B'*Δσ
    end

    F[map] += dF
    return success()
end
