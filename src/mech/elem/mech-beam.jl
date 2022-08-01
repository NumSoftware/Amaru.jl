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

matching_shape_family(::Type{MechBeam}) = LINE_CELL


function elem_init(elem::MechBeam)

    nnodes = length(elem.nodes)
    Dlmn = Array{Float64,2}[]
    C = getcoords(elem)

    for i in 1:nnodes
        Ri = elem.shape.nat_coords[i,:]
        dNdR = elem.shape.deriv(Ri)
        J = C'*dNdR
        R = zeros(3,3)
        set_rot_x_xp(elem, J, R)
        push!(Dlmn, R)
    end
    elem.Dlmn = Dlmn

    return nothing
end

function setquadrature!(elem::MechBeam, n::Int=0)

    ipL = get_ip_coords(LIN2, n) # longitudinal
    ipT = get_ip_coords(LIN2, 2) # transversal
    
    n = size(ipL,1)

    resize!(elem.ips, 4*n)
    for i in 1:n
        for j in 1:2
            for k in 1:2
                R = [ ipT[i,1], ipT[j,1], ipL[k,1] ]
                w = ipL[i,4]*ipT[j,4]*ipT[k,4]
                m = (i-1)*n*2 + (j-1)*2 + k
                elem.ips[m] = Ip(R, w)
                elem.ips[m].id = m
                elem.ips[m].state = ip_state_type(elem.mat)(elem.env)
                elem.ips[m].owner = elem
                # @show m
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
    return mech_line_distributed_forces(elem, key, val)
end


function body_c(elem::MechBeam, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
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
function set_rot_x_xp(elem::MechBeam, J::Matx, R::Matx)
    V1 = normalize(vec(J))

    if V1[1]==0.0
        V2 = [ 1.0, 0.0, 0.0 ]
    elseif V1[2]==0.0
        V2 = [ 0.0, 1.0, 0.0 ]
    else
        V2 = [ 0.0, 0.0, 1.0 ]
    end

    V3 = cross(V1, V2)
    V2 = cross(V3, V1)

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
        dNdx = dNdX[i]

        Bil[1,1] = dNdx;                                                                             Bil[1,5] = dNdx*ζ*th/2;  Bil[1,6] = -dNdx*η*ts/2
                                               Bil[2,3] = dNdx/SR2;  Bil[2,4] = 1/SR2*dNdx*η*ts/2;   Bil[2,5] = 1/SR2*Ni
                          Bil[3,2] = dNdx/SR2;                       Bil[3,4] = -1/SR2*dNdx*ζ*th/2;                           Bil[3,6] = -1/SR2*Ni

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
    Bil = zeros(nstr,6)
    Bi = zeros(nstr,ndof)

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        D = calcD(elem.mat, ip.state)
        
        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)

        detJ′ = dx′dξ*th/2*ts/2
        coef = detJ′*ip.w
        K += coef*B'*D*B

    end

    map = elem_map(elem)
    return K, map, map
end


function elem_update!(elem::MechBeam, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    nnodes = length(elem.nodes)
    th = elem.mat.th
    ts = elem.mat.ts
    ndof = 6
    nstr = 3
    C = getcoords(elem)
    B = zeros(nstr, ndof*nnodes)
    L = zeros(3,3)
    Rθ = zeros(6,ndof)
    Bil = zeros(nstr,6)
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


function elem_vals(elem::MechBeam)
    # get ip average values
    ipvals = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
    sum  = merge(+, ipvals... )
    nips = length(elem.ips)
    vals = OrderedDict( k=>v/nips for (k,v) in sum)
    return vals
end


function elem_extrapolated_node_vals(elem::MechBeam)
    nnodes = length(elem.nodes)

    # displacements
    keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    U = [ node.dofdict[key].vals[key] for node in elem.nodes for key in keys ]
    Rθ = zeros(6,6)
    Uplane = zeros(2*nnodes, 2) # u1 θ1 u2 θ2; u1 θ1 ...
    for i in 1:nnodes
        Rθ[1:3,1:3] .= elem.Dlmn[i]
        Rθ[4:6,4:6] .= elem.Dlmn[i]
        U′i = Rθ*U[(i-1)*6+1:i*6]
        Uplane[(i-1)*2+1, 1] = U′i[3] # uz
        Uplane[(i-1)*2+2, 1] = U′i[5] # θy
        Uplane[(i-1)*2+1, 2] = U′i[2] # uy
        Uplane[(i-1)*2+2, 2] = U′i[6] # θz
    end

    # Coefficients matrix for Hermite interpolation
    ncoefs = nnodes*2
    A = zeros(ncoefs, ncoefs)
    Ξ = elem.shape.nat_coords
    for i in 1:nnodes
        ξ = Ξ[i]
        A[(i-1)*2+1,1] = 1.0
        for j in 2:ncoefs
            A[(i-1)*2+1,j] = ξ^(j-1)
            A[(i-1)*2+2,j] = (j-1)*ξ^(j-2)
        end
    end

    invA = inv(A)
    M = zeros(nnodes,2)
    V = zeros(nnodes,2)
    E = elem.mat.E
    th = elem.mat.th
    ts = elem.mat.ts

    for k in 1:2
        C = invA*Uplane[:,k]
        
        if k==1
            I = ts*th^3/12 #!TODO, check directions
        else
            I = th*ts^3/12
        end

        for i in 1:nnodes
            ξ = Ξ[i]
            d2u = 2*C[3]
            d3u = 0.0
            for j in 4:ncoefs
                d2u += factorial(j-1)/factorial(j-3)*C[j]*ξ^(i-3)
                d3u += factorial(j-1)/factorial(j-4)*C[j]*ξ^(i-4)
            end
            M[i,k] = d2u/(E*I)
            V[i,k] = d3u/(E*I)
        end
    end

    return OrderedDict(
        Symbol("Mx'z'") => M[:,1],
        Symbol("Vx'z'") => V[:,1],
        Symbol("Mx'y'") => M[:,2],
        Symbol("Vx'y'") => V[:,2]
    )

end

