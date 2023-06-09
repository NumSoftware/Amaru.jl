# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct TMShell<:Thermomechanical
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

    function TMShell();
        return new()
    end
end

matching_shape_family(::Type{TMShell}) = BULKCELL

function elem_config_dofs(elem::TMShell)
    nbnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)
        if  i<=(nbnodes)
            add_dof(node, :ut, :ft)
        end
    end
end

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


# DUVIDA!!!!!!
function distributed_bc(elem::TMShell, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_boundary_forces(elem, facet, key, val)
end


function body_c(elem::TMShell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, key, val)
end

function elem_config_dofs(elem::TMShell)
    ndim = elem.env.ndim
    ndim in (1,2) && error("MechShell: Shell elements do not work in $(ndim)d analyses")
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

function elem_map(elem::TMShell)
    keys =(:ux, :uy, :uz, :rx, :ry, :rz, :ut)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
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
    th = elem.mat.thickness
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

function elem_stiffness(elem::TMShell)
    nnodes = length(elem.nodes)
    th     = elem.mat.thickness
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

        D = calcD(elem.mat, ip.state)
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

    map = elem_map(elem)
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
    #th     = elem.env.thickness
    #nnodes = length(elem.nodes)
    #nbnodes = elem.shape.basic_shape.npoints
    #C   = getcoords(elem)
    #Bu  = zeros(6, nnodes*ndim)
    #Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    #J    = Array{Float64}(undef, ndim, ndim)
    #dNdX = Array{Float64}(undef, nnodes, ndim)
    #m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    #β    = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu) # thermal stress modulus

    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    ndim   = elem.env.ndim

    th     = elem.mat.thickness
    ndof   = 6
    nstr   = 6
    C      = getcoords(elem)
    K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rrot   = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)

    C   = getcoords(elem)
    Cut = zeros(ndof*nnodes, nbnodes) # u-t coupling matrix
    m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu) # ther

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

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
        Nt    = elem.shape.basic_shape.func(ip.R)
        coef  = β
        coef *= detJ′*ip.w*th
        mNt   = m*Nt'
        @gemm Cut -= coef*B'*mNt
    end
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cut, map_u, mat_t
end

# thermal conductivity
function elem_conductivity_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.mat.thickness   # elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    H      = zeros(nnodes, nnodes)
    dNtdX  = zeros(nnodes, ndim)
    Bt     = zeros(ndim, nnodes)
    KBt    = zeros(ndim, nnodes)
    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @gemm dNtdX = dNtdR*inv(J)
        Bt .= dNtdX'

        # compute H
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @gemm KBt = K*Bt
        @gemm H  -= coef*Bt'*KBt
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return H, map, map
end

function elem_mass_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = thickness   # elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)
    M      = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        Nt   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Cut
        coef  = elem.mat.ρ*elem.mat.cv
        coef *= detJ*ip.w*th
        M    -= coef*Nt*Nt'
    end

    # map
    map = [  node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes]  ]

    return M, map, map
end

#=
function elem_internal_forces(elem::TMShell, F::Array{Float64,1}, DU::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = thickness   # elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    T0     = elem.env.T0 + 273.15
    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]
    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nbnodes)
    Bt  = zeros(ndim, nbnodes)
    m = [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ] # = tI
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Jp  = Array{Float64}(undef, ndim, nbnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    dUt = DU[mat_t] # nodal temperature increments
    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)
        # compute Bu matrix and Bt
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        @gemm dNdX = dNdR*inv(J)
        set_Bu(elem, ip, dNdX, Bu)
        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        Jp = dNtdR*Ct
        @gemm dNtdX = inv(Jp)*dNtdR
        Bt = dNtdX
        # compute N
        # internal force
        ut   = ip.state.ut + 273
        β   = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu)
        σ    = ip.state.σ - β*ut*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ
        # internal volumes dFt
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = β*detJ*ip.w*th
        dFt  -= coef*Nt*εvol
        coef = detJ*ip.w*elem.mat.ρ*elem.mat.cv*th/T0
        dFt -= coef*Nt*ut
        QQ   = ip.state.QQ
        coef = detJ*ip.w*th/T0
        @gemv dFt += coef*Bt'*QQ
    end
    F[map_u] += dF
    F[mat_t] += dFt
end
=#


function elem_update!(elem::TMShell, DU::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = thickness   # elem.env.thickness
    T0     = elem.env.T0 + 273.15
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = getcoords(elem)

    E = elem.mat.E
    α = elem.mat.α
    ρ = elem.mat.ρ
    nu = elem.mat.nu
    cv = elem.mat.cv
    β = E*α/(1-2*nu)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUt = DU[mat_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1
    m   = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]  #

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, nnodes, ndim)
    dNtdX = Array{Float64}(undef, nbnodes, ndim)
    Δε = zeros(6)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @gemm dNdX = dNdR*invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm dNtdX = dNtdR*invJ

        # compute Nt
        Nt = elem.shape.basic_shape.func(ip.R)

        # compute Δε
        @gemv Δε = Bu*dU

        # compute Δut
        Δut = Nt'*dUt # interpolation to the integ. point

        # compute thermal gradient G
        Bt .= dNtdX'
        G  = Bt*Ut

        # internal force dF
        Δσ, q = stress_update(elem.mat, ip.state, Δε, Δut, G, Δt)
        Δσ -= β*Δut*m # get total stress

        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt

        coef  = ρ*cv
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt*Δut

        coef  = Δt
        coef *= detJ*ip.w*th
        @gemv dFt += coef*Bt'*q
    end

    DF[map_u] += dF
    DF[mat_t] += dFt
    return success()
end