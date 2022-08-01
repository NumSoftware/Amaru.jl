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
    env::ModelEnv

    function TMShell();
        return new()
    end
end

matching_shape_family(::Type{TMShell}) = SOLID_CELL

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
    nothing
end


function distributed_bc(elem::TMShell, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th     = thickness   # elem.env.thickness
    suitable_keys = (:tx, :ty, :tz, :tn, :tq)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet!=nothing ? facet : elem
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

    if key == :tq # energy per area
        for i in 1:size(ips,1)
            R = vec(ips[i,:])
            w = R[end]
            N = shape.func(R)
            D = shape.deriv(R)
            J = C'*D
            nJ = norm2(J)
            X = C'*N
            if ndim==2
                x, y = X
                vip = eval_arith_expr(val, t=t, x=x, y=y)
            else
                x, y, z = X
                vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            end
            coef = vip*nJ*w
            F .+= N*coef # F is a vector
        end

        # generate a map
        map  = [ node.dofdict[:ut].eq_id for node in target.nodes ]

        return F, map
    end

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
            if elem.env.modeltype=="axisymmetric"
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            if key == :tx
                Q = [vip, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, vip, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, vip]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = vip*normalize(n)
            end
        end
        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end

# the strain-displacement matrix for membrane forces
function Dm_maxtrix(elem::TMShell)

    coef1 = elem.mat.t*elem.mat.E/(1-elem.mat.nu^2)
    coef2 = elem.mat.nu*coef1
    coef3 = coef1*(1-elem.mat.nu)/2

        Dm = [coef1  coef2 0
                  coef2  coef1 0
                  0      0     coef3]
    return Dm
end

# the strain-displacement matrix for bending moments
function Db_maxtrix(elem::TMShell)

    Dm = Dm_maxtrix(elem)

    Db = Dm*(elem.mat.t^2/12)

    return Db
end

# the strain-displacement matrix for shear forces

function Ds_maxtrix(elem::TMShell)

    coef = elem.mat.t*(5/6)*elem.mat.E/(2*(1+elem.mat.nu))

            Ds = [coef    0
                        0     coef]
    return Ds
end

# Rotation Matrix
function RotMatrix(elem::TMShell, J::Matrix{Float64})
    
    Z = zeros(1,2) # zeros(2,1)

    if size(J,1)==2
        J = [J
             Z]
    else
        J = J
    end
    
    L1 = vec(J[:,1])
    L2 = vec(J[:,2])
    L3 = cross(L1, L2)  # L1 is normal to the first element face
    L2 = cross(L1, L3)
    normalize!(L1)
    normalize!(L2)
    normalize!(L3)

    Z1 = zeros(1,2) # Z = zeros(1,3)

    Rot = [ L2' Z1
    L1' Z1
    L3' Z1
    Z1   L2'
    Z1   L1']

    return Rot
             
end

function setBb(elem::TMShell, N::Vect, dNdX::Matx, Bb::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)

    ndof = 5
    Bb .= 0.0
   
    for i in 1:nnodes
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        j    = i-1

        Bb[1,4+j*ndof] = -dNdx  
        Bb[2,5+j*ndof] = -dNdy   
        Bb[3,4+j*ndof] = -dNdy 
        Bb[3,5+j*ndof] = -dNdx 

    end
end

function setBm(elem::TMShell, N::Vect, dNdX::Matx, Bm::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)
    ndof = 5
    Bm .= 0.0
   
    for i in 1:nnodes
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        j    = i-1

        Bm[1,1+j*ndof] = dNdx  
        Bm[2,2+j*ndof] = dNdy   
        Bm[3,1+j*ndof] = dNdy 
        Bm[3,2+j*ndof] = dNdx 

    end
end

function setBs_bar(elem::TMShell, N::Vect, dNdX::Matx, Bs_bar::Matx)
    nnodes = length(elem.nodes)

    cx = [ 0 1 0 -1]
    cy = [-1 0 1 0 ]

    Bs_bar .= 0.0
    Ns= zeros(4,1)
    
    for i in 1:nnodes
      Ns[1] = (1-cx[i])*(1-cy[i])/4
      Ns[2] = (1+cx[i])*(1-cy[i])/4
      Ns[3] = (1+cx[i])*(1+cy[i])/4
      Ns[4] = (1-cx[i])*(1+cy[i])/4

      bs1  = [ dNdX[1,1] -Ns[1]    0
               dNdX[1,2]     0 -Ns[1]];

      bs2  = [ dNdX[2,1] -Ns[2]    0
               dNdX[2,2]     0 -Ns[2]]

      bs3  = [ dNdX[3,1] -Ns[3]    0
               dNdX[3,2]     0 -Ns[3]]

      bs4  = [ dNdX[4,1] -Ns[4]    0
               dNdX[4,2]     0 -Ns[4]]

          bs = [bs1 bs2 bs3 bs4]

          Bs_bar[2*i-1:2*i,:] = bs[1:2,:]
    end
end

#=
function elem_map(elem::TMShell)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:ux, :uy, :uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end
=#

function elem_stiffness(elem::TMShell)
    ndim   = elem.env.ndim

    nnodes = length(elem.nodes)

    Db = Db_maxtrix(elem)
    Dm = Dm_maxtrix(elem)
    Ds = Ds_maxtrix(elem)

    Bb = zeros(3, nnodes*5)
    Bm = zeros(3, nnodes*5)
    Bs_bar = zeros(8,nnodes*3)

    c  = zeros(8,8)
    nr = 5   
    nc = 5
    R = zeros(nnodes*nr, nnodes*nc)
    K = zeros( nnodes*5 , nnodes*5 )

    C = getcoords(elem)

    if size(C,2)==2
        cxyz  = zeros(4,3)
        cxyz[:,1:2]  = C
    else
        cxyz = C
    end

    for ip in elem.ips      
        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        J = cxyz'*dNdR

        Ri = RotMatrix(elem, J)
        Ri′ = Ri[1:2, 1:3]
    
        ctxy = cxyz*Ri[1:3, 1:3]' # Rotate coordinates to element mid plane
      
        dNdX = dNdR*pinv(J)

        dNdX′ = dNdX*(Ri′)'
              
        for i in 1:nnodes
            R[(i-1)*nr+1:i*nr, (i-1)*nc+1:i*nc] = Ri
        end     
   
        J1 = ctxy'*dNdR
        invJ1  =  pinv(J1)
        detJ1 = norm2(J1)

        for i in 1:nnodes
            c[(i-1)*2+1:i*2, (i-1)*2+1:i*2] = J1[1:2,1:2]
        end

        setBb(elem, N, dNdX′, Bb)
        setBm(elem, N, dNdX′, Bm)
        setBs_bar(elem, N, dNdX′, Bs_bar)

                    T_mat = [ 1  0  0  0  0  0  0  0
                              0  0  0  1  0  0  0  0
                              0  0  0  0  1  0  0  0
                              0  0  0  0  0  0  0  1 ]

                    P_mat = [ 1  -1   0   0
                              0   0   1   1
                              1   1   0   0
                              0   0   1  -1 ]
       
                    A_mat = [ 1  ip.R[2] 0    0
                              0    0  1  ip.R[1]]
      
                    bmat_ss = invJ1[1:2, 1:2]* A_mat * inv(P_mat) * T_mat * c * Bs_bar

                    bmat_s1 = [0  0 bmat_ss[1, 1]
                               0  0 bmat_ss[2, 1]]

                    bmat_s2 = [0  0 bmat_ss[1, 4]
                               0  0 bmat_ss[2, 4]]

                    bmat_s3 = [0  0 bmat_ss[1, 7]
                               0  0 bmat_ss[2, 7]]

                    bmat_s4 = [0  0 bmat_ss[1,10]
                               0  0 bmat_ss[2,10]]

                    bmat_s1 = [bmat_s1*Ri[1:3, 1:3]  bmat_ss[:,2:3]]
                    bmat_s2 = [bmat_s2*Ri[1:3, 1:3] bmat_ss[:,5:6]]
                    bmat_s3 = [bmat_s3*Ri[1:3, 1:3] bmat_ss[:,8:9]]
                    bmat_s4 = [bmat_s4*Ri[1:3, 1:3] bmat_ss[:,11:12]]

                    bmat_s = [bmat_s1 bmat_s2 bmat_s3 bmat_s4]

                    coef = detJ1*ip.w

                     Kb =    Bb'*Db*Bb*coef 
                     Km = R'*Bm'*Dm*Bm*R*coef 
                     Ks = bmat_s'*Ds*bmat_s*coef 

                     K += (Kb + Km + Ks)
            end
            
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    #map = elem_map(elem) 

    return K, map, map
end


@inline function set_Bu1(elem::Element, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end


# matrix C
function elem_coupling_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = t   # elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = getcoords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu) # thermal stress modulus

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        set_Bu1(elem, ip, dNdX, Bu)

        # compute Cut
        Nt    = elem.shape.basic_shape.func(ip.R)
        coef  = β
        coef *= detJ*ip.w*th
        mNt   = m*Nt'
        @gemm Cut -= coef*Bu'*mNt
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
    th     = thickness   # elem.env.thickness
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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
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


function elem_update!(elem::TMShell, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
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
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
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