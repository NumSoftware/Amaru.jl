# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    ShellDegenerated
A bulk finite element for mechanical equilibrium analyses.
"""
mutable struct ShellDegenerated<:Mech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    active::Bool
    linked_elems::Array{Element,1}
    ctx::Context
    Dlmn::Array{ Array{Float64,2}, 1}

    function ShellDegenerated();
        return new()
    end
end

compat_shape_family(::Type{ShellDegenerated}) = BULKCELL


function elem_init(elem::ShellDegenerated)
    elem.shape==QUAD8 || error("elem_init: ShellDegenerated only works with shape QUAD8.")

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
        #V2 = cross(V3, V1)
        V2 = cross(V1, V3)
        normalize!(V1)
        normalize!(V2)
        normalize!(V3)

        push!(Dlmn, [V1 V2 V3 ])
    end
    elem.Dlmn = Dlmn
   
    return nothing
end


function setquadrature!(elem::ShellDegenerated, n::Int=0)

    # if !(n in keys(elem.shape.quadrature))
    #     alert("setquadrature!: cannot set $n integration points for shape $(elem.shape.name)")
    #     return
    # end

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
            elem.ips[j].state = compat_state_type(elem.mat)(elem.ctx)
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


function distributed_bc(elem::ShellDegenerated, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.ctx.ndim
    th    = elem.ctx.t
    suitable_keys = (:tx, :ty, :tz, :tn)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet !== nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.ctx.t

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

    for i in 1:size(ips,1)
        R = ips[i].coord
        w = ips[i].w
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
        X = C'*N
        if ndim==2
            x, y = X
            vip = evaluate(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
            if elem.ctx.stressmodel==:axisymmetric
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = evaluate(val, t=t, x=x, y=y, z=z)
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
        @mul F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end

# Rotation Matrix
function set_dir_matrix(elem::ShellDegenerated, J::Matx, dir::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    #V2 = cross(V3, V1)
    V2 = cross(V1, V3)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    dir[:,1] = V1
    dir[:,2] = V2
    dir[:,3] = V3
end


function set_jacobian_matrix(elem::ShellDegenerated, C::Matx, R::Matx, J::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    #V2 = cross(V3, V1)
    V2 = cross(V1, V3)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    dir[:,1] = V1
    dir[:,2] = V2
    dir[:,3] = V3
end

# Rotation Matrix
function set_trans_matrix(elem::ShellDegenerated, dir::Matx, T::Matx)


    lx, ly, lz = dir[:,1]
    mx, my, mz = dir[:,2]
    nx, ny, nz = dir[:,3]

    T[1,1] =     lx*lx;  T[1,2] =     ly*ly;  T[1,3] =     lz*lz;    T[1,4] =       lx*ly; T[1,5] =       lz*lx;   T[1,6] =       ly*lz; 
    T[2,1] =     mx*mx;  T[2,2] =     my*my;  T[2,3] =     mz*mz;    T[2,4] =       mx*my; T[2,5] =       mz*mx;   T[2,6] =       my*mz; 
   #T[3,1] =     nx*nx;  T[3,2] =     ny*ny;  T[3,3] =     nz*nz;    T[3,4] =       nx*ny; T[3,5] =       nz*nx;   T[3,6] =       ny*nz; 
    T[3,1] =   2*lx*mx;  T[3,2] =   2*ly*my;  T[3,3] =   2*lz*mz;    T[3,4] = lx*my+mx*ly; T[3,5] = lx*mz+mx*lz;   T[3,6] = ly*mz+my*lz; 
    T[4,1] =   2*nx*lx;  T[4,2] =   2*ny*ly;  T[4,3] =   2*nz*lz;    T[4,4] = nx*ly+lx*ny; T[4,5] = nx*lz+lx*nz;   T[4,6] = ny*lz+ly*nz; 
    T[5,1] =   2*mx*nx;  T[5,2] =   2*my*ny;  T[5,3] =   2*mz*nz;    T[5,4] = mx*ny+nx*my; T[5,5] = mx*nz+nx*mz;   T[5,6] = my*nz+ny*mz;


    #indexL = [1,2,3,5,4]
    #indexC = [1,2,3,4,6,5]
    #T = T[indexL,indexC]
#=
    
    l1, m1, n1 = dir[:,1]
    l2, m2, n2 = dir[:,2]
    l3, m3, n3 = dir[:,3]
    
    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;   T[1,4] =       l1*m1;  T[1,5] =       m1*n1;  T[1,6] =       n1*l1;
    T[2,1] =     l2*l2;  T[2,2] =     m2*m2;  T[2,3] =     n2*n2;   T[2,4] =       l2*m2;  T[2,5] =       m2*n2;  T[2,6] =       n2*l2;
    T[3,1] =   2*l1*l2;  T[3,2] =   2*m1*m2;  T[3,3] =   2*n1*n2;   T[3,4] = l1*m2+l2*m1;  T[3,5] = m1*n2+m2*n1;  T[3,6] = n1*l2+n2*l1;
    T[4,1] =   2*l2*l3;  T[4,2] =   2*m2*m3;  T[4,3] =   2*n2*n3;   T[4,4] = l2*m3+l3*m2;  T[4,5] = m2*n3+m3*n2;  T[4,6] = n2*l3+n3*l2;
    T[5,1] =   2*l3*l1;  T[5,2] =   2*m3*m1;  T[5,3] =   2*n3*n1;   T[5,4] = l3*m1+l1*m3;  T[5,5] = m3*n1+m1*n3;  T[5,6] = n3*l1+n1*l3;
   =# 
     #@showm T
     #error()
end


function setB(elem::ShellDegenerated, R::Matx, J::Matx , ip::Ip, dNdR::Matx, dNdX::Matx, N::Vect, B::Matx)
    nnodes, ndim = size(dNdX)
    ndof = 5
    B .= 0.0
    t = elem.mat.t
    ζ = ip.R[3]
    
       for i in 1:nnodes


        #=
        lx, ly, lz = elem.Dlmn[i][:,2]
        mx, my, mz = elem.Dlmn[i][:,1]
        =#
           
        lx, ly, lz = elem.Dlmn[i][:,1]
        mx, my, mz = elem.Dlmn[i][:,2]

        #dNdx = dNdX[i,1]
        #dNdy = dNdX[i,2]

        dNdxi  = dNdR[i,1]
        dNdeta = dNdR[i,2]

          J_inv = inv(J)'
          #@show J_inv
          #error()
          H1 = J_inv[1,1]*dNdxi + J_inv[1,2]*dNdeta
          H2 = J_inv[2,1]*dNdxi + J_inv[2,2]*dNdeta
          H3 = J_inv[3,1]*dNdxi + J_inv[3,2]*dNdeta

          G1 = (J_inv[1,1]*dNdxi  + J_inv[1,2]*dNdeta)*ζ +  J_inv[1,3]*N[i]
          G2 = (J_inv[2,1]*dNdxi  + J_inv[2,2]*dNdeta)*ζ +  J_inv[2,3]*N[i]
          G3 = (J_inv[3,1]*dNdxi  + J_inv[3,2]*dNdeta)*ζ +  J_inv[3,3]*N[i]
          
          #=
          g11 = -t/2*mx
          g12 = -t/2*my
          g13 = -t/2*mz

          g21 = -t/2*lx
          g22 = -t/2*ly
          g23 = -t/2*lz
          =#
          
          g11 = -t/2*lx
          g12 = -t/2*ly
          g13 = -t/2*lz

          g21 = -t/2*mx
          g22 = -t/2*my
          g23 = -t/2*mz
    

          j    = i-1

          B[1,1+j*ndof] = H1;                                                B[1,4+j*ndof] = g11*G1;         B[1,5+j*ndof] = g21*G1

                               B[2,2+j*ndof] = H2;                           B[2,4+j*ndof] = g12*G2;         B[2,5+j*ndof] = g22*G2

                                                     B[3,3+j*ndof] = H3;     B[3,4+j*ndof] = g13*G3;         B[3,5+j*ndof] = g23*G3

          B[4,1+j*ndof] = H2;  B[4,2+j*ndof] = H1;                           B[4,4+j*ndof] = g11*G2+g12*G1;  B[4,5+j*ndof] = g21*G2+g22*G1

          B[5,1+j*ndof] = H3;                        B[5,3+j*ndof] = H1;     B[5,4+j*ndof] = g11*G3+g13*G1;  B[5,5+j*ndof] = g21*G3+g23*G1

                               B[6,2+j*ndof] = H3;   B[6,3+j*ndof] = H2;     B[6,4+j*ndof] = g12*G3+g13*G2;  B[6,5+j*ndof] = g22*G3+g23*G2
    
    #=
          B[1,1+j*ndof] = H1;                                                B[1,4+j*ndof] = g11*G1;         B[1,5+j*ndof] = g21*G1
                               B[2,2+j*ndof] = H2;                           B[2,4+j*ndof] = g12*G2;         B[2,5+j*ndof] = g22*G2
                                                     B[3,3+j*ndof] = H3;     B[3,4+j*ndof] = g13*G3;         B[3,5+j*ndof] = g23*G3
                               B[4,2+j*ndof] = H3;   B[4,3+j*ndof] = H2;     B[4,4+j*ndof] = g12*G3+g13*G2;  B[4,5+j*ndof] = g22*G3+g23*G2
          B[5,1+j*ndof] = H3;                        B[5,3+j*ndof] = H1;     B[5,4+j*ndof] = g11*G3+g13*G1;  B[5,5+j*ndof] = g21*G3+g23*G1
                                  
          B[6,1+j*ndof] = H2;  B[6,2+j*ndof] = H1;                           B[6,4+j*ndof] = g11*G2+g12*G1;  B[6,5+j*ndof] = g21*G2+g22*G1
     =#
           end 
    #@showm B
    #error()
end


function setD(elem::ShellDegenerated, D::Matx)

    nu = elem.mat.ν
    E1 = elem.mat.E/(1-elem.mat.ν^2)
    G  = elem.mat.E/(2*(1+elem.mat.ν))
    G1 = 5/6*G

              D .=   [E1      nu*E1  0  0   0
                      nu*E1   E1     0  0   0
                      0       0      G  0   0
                      0       0      0  G1  0
                      0       0      0  0  G1 ]

end


function elem_config_dofs(elem::ShellDegenerated)
    ndim = elem.ctx.ndim
    ndim in (1,2) && error("ShellDegenerated: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        #add_dof(node, :rz, :mz)
    end
end


function elem_map(elem::ShellDegenerated)::Array{Int,1}

    #dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)

    keys =(:ux, :uy, :uz, :rx, :ry)
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

end


function elem_stiffness(elem::ShellDegenerated)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    K = zeros(5*nnodes, 5*nnodes)
    B = zeros(6, 5*nnodes)
    J  = Array{Float64}(undef, ndim, ndim)
    R = zeros(3,3)
    L = zeros(3,3)
    T = zeros(5,6)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    
    D  = Array{Float64}(undef, 5, 5)
    setD(elem, D)

    t = elem.mat.t
    C = getcoords(elem)

    Dn = [ elem.Dlmn[i][j,3] for i in 1:nnodes, j in 1:3 ] # nx3

    for ip in elem.ips
        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        N    = elem.shape.func(ip.R)
        J = [ C'*dNdR + t/2*ip.R[3]*Dn'*dNdR   t/2*Dn'*N ] # 3x3
        J2D = C'*dNdR

        set_dir_matrix(elem, J2D, L)
        set_trans_matrix(elem, L, T)

        #indexL = [1,2,3,5,4]
        #indexC = [1,2,3,4,6,5]
        #T = T[indexL,indexC]
        #@showm T
        #error()

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX = dNdR*inv(J)
        #@showm  dNdX
        #@show ζ
        #error()
      

        setB(elem, R, J, ip, dNdR, dNdX, N, B)  # 6x40

          #=
        @showm L
        @showm inv(J)
        @showm B[1:6,1:5]
        error()  
        =#

        #Bs = B[1:3,:]
        #Bt = B[4:5,:]

        #Ds = D[1:3,1:3]
        #Dt = D[4:5,4:5]

        #Ts = T[1:3,:]
        #Tt = T[4:5,:]

        #Bs = Ts*B
        #Bt = Tt*B

        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

  
        #coef1 = detJ*ip.w  # check *0.5
        #coef2 = detJ*ip.w   # check *0.5

        #Ks = Bs'*Ds*Bs*coef1
        #Kt = Bt'*Dt*Bt*coef2
        #K += Ks + Kt
    
        coef = detJ*ip.w 
        K += coef*B'*T'*D*T*B
     

    end
    
     map = elem_map(elem)
    return K, map, map
    
end


function update_elem!(elem::ShellDegenerated, U::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end


#=
function update_elem!(elem::ShellDegenerated, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    dU = U[map]
    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)
    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Δε = zeros(6)
    C = getcoords(elem)
    for ip in elem.ips
        if elem.ctx.stressmodel==:axisymmetric
            th = 2*pi*ip.coord.x
        end
        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || er
        =#