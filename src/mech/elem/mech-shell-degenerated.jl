# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    ShellDegenerated
A bulk finite element for mechanical equilibrium analyses.
"""
mutable struct ShellDegenerated<:Mechanical
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv

    function ShellDegenerated();
        return new()
    end
end

matching_shape_family(::Type{ShellDegenerated}) = SOLID_CELL


function elem_init(elem::ShellDegenerated)
    elem.shape==QUAD8 || error("elem_init: ShellDegenerated only works with shape QUAD8.")
    
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


function distributed_bc(elem::ShellDegenerated, facet::Union{Facet, Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.t
    suitable_keys = (:tx, :ty, :tz, :tn)

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

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
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

# Rotation Matrix
function rot_matrix_R(elem::ShellDegenerated, J::Matx, R::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    V2 = cross(V3, V1)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    R[:,1] = V1
    R[:,2] = V2
    R[:,3] = V3
end

# Rotation Matrix
function rot_matrix_T(elem::ShellDegenerated, R::Matx, T::Matx)
    
    l1, m1, n1 = R[:,1]
    l2, m2, n2 = R[:,2]
    l3, m3, n3 = R[:,3]
    
    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;   T[1,4] =       l1*m1;  T[1,5] =       m1*n1;  T[1,6] =       n1*l1;
    T[2,1] =     l2*l2;  T[2,2] =     m2*m2;  T[2,3] =     n2*n2;   T[2,4] =       l2*m2;  T[2,5] =       m2*n2;  T[2,6] =       n2*l2;
    T[3,1] =     l3*l3;  T[3,2] =     m3*m3;  T[3,3] =     n3*n3;   T[3,4] =       l3*m3;  T[3,5] =       m3*n3;  T[3,6] =       n3*l3;
    T[4,1] =   2*l1*l2;  T[4,2] =   2*m1*m2;  T[4,3] =   2*n1*n2;   T[4,4] = l1*m2+l2*m1;  T[4,5] = m1*n2+m2*n1;  T[4,6] = n1*l2+n2*l1;
    T[5,1] =   2*l2*l3;  T[5,2] =   2*m2*m3;  T[5,3] =   2*n2*n3;   T[5,4] = l2*m3+l3*m2;  T[5,5] = m2*n3+m3*n2;  T[5,6] = n2*l3+n3*l2;
    T[6,1] =   2*l3*l1;  T[6,2] =   2*m3*m1;  T[6,3] =   2*n3*n1;   T[6,4] = l3*m1+l1*m3;  T[6,5] = m3*n1+m1*n3;  T[6,6] = n3*l1+n1*l3;
end



function setB2(elem::ShellDegenerated, R::Matx, ip::Ip, dNdX::Matx, N::Vect, B::Matx)
    nnodes, ndim = size(dNdX)
    B .= 0.0
    t = elem.mat.t
    #ζ = ip.R[3] #ip.coord[3] não sei se está certo
    #@show ζ
  
    # C = getcoords(elem)

    # V1=zeros(3,nnodes)
    # V2=zeros(3,nnodes)
    # V3=zeros(3,nnodes)

    for i in 1:nnodes

        # J (3x2)
        # V3[:,i] = cross(vec(J[:,1]), vec(J[:,2]))
        # normalize!(V3[:,i])
        # V1[:,i] = cross(vec([1 0 0]),V3[:,i])
        # V2[:,i] = cross(V3[:,i],V1[:,i])

        l1, m1, n1 = R[:,1]
        l2, m2, n2 = R[:,2]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]

        ζ = ip.R[3] #ip.coord[3] não sei se está certo
        
        j    = i-1

        B[1,1+j*ndim] = dNdx;  B[1,4+j*ndim] = -ζ*dNdx*t/2*l2;  B[1,5+j*ndim] = -ζ*dNdx*t/2*l1

        B[2,1+j*ndim] = dNdy;  B[2,4+j*ndim] = -ζ*dNdy*t/2*l2;  B[2,5+j*ndim] = -ζ*dNdy*t/2*l1

                               B[3,4+j*ndim] = -N[i]*t/2*l2;    B[3,5+j*ndim] =  N[i]*t/2*l1

        B[4,2+j*ndim] = dNdx;  B[4,4+j*ndim] = -ζ*dNdx*t/2*m2;  B[4,5+j*ndim] = -ζ*dNdx*t/2*m1                      

        B[5,2+j*ndim] = dNdy;  B[5,4+j*ndim] = -ζ*dNdy*t/2*m2;  B[5,5+j*ndim] = -ζ*dNdy*t/2*m1  

                               B[6,4+j*ndim] = -N[i]*t/2*m2;    B[6,5+j*ndim] =  N[i]*t/2*m1   

        B[7,2+j*ndim] = dNdx;  B[7,4+j*ndim] = -ζ*dNdx*t/2*n2;  B[7,5+j*ndim] = -ζ*dNdx*t/2*n1                     

        B[8,2+j*ndim] = dNdy;  B[8,4+j*ndim] = -ζ*dNdy*t/2*n2;  B[8,5+j*ndim] = -ζ*dNdy*t/2*n1

                               B[9,4+j*ndim] = -N[i]*t/2*n2;    B[9,5+j*ndim] =  N[i]*t/2*n1
         
    end

end

function setB(elem::ShellDegenerated, R::Matx, ip::Ip, dNdX::Matx, B::Matx)
    nnodes, ndim = size(dNdX)
    t = elem.mat.t
    # B .= 0.0

    l1, m1, n1 = R[:,1]
    l2, m2, n2 = R[:,2]
    ndof = 5

    # @show dNdX
    # error()

    for i in 1:nnodes
        c = (i-1)*ndof
        ζ = ip.R[3]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        dNdz = dNdX[i,3]

        # dNdz = N[i]*(2/t)  # t/2 is the jacobian in the normal direction
        # dNdz = dot(dNdR[i,:], pinv(R[:,3]*t/2))

        B[1,1+c] = dNdx;                                    B[1,4+c] = -ζ*dNdx*t/2*l2;            B[1,5+c] = ζ*dNdx*t/2*l1
                         B[2,2+c] = dNdy;                   B[2,4+c] = -ζ*dNdy*t/2*m2;            B[2,5+c] = ζ*dNdy*t/2*m1
                                          B[3,3+c] = dNdz;  B[3,4+c] = -ζ*dNdz*t/2*n2;            B[3,5+c] = ζ*dNdz*t/2*n1
                         B[4,2+c] = dNdz; B[4,3+c] = dNdy;  B[4,4+c] = -ζ*t/2*(dNdy*l2+dNdx*m2);  B[4,5+c] = ζ*t/2*(dNdy*l1+dNdx*m1)
        B[5,1+c] = dNdz;                  B[5,3+c] = dNdx;  B[5,4+c] = -ζ*t/2*(dNdz*m2+dNdy*n2);  B[5,5+c] = ζ*t/2*(dNdz*m1+dNdy*n1)
        B[6,1+c] = dNdy; B[6,2+c] = dNdx;                   B[6,4+c] = -ζ*t/2*(dNdz*l2+dNdx*n2);  B[6,5+c] = ζ*t/2*(dNdz*l1+dNdx*n1)

    end

end


function Dmatrix(elem::ShellDegenerated)

    nu = elem.mat.nu
    E1 = elem.mat.E/(1-elem.mat.nu^2)
    G  = elem.mat.E/(2*(1+elem.mat.nu))
    G1 = 5/6*G

              D =   [E1   nu*E1  0  0   0
                    nu*E1    E1  0  0   0
                     0        0  G  0   0
                     0        0  0  G1  0
                     0        0  0  0  G1 ]

    return D
end


function setD(elem::ShellDegenerated, D::Matx)

    nu = elem.mat.nu
    E1 = elem.mat.E/(1-elem.mat.nu^2)
    G  = elem.mat.E/(2*(1+elem.mat.nu))
    G1 = 5/6*G

    D .=  [   E1  nu*E1  0  0  0  0
           nu*E1     E1  0  0  0  0
            0         0  0  0  0  0
            0         0  0  G  0  0
            0         0  0  0  G1 0
            0         0  0  0   0 G1 ]

end


function elem_config_dofs(elem::ShellDegenerated)
    ndim = elem.env.ndim
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
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    t = elem.mat.t

    C = getcoords(elem)
    K = zeros(5*nnodes, 5*nnodes)
    B = zeros(6, 5*nnodes)

    D  = Array{Float64}(undef, 6, 6)
    # DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    R = zeros(3,3)
    T = zeros(6,6)

    # D = Dmatrix(elem)
    setD(elem, D)


    for ip in elem.ips
        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        J2D = C'*dNdR 
        rot_matrix_R(elem, J2D, R)
        rot_matrix_T(elem, R, T)

        J = [ J2D  R[:,3]*t/2 ]

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX = dNdR*inv(J)

        setB(elem, R, ip, dNdX, B)
        detJ = det(J)

        coef = detJ*ip.w*0.5  # check 0.5
        K += coef*B'*T'*D*T*B

    end

    keys =(:ux, :uy, :uz, :rx, :ry)
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    map = elem_map(elem)
    return K, map, map
end

function elem_stiffness2(elem::ShellDegenerated)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    K = zeros(5*nnodes, 5*nnodes)
    B = zeros(9, 5*nnodes)
    JJ  = zeros(9,9)
    J  = Array{Float64}(undef, ndim, ndim)
    V  = Array{Float64}(undef, ndim, ndim)
    R = zeros(3,3)
    T = zeros(6,6)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    D = Dmatrix(elem)

    C = getcoords(elem)
    if size(C,2)==2
        cxyz  = zeros(8,3)
        cxyz[:,1:2]  = C
    else
        cxyz = C
    end
    #@show cxyz
    # @show  elem.ips

    for ip in elem.ips
        # elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute B matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        #@gemm J = C'*dNdR
        #@gemm dNdX = dNdR*inv(J)
        J = cxyz'*dNdR  #J (3x2)
        #dNdX = dNdR*inv(J)
        dNdX = dNdR*pinv(J) #! correto

        rot_matrix_R(elem, J, R)
        rot_matrix_T(elem, R, T)
        # @show T
        #@show V
        #@show 
        
        # aux = (elem.mat.t/2)*cross(J[:,1],J[:,2])
        #@show aux

        # normalize!(aux)
        #@show aux
        #falta normalizar o vetor     


        detJ = norm2(J)
        #detJ = det(J)
        #@show detJ

        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

           
        setB(elem, J, ip, dNdX, N, B)  # 6x40

        H = [1 0 0 0 0 0 0 0 0
             0 0 0 0 1 0 0 0 0
             0 0 0 0 0 0 0 0 1
             0 1 0 1 0 0 0 0 0
             0 0 0 0 0 1 0 1 0
             0 0 1 0 0 0 1 0 0]

        #J3= [J (elem.mat.t/2)*cross(vec(J[:,1]), vec(J[:,2])) ]
        J3= [J normalize(cross(vec(J[:,1]), vec(J[:,2]))) ]

           for i in 1:3
                JJ[(i-1)*3+1:i*3, (i-1)*3+1:i*3] = inv(J3)
           end
        
           @show J3
           @show JJ
           error()

        B1 = H*JJ*B
        B2 = T*B1 # 6x40
        B3 = [B2[1:3,:]
              B2[5:6,:]]  # 5x40

        #@show T
        #@show B3

       # coef = detJ*ip.w*th
        coef = 0.5*detJ*ip.w  # 0.5 because in zeta goes from -1 to +1
        #@show coef

        K += (B3'*D*B3)*coef

        # @showm D
        # @showm B3'*B3
    end
    
    # @show K

     map = elem_map(elem)
    return K, map, map
    
end


function elem_update!(elem::ShellDegenerated, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end


#=
function elem_update!(elem::ShellDegenerated, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
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
        if elem.env.modeltype=="axisymmetric"
            th = 2*pi*ip.coord.x
        end
        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        setB(elem, ip, dNdX, B)
        @gemv Δε = B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return failure("ShellDegenerated: Error at integration point $(ip.id)")
        #if failed(status)
            #status.message = "ShellDegenerated: Error at integration point $(ip.id)\n" * status.message
            #return status
        #end
        coef = detJ*ip.w*th
        @gemv dF += coef*B'*Δσ
    end
    F[map] += dF
    return success()
end
function elem_vals(elem::ShellDegenerated)
    vals = OrderedDict{Symbol,Float64}()
    if haskey(ip_state_vals(elem.mat, elem.ips[1].state), :damt)
        mean_dt = mean( ip_state_vals(elem.mat, ip.state)[:damt] for ip in elem.ips )
        vals[:damt] = mean_dt
        mean_dc = mean( ip_state_vals(elem.mat, ip.state)[:damc] for ip in elem.ips )
        vals[:damc] = mean_dc
    end
    #vals = OrderedDict{String, Float64}()
    #keys = elem_vals_keys(elem)
#
    #dicts = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
    #nips = length(elem.ips)
#
    #for key in keys
        #s = 0.0
        #for dict in dicts
            #s += dict[key]
        #end
        #vals[key] = s/nips
    #end
    return vals
end
=#
