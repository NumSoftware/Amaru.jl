# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ShellQUAD4

mutable struct ShellQUAD4<:Mechanical
    id    ::Int
    shape ::ShapeType
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function ShellQUAD4()
        return new()
    end
end

matching_shape_family(::Type{ShellQUAD4}) = SOLID_SHAPE

function distributed_bc(elem::ShellQUAD4, facet::Union{Facet, Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
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
    C = get_coords(nodes, ndim)

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

# the strain-displacement matrix for membrane forces
function D_matrixm(elem::ShellQUAD4)

    coef1 = elem.mat.t*elem.mat.E/(1-elem.mat.nu^2)
    coef2 = elem.mat.nu*coef1
    coef3 = coef1*(1-elem.mat.nu)/2

        D_matm = [coef1  coef2 0
                  coef2  coef1 0
                  0      0     coef3];
    return D_matm
end

# the strain-displacement matrix for bending moments
function D_matrixb(elem::ShellQUAD4)

    D_matm = D_matrixm(elem)

    D_matb = D_matm*(elem.mat.t^2/12)

    return D_matb
end

# the strain-displacement matrix for shear forces

function D_matrixs(elem::ShellQUAD4)

    coef = elem.mat.t*(5/6)*elem.mat.E/(2*(1+elem.mat.nu))

            D_mats = [coef    0
                        0     coef];
    return D_mats
end

# Rotation Matrix

function RotMatrix(elem::ShellQUAD4, J::Matrix{Float64})
    
    Z = zeros(2,1)

    if size(J,2)==2
        J = [J Z]
    else
        J = J
    end
    
    L1 = vec(J[1,:])
    L2 = vec(J[2,:])
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

function elem_config_dofs(elem::ShellQUAD4)
    ndim = elem.env.ndim
    ndim == 1 && error("ShellQUAD4: Shell elements do not work in 1d analyses")
    #if ndim==2
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
        end
    #else
        #error("ShellQUAD4: Shell elements do not work in this analyses")
        #=
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
        =#
    #end
end

function elem_map(elem::ShellQUAD4)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:ux, :uy, :uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function setBb(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bb::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)
    ndof = 5
    Bb .= 0.0
   
    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        j    = i-1

        Bb[1,4+j*ndof] = -dNdx  

        Bb[2,5+j*ndof] = -dNdy   

        Bb[3,4+j*ndof] = -dNdy 

        Bb[3,5+j*ndof] = -dNdx 

    end
end

function setBm(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bm::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)
    ndof = 5
    Bm .= 0.0
   
    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        j    = i-1

        Bm[1,1+j*ndof] = dNdx  

        Bm[2,2+j*ndof] = dNdy   

        Bm[3,1+j*ndof] = dNdy 

        Bm[3,2+j*ndof] = dNdx 

    end
end

function setBs(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bs::Matx)
    nnodes = length(elem.nodes)
    ndof = 5
    Bs .= 0.0

    for i in 1:nnodes

        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        j    = i-1
        
        Bs[1,3+j*ndof] = dNdx 
        Bs[1,4+j*ndof] = -N[i]  
        Bs[2,3+j*ndof] = dNdy  
        Bs[2,5+j*ndof] = -N[i]

    end
end

function elem_stiffness(elem::ShellQUAD4)

    nnodes = length(elem.nodes)

    D_matm = D_matrixm(elem)
    D_mats = D_matrixs(elem)
    D_matb = D_matrixb(elem)

    Bb = zeros(3, nnodes*5)
    Bm = zeros(3, nnodes*5)
    Bs = zeros(2, nnodes*5)

    nr = 5   
    nc = 5
    R = zeros(nnodes*nr, nnodes*nc)
    c     = zeros(8,8)
    b_bar = zeros(8,12)
    K_elem = zeros( nnodes*5 , nnodes*5 )

    C = get_coords(elem)

    for ip in elem.ips

        if size(C,2)==2
            cxyz  = zeros(4,3)
            cxyz[:,1:2]  = C
        else
            cxyz = C
        end

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        #@show size(dNdR)

        J = dNdR*cxyz
        #@show size(J)
        #@showm J
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        Ri = RotMatrix(elem, J)
        Ri′ = Ri[1:2, 1:3]
  
        ctxy = cxyz*Ri[1:3, 1:3]' # Rotate coordinates to element mid plane

        J1 = dNdR*ctxy
        #@showm J1
        
        dNdX = pinv(J)*dNdR
        #@show size(dNdX)

        dNdX′ = Ri′*dNdX
        #@show size(dNdX′)
              
        for i in 1:nnodes
            R[(i-1)*nr+1:i*nr, (i-1)*nc+1:i*nc] = Ri
        end

        for igaus = 1 : 4
            
        setBb(elem, N, dNdX′, Bb)
        #@show size(Bb)

        setBm(elem, N, dNdX′, Bm)
        #@showm Bm
        #@show size(Bm)

        setBs(elem, N, dNdX′, Bs)
        #@showm Bs
        #@show size(Bs)
        
        invJ1  =  pinv(J1) 
        #@showm invJ1
        detJ1 = norm2(J1)
        detJ1 > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        for i in 1:nnodes
            c[(i-1)*2+1:i*2, (i-1)*2+1:i*2] = J1[1:2,1:2]
        end

        bmat_s = [Bs[1:2, 3:5] Bs[1:2, 8:10] Bs[1:2, 13:15] Bs[1:2, 18:20]]

        for i in 1:nnodes
            b_bar[2*i-1:2*i,:] = bmat_s[1:2,:]
        end

        T_mat = [ 1  0  0  0  0  0  0  0
                  0  0  0  1  0  0  0  0
                  0  0  0  0  1  0  0  0
                  0  0  0  0  0  0  0  1 ]

        P_mat = [ 1  -1   0   0
                  0   0   1   1
                  1   1   0   0
                  0   0   1  -1 ];

        A_mat = [ 1  ip.w  0    0   
                  0    0  1  ip.w ];  

        bmat_ss = invJ1[1:2, 1:2] * A_mat * inv(P_mat) * T_mat * c * b_bar

        bmat_s1 = [0  0 bmat_ss[1, 1]
                   0  0 bmat_ss[2, 1]];

        bmat_s2 = [0  0 bmat_ss[1, 4]
                   0  0 bmat_ss[2, 4]];

        bmat_s3 = [0  0 bmat_ss[1, 7]
                   0  0 bmat_ss[2, 7]];

        bmat_s4 = [0  0 bmat_ss[1,10]
                   0  0 bmat_ss[2,10]];

        bmat_s1 = [bmat_s1*Ri[1:3, 1:3] bmat_ss[:,2:3]];
        bmat_s2 = [bmat_s2*Ri[1:3, 1:3] bmat_ss[:,5:6]];
        bmat_s3 = [bmat_s3*Ri[1:3, 1:3] bmat_ss[:,8:9]];
        bmat_s4 = [bmat_s4*Ri[1:3, 1:3] bmat_ss[:,11:12]];
        
        bmat_s = [bmat_s1 bmat_s2 bmat_s3 bmat_s4]
        #@show size(bmat_s)

        coef = detJ1*ip.w

        K_b = Bb'*D_matb*Bb*coef;                
        K_m = R'*Bm'*D_matm*Bm*R*coef;
        K_s = bmat_s'*D_mats*bmat_s*coef;

        K_elem += (K_b + K_m + K_s)
    end

        map = elem_map(elem) 

    return K_elem, map, map
  end
end

function elem_update!(elem::ShellQUAD4, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end
