# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ShellQUAD4

mutable struct ShellQUAD4<:Mechanical
    id    ::Int
    shape ::CellShape
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

matching_shape_family(::Type{ShellQUAD4}) = SOLID_CELL

function distributed_bc(elem::ShellQUAD4, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
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

# the strain-displacement matrix for membrane forces
function Dm_maxtrix(elem::ShellQUAD4)

    coef1 = elem.mat.t*elem.mat.E/(1-elem.mat.nu^2)
    coef2 = elem.mat.nu*coef1
    coef3 = coef1*(1-elem.mat.nu)/2

        Dm = [coef1  coef2 0
                  coef2  coef1 0
                  0      0     coef3]
    return Dm
end

# the strain-displacement matrix for bending moments
function Db_maxtrix(elem::ShellQUAD4)

    Dm = Dm_maxtrix(elem)

    Db = Dm*(elem.mat.t^2/12)

    return Db
end

# the strain-displacement matrix for shear forces
function Ds_maxtrix(elem::ShellQUAD4)

    coef = elem.mat.t*(5/6)*elem.mat.E/(2*(1+elem.mat.nu))

            Ds = [coef    0
                        0     coef]
    return Ds
end

# Rotation Matrix
function RotMatrix(elem::ShellQUAD4, J::Matrix{Float64})
    
    Z = zeros(1,2) # zeros(2,1)

    # artifice for mounting the rotation matrix for flat elements
    if size(J,1)==2
        J = [J
             Z]
    else
        J = J
    end
    
    L1 = vec(J[:,1])
    L2 = vec(J[:,2])
    L3 = cross(L1, L2) 
    L2 = cross(L1, L3)
    normalize!(L1)
    normalize!(L2)
    normalize!(L3)

    Z1 = zeros(1,2)

    Rot = [ L1' Z1 # changing the position of L2 by L1 in the rotation matrix to be equivalent that one used by Oñate and Aldo
    L2' Z1
    L3' Z1
    Z1   -L2'
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
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) 
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
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
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
        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        j    = i-1

        Bm[1,1+j*ndof] = dNdx  
        Bm[2,2+j*ndof] = dNdy   
        Bm[3,1+j*ndof] = dNdy 
        Bm[3,2+j*ndof] = dNdx 

    end
end

function setBs_bar(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bs_bar::Matx)
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
            # dimension 8x12
    end
end

function elem_stiffness(elem::ShellQUAD4)

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
    Kelem = zeros( nnodes*5 , nnodes*5 )

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

        #Assembly of substitutive B matrix

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

                    # Rotação dos termos  bmat_ss[1:2,1] bmat_ss[1:2,4] bmat_ss[1:2,7] bmat_ss[1:2,10]
                    bmat_s  = zeros(2,20)
                    for i in 1:nnodes
                        bb = [0  0 bmat_ss[1, (i-1)*3+1]
                              0  0 bmat_ss[2, (i-1)*3+1]]

                        bmat_s[:, (i-1)*5+1:i*5] = [ bb*Ri[1:3, 1:3]   bmat_ss[:,(i-1)*3+2:(i-1)*3+3] ]
                    end

                     coef = detJ1*ip.w

                     Kb =    Bb'*Db*Bb*coef # Pq não rotaciona?
                     Km = R'*Bm'*Dm*Bm*R*coef 
                     Ks = bmat_s'*Ds*bmat_s*coef 

                     Kelem += (Kb + Km + Ks)
            end
        map = elem_map(elem) 
    return Kelem, map, map
end

function elem_update!(elem::ShellQUAD4, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end
