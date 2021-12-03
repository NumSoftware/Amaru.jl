# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ShellQUAD4

mutable struct ShellQUAD4<:Mechanical
    id    ::Int
    shape ::ShapeType # CellShape
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
    th    = elem.env.thickness
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



function matrixR(elem::ShellQUAD4, J::Matrix{Float64})
    L1 = vec(J[1,:])
    L2 = vec(J[2,:])
    L3 = cross(L1, L2)  # L1 is normal to the first element face
    L2 = cross(L3, L1)
    normalize!(L1)
    normalize!(L2)
    normalize!(L3)
    Z = zeros(1,3)

    return [ L1' Z
             L2' Z
             L3' Z
             Z   L1'
             Z   L2'
             Z   L3']

             
end


function elem_config_dofs(elem::ShellQUAD4)
    ndim = elem.env.ndim
    ndim == 1 && error("ShellQUAD4: Shell elements do not work in 1d analyses")
    
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end


function setBm(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bm::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)
    ndof = 6 #5
    Bm .= 0.0
    k =  1e-8

    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        j    = i-1

        Bm[1,1+j*ndof] = dNdx

        Bm[2,2+j*ndof] = dNdy

        Bm[3,1+j*ndof] = 0.5*dNdy
        Bm[3,2+j*ndof] = 0.5*dNdx

        Bm[4,3+j*ndof] = 0.5*dNdx
        Bm[4,5+j*ndof] = 0.5*N[i]  

        Bm[5,3+j*ndof] = 0.5*dNdy
        Bm[5,4+j*ndof] = -0.5*N[i] 

        Bm[6,6+j*ndof] = k*N[i]

    end
end

function setBb(elem::ShellQUAD4, N::Vect, dNdX::Matx, Bb::Matx)
    nnodes = length(elem.nodes)
    # ndim, nnodes = size(dNdX)
    ndof = 6 #5
    Bb .= 0.0
   
    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        j    = i-1

        Bb[1,5+j*ndof] = dNdx  

        Bb[2,4+j*ndof] = -dNdy   

        Bb[3,4+j*ndof] =  -0.5*dNdy 

        Bb[3,5+j*ndof] = 0.5*dNdx 

    end
end


function setDm(elem::ShellQUAD4, Dm::Matx)
    E  = elem.mat.E
    nu = elem.mat.nu
    d = 1 #1e-5
    
    Dm .= E/(1-nu^2)*[ 1 nu 0 0 0 0
                       nu 1 0 0 0 0
                      0 0 1-nu 0 0 0
                      0 0 0 5/6*(1-nu) 0 0
                      0 0 0 0 5/6*(1-nu) 0
                      0 0 0 0 0 d]

end

function setDb(elem::ShellQUAD4, Db::Matx)
    E  = elem.mat.E
    nu = elem.mat.nu

    Db .= E/(1-nu^2)*[ 1 nu 0
                       nu 1 0
                        0 0 1-nu]

end

function elem_stiffness(elem::ShellQUAD4)
    nnodes = length(elem.nodes)
    ndofg   = 6
    ndofl   = 6 #5

    t  = elem.mat.t

    C = get_coords(elem)
    # R = RotMatrix(elem)
    #@show size(C)


    Bm = zeros(6, nnodes*ndofl)
    Bb = zeros(3, nnodes*ndofl)
    #@showm Bm
    #@showm Bb
    Dm = zeros(6, 6)
    Db = zeros(3, 3)
    K = zeros(nnodes*ndofg, nnodes*ndofg)
    nr = 6 #5  
    nc = 6 
    R = zeros(nnodes*nr, nnodes*nc)

    for ip in elem.ips
        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        #@show size(dNdR)
        
        J = dNdR*C  #dNdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        
        Ri = matrixR(elem, J)
        #display(Ri)

        for i in 1:nnodes
            R[(i-1)*nr+1:i*nr, (i-1)*nc+1:i*nc] = Ri
        end


        dNdX = pinv(J)*dNdR
        # (3x4) = (3x2) (2x4)

        #@showm J
        #@showm pinv(J)

        #@showm dNdR
        #@showm dNdX

        Ri′ = Ri[1:2, 1:3]
        #(2X3)
        #@show size(Ri′)

        dNdX′ = Ri′*dNdX
        #@showm dNdX′

        #(2x4) = (2x3)*(3x4)

        setBm(elem, N, dNdX′, Bm)
        #@showm Bm
        setBb(elem, N, dNdX′, Bb)
        #@showm Bb

        setDm(elem, Dm)
        #@showm Dm
        setDb(elem, Db)
        #@showm Db

        coef = detJ*ip.w
        # @show size(K)
        # @show size(B)

        K1 = t*(Bm'*Dm*Bm*coef)
        K2 = (t^3/12)*(Bb'*Db*Bb*coef)
        #@show K1
        #@show K2

        K += R'*(K1+K2)*R
        #K += R'*(K1)*R
        #K += R'*(K2)*R
        # (24x24) = (24x24) (24x24 + 24x24) 24x24
        #@show K
    end

    #@showm K

    keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
     
    return K, map, map
end
 
function elem_update!(elem::ShellQUAD4, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)

    dU  = U[map]
    F[map] += K*dU
    return success()
end