# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ShellQUAD4

mutable struct ShellQUAD4<:Mechanical
    id    ::Int
    shape ::ShapeType   #CellShape
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

#=
function RotMatrix(elem::ShellQUAD4)

    v12 = zeros(3,1)
    v13 = zeros(3,1)
    vxe = zeros(3,1)
    vye = zeros(3,1)
    vze = zeros(3,1)

    cxyz = get_coords(elem)

        if size(cxyz,2)==2
            v12[3] = 0
            v13[3] = 0
        else
            v12[3] = cxyz[2,3] - cxyz[1,3]
            v13[3] = cxyz[3,3] - cxyz[1,3]
        end

    v12[1] = cxyz[2,1] - cxyz[1,1]
    v12[2] = cxyz[2,2] - cxyz[1,2]

    v13[1] = cxyz[3,1] - cxyz[1,1]
    v13[2] = cxyz[3,2] - cxyz[1,2]

    vze[1] = v12[2]*v13[3] - v12[3]*v13[2]
    vze[2] = v12[3]*v13[1] - v12[1]*v13[3]
    vze[3] = v12[1]*v13[2] - v12[2]*v13[1]

    dz = sqrt(vze[1]^2 + vze[2]^2 + vze[3]^2);

  # Unit vector normal to element surface
    vze[1] = vze[1]/dz
    vze[2] = vze[2]/dz
    vze[3] = vze[3]/dz

  # XZ plane intesection with element surface
    vxe[1] =  1/sqrt(1+(vze[1]/vze[3])^2)
    vxe[2] =  0
    vxe[3] = -1/sqrt(1+(vze[3]/vze[1])^2)

    dd = vxe[1]*vze[1] + vxe[3]*vze[3];
    if (abs(dd) > 1e-8)
      vxe[3] = -vxe[3]
    end

    if ((vze[3] == 0) && (vze[1] == 0))
      vxe[1] =  1
      vxe[2] =  0
      vxe[3] =  0
    end

  # Vector product
    vye[1] = vze[2]*vxe[3] - vxe[2]*vze[3]
    vye[2] = vze[3]*vxe[1] - vxe[3]*vze[1]
    vye[3] = vze[1]*vxe[2] - vxe[1]*vze[2]

    dy = sqrt(vye[1]^2 + vye[2]^2 + vye[3]^2)
    vye[1] = vye[1]/dy
    vye[2] = vye[2]/dy
    vye[3] = vye[3]/dy

    if (vye[2] < 0 )
      vye[1] = -vye[1]
      vye[2] = -vye[2]
      vye[3] = -vye[3]
      vxe[1] = -vxe[1]
      vxe[2] = -vxe[2]
      vxe[3] = -vxe[3]
    end

    Rot = [ vxe[1] vxe[2] vxe[3]
           vye[1] vye[2] vye[3]
           vze[1] vze[2] vze[3] ]

    return Rot
end
=#

function RotMatrix(elem::ShellQUAD4)

    C = get_coords(elem)

    R = [ 0.0, 0.0, 0.0 ]
    dNdR = elem.shape.deriv(R)
    J = dNdR*C

    if elem.env.ndim==2
        return Matrix{Float64}(I, 3,3) # identity matrix
    else
        L1 = vec(J[1,:])
        L2 = vec(J[2,:])
        L3 = cross(L1, L2)  # L3 is normal to the first element facet
        L2 = cross(L3, L1)
        normalize!(L1)
        normalize!(L2)
        normalize!(L3)
        return [L1 L2 L3]
    end
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

function setB(elem::ShellQUAD4, N::Vect, dNdX::Matx, B::Matx)
    ndim, nnodes = size(dNdX)
    ndof = 5
    B .= 0.0

    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]

        j    = i-1

        # matrix Bm
        B[1,1+j*ndof] = dNdx
        B[2,2+j*ndof] = dNdy
        B[3,1+j*ndof] = dNdy
        B[3,2+j*ndof] = dNdx

        # matrix Bb
        B[4,4+j*ndof] =  dNdx  # -dNdx
        B[5,5+j*ndof] =  dNdy  # -dNdy
        B[6,4+j*ndof] =  dNdy  # -dNdy
        B[6,5+j*ndof] =  dNdx  # B[7,5+j*ndof]
        
        # matrix Bs
        B[7,3+j*ndof] = dNdx
        B[7,4+j*ndof] = -N[i]
        B[8,3+j*ndof] = dNdy
        B[8,5+j*ndof] = -N[i]   # B[8,4+j*ndof] = -N[i]

    end
end

function setD(elem::ShellQUAD4, D::Matx)
    E  = elem.mat.E
    nu = elem.mat.nu
    t  = elem.mat.t
    α  = 5/6

    d11 = d22 = E/(1-nu^2)
    d12 = d21 = E*nu/(1-nu^2)
    d33 = E/(1-nu^2)*(1-nu)/2
    d44 = d55 = E/(2*(1+nu))*α  #E/2/(1+nu)*α

    D .= [ d11*t d12*t 0     0 0 0 0 0
          d21*t d22*t 0     0 0 0 0 0
          0     0     d33*t 0 0 0 0 0
          0     0     0     d11*t^3/12 d12*t^3/12 0          0 0
          0     0     0     d21*t^3/12 d22*t^3/12 0          0 0
          0     0     0     0          0          d33*t^3/12 0 0
          0     0     0     0          0          0          d44*t 0
          0     0     0     0          0          0          0     d55*t ]

end

function elem_stiffness(elem::ShellQUAD4)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    ndof   = 5
    
    C = get_coords(elem)
    R = RotMatrix(elem)
    C = (C*R')[:,1:2]
    B = zeros(8, nnodes*ndof)
    D = zeros(8, 8)
    K = zeros(nnodes*ndof, nnodes*ndof)


    for ip in elem.ips
    
        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        J = dNdR*C
        detJ = norm2(J)
        dNdX = inv(J)*dNdR

        setB(elem, N, dNdX, B)
        setD(elem, D)

        coef = detJ*ip.w
         
        # @show size(B)
        # @show size(D)
        
        K += B'*D*B*coef
        # @show size(K)

    end

    #K = R*D*R'

    keys = (:ux, :uy, :uz, :rx, :ry)

    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
     
    return K, map, map
end

function elem_update!(elem::ShellQUAD4, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)

    dU  = U[map]
    F[map] += K*dU
    return success()
end