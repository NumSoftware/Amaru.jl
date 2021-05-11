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


# the strain-displacement matrix for membrane forces
function D_matrixm(elem::ShellQUAD4)

    coef1 = elem.mat.thick*elem.mat.E/(1-elem.mat.nu^2)
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

    D_matb = D_matm*(elem.mat.thick^2/12)

    return D_matb
end

# the strain-displacement matrix for shear forces

function D_matrixs(elem::ShellQUAD4)

    coef = elem.mat.thick*(5/6)*elem.mat.E/(2*(1+elem.mat.nu))

            D_mats = [coef    0
                        0     coef];
    return D_mats
end


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
    B .= 0.0

    for i in 1:nnodes
        dNdx = dNdX[1,i]
        dNdy = dNdX[2,i]
        Ni = N[i]

        j    = i-1

        # matrix Bm
        B[1,1+j*ndim] = dNdx
        B[2,2+j*ndim] = dNdy
        B[3,1+j*ndim] = dNdy
        B[3,2+j*ndim] = dNdx

        # matrix Bb
        B[4,4+j*ndim] = -dNdx
        B[5,5+j*ndim] = -dNdy
        B[6,4+j*ndim] = -dNdy
        B[7,5+j*ndim] =  dNdx

        # matrix Bs
        B[7,3+j*ndim] = dNdx
        B[7,4+j*ndim] = -Ni
        B[8,3+j*ndim] = dNdy
        B[8,4+j*ndim] = -Ni

    end
end



function elem_stiffness(elem::ShellQUAD4)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    ndof   = 5

    C = get_coords(elem)
    R = RotMatrix(elem)
    C = (C*R')[:,1:2]
    # B = zeros(ndim, nnodes*ndim)
    K = zeros(nnodes*ndof, nnodes*ndof)

    # DB = zeros(ndim, nnodes*ndim)
    # J  = zeros(ndim-1, ndim)

    D_matm = D_matrixm(elem)
    D_matb = D_matrixb(elem)
    D_mats = D_matrixs(elem)

        for ip in elem.ips
        
            # compute shape Jacobian
            N    = elem.shape.func(ip.R)
            dNdR = elem.shape.deriv(ip.R)

            J = dNdR*C
            detJ = norm2(J)
            dNdX = inv(J)*dNdR

            B = zeros(8,ndof*nnodes)
            setB(elem, N, dNdX, B)


            # COMO MONTAR A MATRIZ D?

            #D = [D_matm D_matb D_mats]
            println(B)
            
        
         end

        

        #keys = (:ux, :uy, :uz, :rx, :ry)[1:ndim]
         keys = (:ux, :uy, :uz, :rx, :ry)

        map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
        #map  = elem_map(elem)
     
        return K, map, map
end

function elem_update!(elem::ShellQUAD4, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)

    dU  = U[map]
    F[map] += K*dU
    return success()
end
