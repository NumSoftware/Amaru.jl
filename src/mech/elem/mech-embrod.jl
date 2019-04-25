# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechEmbRod<:Mechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    # specific fields
    cache_NN::Array{Float64,2}

    MechEmbRod() = new()
end

matching_shape_family(::Type{MechEmbRod}) = LINE_SHAPE

# Returns the corresponding embedded element type
function embedded_elem_type(::Type{MechRod})
    return MechEmbRod
end

function elem_config_dofs(elem::MechEmbRod)
    # No-op function.
    # The solid linked element will set the required dofs.
end

function elem_map(elem::MechEmbRod)
    ndim = elem.env.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    solid = elem.linked_elems[1]
    return [ node.dofdict[key].eq_id for node in solid.nodes for key in keys ]
end

function mountNN(elem::MechEmbRod)
    ndim = elem.env.ndim
    solid = elem.linked_elems[1]
    n  = length(solid.nodes)
    m  = length(elem.nodes)
    NN = zeros(ndim*n, ndim*m)
    Cs = elem_coords(solid)

    for j=1:m 
        R = inverse_map(solid.shape, Cs, elem.nodes[j].X)
        N = solid.shape.func(R)
        for i=1:n
            for k=1:ndim
                NN[(i-1)*ndim+k, (j-1)*ndim+k] = N[i]
            end
        end
    end
    return NN
end


function elem_init(elem::MechEmbRod)
    elem.cache_NN = mountNN(elem)
    return nothing
end

function elem_stiffness(elem::MechEmbRod)
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A = elem.mat.A
    C = elem_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef,1, ndim)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j=1:ndim
                B[1,j+(i-1)*ndim] = dNdR[1,i]*J[j]/detJ^2.0
            end
        end

        E    = calcD(elem.mat, ip.data)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end

    NN = elem.cache_NN
    map = elem_map(elem)
    return NN*K*NN', map, map
end

function elem_update!(elem::MechEmbRod, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    NN  = elem.cache_NN
    map = elem_map(elem)
    dU  = U[map]
    dUr = NN'*dU

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A = elem.mat.A

    dF = zeros(nnodes*ndim)
    C  = elem_coords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array{Float64}(undef,1, ndim)
    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j=1:ndim
                B[1,j+(i-1)*ndim] = dNdR[1,i]*J[j]/detJ^2.0
            end
        end

        deps = (B*dUr)[1]
        dsig = stress_update(elem.mat, ip.data, deps)
        coef = A*detJ*ip.w
        dF  += coef*B'*dsig
    end

    F[map] += NN*dF
end

