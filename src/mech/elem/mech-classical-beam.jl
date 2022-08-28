# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechClassicalBeam

mutable struct MechClassicalBeam<:Mechanical
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function MechClassicalBeam()
        return new()
    end
end

matching_shape_family(::Type{MechClassicalBeam}) = LINECELL

function beam_shape_func(x::Float64, L::Float64)
    N = Array{Float64}(undef,6)
    N[1] = 1 - x/L
    N[2] = 1 - 3*x^2/L^2 + 2*x^3/L^3
    N[3] = x - 2*x^2/L + x^3/L^2
    N[4] = x/L
    N[5] = 3*x^2/L^2 - 2*x^3/L^3
    N[6] = x^3/L^2 - x^2/L

    return N
end

function beam_second_deriv(ξ::Float64, nnodes::Int)
    if nnodes==2
        DD = Array{Float64}(undef,4)
    else
        DD = Array{Float64}(undef,6)
    end
    return DD
end

function elem_config_dofs(elem::MechClassicalBeam)
    ndim = elem.env.ndim
    ndim == 1 && error("MechClassicalBeam: Beam elements do not work in 1d analyses")
    if ndim==2
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :rz, :mz)
        end
    else
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
    end
end

function elem_map(elem::MechClassicalBeam)::Array{Int,1}
    if elem.env.ndim==2
        dof_keys = (:ux, :uy, :rz)
    else
        dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    end
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

# Return the class of element where this material can be used
#client_shape_class(mat::MechClassicalBeam) = LINECELL

function calcT(elem::MechClassicalBeam, C)
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,1])/L
    return
end

function distributed_bc(elem::MechClassicalBeam, facet::Nothing, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim

    # Check bcs
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")
    !(key in (:tx, :ty, :tl, :tn)) && error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")

    # target = facet!=nothing ? facet : elem
    nodes  = elem.nodes
    nnodes = length(nodes)
    t      = elem.env.t
    A      = elem.mat.A
    

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = getcoords(nodes, ndim)
    L = norm(C[2,:]-C[1,:])

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(6)
    shape = LIN2
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        X = C'*LIN2.func(R)

        l = (C[2,:]-C[1,:])./L
        n = [-l[2], l[1]]

        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)

            if key == :tx
                tl = vip*l[1]
                tn = vip*n[1]
                # tl = vip*dot([1,0], l)
                # tn = vip*dot([1,0], n)
            elseif key == :ty
                tl = vip*l[2]
                tn = vip*n[2]
                # tl = vip*dot([0,1], l)
                # tn = vip*dot([0,1], n)
            elseif key == :tl
                tl = vip
                tn = 0.0
            elseif key == :tn
                tl = 0.0
                tn = vip
            end
        else
            error("This beam element is for 2D only")
        end

        N = beam_shape_func(R[1]*L/2+L/2, L)
        Nl = [ N[1], 0, 0, N[4], 0, 0 ]
        Nn = [ 0, N[2], N[3], 0, N[5], N[6] ]
        
        F += (Nl*tl + Nn*tn)*L/2*w # F is a vector
    end

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L
    T = [ c s 0  0 0 0
          -s c 0  0 0 0
          0 0 1  0 0 0
          0 0 0  c s 0
          0 0 0 -s c 0
          0 0 0  0 0 1 ]
            
    F = T'*F

    # generate a map
    keys = [:ux, :uy, :rz]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return F, map
end


function elem_stiffness(elem::MechClassicalBeam)
    C  = getcoords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    L3 = L*L*L
    mat = elem.mat
    EA = mat.E*mat.A
    EI = mat.E*mat.I

    K0 = [ EA/L     0         0         -EA/L    0         0
           0       12*EI/L3   6*EI/L2    0     -12*EI/L3   6*EI/L2
           0        6*EI/L2   4*EI/L     0      -6*EI/L2   2*EI/L
          -EA/L     0          0         EA/L     0        0
           0      -12*EI/L3  -6*EI/L2    0      12*EI/L3  -6*EI/L2
           0        6*EI/L2   2*EI/L     0      -6*EI/L2   4*EI/L  ]

    display(K0)

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L

    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]

    map = elem_map(elem)

    return T'*K0*T, map, map
end

function elem_mass(elem::MechClassicalBeam)
    C  = getcoords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    mat = elem.mat
    EA = mat.E*mat.A
    EI = mat.E*mat.I


    M0 = mat.ρ*L/420.0*[ 140   0      0      70    0      0
                         0     156    22*L   0     54    -13*L
                         0     22*L   4*L2   0     13*L  -3*L2
                         70    0      0      140   0      0
                         0     54     13*L   0     156   -22*L
                         0    -13*L  -3*L2   0    -22*L   4*L2 ]

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L
    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]

    map = elem_map(elem)
    return T'*M0*T, map, map
end


function elem_update!(elem::MechClassicalBeam, U::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end

