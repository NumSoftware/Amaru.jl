# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function mech_line_distributed_forces(elem::Element, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim = elem.env.ndim
    suitable_keys = (:qx, :qy, :qz, :qn, :wx, :wy, :wz)
    isedgebc = key in (:qx, :qy, :qz, :qn) 
    
    # Check keys
    key in suitable_keys || error("mech_line_distributed_forces: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key in (:wz,:qz) && ndim==2) && error("mech_line_distributed_forces: boundary condition $key is not applicable in a 2D analysis")
    (key == :qn && ndim==3) && error("mech_line_distributed_forces: boundary condition $key is not applicable in a 3D analysis")

    nodes  = elem.nodes
    nnodes = length(nodes)
    t      = elem.env.t
    A      = isedgebc ? 1.0 : elem.mat.A

    # Calculate the elem coordinates matrix
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = elem.shape
    ips   = get_ip_coords(shape)

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
            Q = zeros(2)
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            Q = zeros(3)
        end

        if key == :qx
            Q[1] = vip
        elseif key == :qy
            Q[2] = vip
        elseif key == :qz
            Q[3] = vip
        elseif key == :qn
            if  ndim==2
                n = [J[1,2], -J[1,1]]
            else
                n = cross(J[1,:], J[2,:])
            end
            Q = vip*normalize(n)
        end

        F += N*Q'*(A*nJ*w) # F is a matrix
    end

    # generate a map
    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


# Distributed natural boundary conditions for faces and edges
function mech_solid_boundary_forces(elem::Element, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    suitable_keys = (:qx, :qy, :qz, :qn, :tx, :ty, :tz, :tn)
    isedgebc = key in (:qx, :qy, :qz, :qn) 
    
    # Check keys
    key in suitable_keys || error("mech_solid_boundary_forces: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    key in (:tz,:qz) && ndim==2 && error("mech_solid_boundary_forces: boundary condition $key is not applicable in a 2D analysis")
    isedgebc && elem.env.modeltype=="axisymmetric" && error("mech_solid_boundary_forces: boundary condition $key is not applicable in a axisymmetric analysis")
    isedgebc && facet.shape.ndim==2 && error("mech_solid_boundary_forces: boundary condition $key is not applicable on surfaces")
    !isedgebc && facet.shape.ndim==1 && ndim==3 && error("mech_solid_boundary_forces: boundary condition $key is not applicable on 3D edges")

    th     = isedgebc ? 1.0 : elem.env.thickness
    nodes  = facet.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Calculate the facet coordinates 
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = facet.shape
    ips   = get_ip_coords(shape)

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
            Q = zeros(2)
            elem.env.modeltype=="axisymmetric" && (th = 2*pi*X[1])
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            Q = zeros(3)
        end

        if key in (:tx, :qx)
            Q[1] = vip
        elseif key in (:ty, :qy)
            Q[2] = vip
        elseif key in (:tz, :qz)
            Q[3] = vip
        elseif key in (:tn , :qn)
            if  ndim==2
                n = [J[1,2], -J[1,1]]
            else
                n = cross(J[1,:], J[2,:])
            end
            Q = vip*normalize(n)
        end

        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in facet.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


function mech_solid_body_forces(elem::Element, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.thickness
    suitable_keys = (:wx, :wy, :wz)

    # Check keys
    key in suitable_keys || error("mech_solid_body_forces: condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key == :wx && ndim==2) && error("mech_solid_body_forces: key $key is not applicable in a 2D analysis")

    nodes  = elem.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Calculate the elem coordinates matrix
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = elem.shape
    ips   = get_ip_coords(shape)

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
            Q = zeros(2)
            elem.env.modeltype=="axisymmetric" && (th = 2*pi*X[1])
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            Q = zeros(3)
        end

        if key == :wx
            Q[1] = vip
        elseif key == :wy
            Q[2] = vip
        elseif key == :wz
            Q[3] = vip
        end

        coef = det(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end

function mech_shell_boundary_forces(elem::Element, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.mat.th
    suitable_keys = (:tx, :ty, :tz, :tn)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    nodes  = facet.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the facet coordinates matrix
    C = getcoords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = facet.shape
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
        X = C'*N

        x, y, z = X
        vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        if key == :tx
            Q = [vip, 0.0, 0.0]
        elseif key == :ty
            Q = [0.0, vip, 0.0]
        elseif key == :tz
            Q = [0.0, 0.0, vip]
        elseif key == :tn
            n = cross(J[:,1], J[:,2])
            Q = vip*normalize(n)
        end

        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in facet.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end


function mech_shell_body_forces(elem::Element, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.mat.th
    suitable_keys = (:wx, :wy, :wz, :wn)

    # Check keys
    key in suitable_keys || error("mech_shell_body_forces: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")

    nodes  = elem.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    C = getcoords(nodes, ndim)
    Q = zeros(ndim) # Vector with values to apply
    F     = zeros(nnodes, ndim) # Nodal values
    shape = elem.shape
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = C'*D
        X = C'*N

        x, y, z = X
        vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        if key == :wx
            Q = [vip, 0.0, 0.0]
        elseif key == :wy
            Q = [0.0, vip, 0.0]
        elseif key == :wz
            Q = [0.0, 0.0, vip]
        elseif key == :wn
            n = cross(J[:,1], J[:,2])
            Q = vip*normalize(n)
        end

        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end
