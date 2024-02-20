# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Distributed natural boundary conditions for faces and edges of bulk elements
function acoustic_mech_solid_bc(elem::AcousticMech, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.ana.thickness
    suitable_keys = (:tq,) # tq: mass acceleration per area?, ax: x acceleration

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a AcousticFluid element")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = getcoords(nodes, ndim)
    

    # Calculate the nodal values
    F     = zeros(nnodes)
    shape = target.shape
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = ips[i].coord
        w = ips[i].w
        N = shape.func(R)
        D = shape.deriv(R)

        J = C'*D
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if elem.env.ana.stressmodel=="axisymmetric"
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
        end

        if  ndim==2
            n = [J[2], -J[1]]
        else
            n = cross(J[:,1], J[:,2])
        end
        normalize!(n)

        coef = vip*norm(J)*w*th
        F .+= coef*N # F is a vector
    end

    # generate a map
    map  = [ node.dofdict[:up].eq_id for node in target.nodes ]

    return F, map
end