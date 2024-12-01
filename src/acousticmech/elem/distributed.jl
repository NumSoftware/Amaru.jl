# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Distributed natural boundary conditions for faces and edges of bulk elements
function acoustic_mech_bc(elem::Element, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.ctx.ndim
    th    = elem.ctx.thickness
    suitable_keys = (:tq,) # tq: mass acceleration per area?, ax: x acceleration

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable in a AcousticFluid element")

    nodes  = facet.nodes
    nnodes = length(nodes)

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the facet coordinates matrix
    C = getcoords(nodes, ndim)

    # Calculate the nodal values
    F     = zeros(nnodes)
    shape = facet.shape
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
            vip = evaluate(val, t=t, x=x, y=y)
            if elem.ctx.stressmodel==:axisymmetric
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = evaluate(val, t=t, x=x, y=y, z=z)
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
    map  = [ node.dofdict[:up].eq_id for node in facet.nodes ]

    return F, map
end