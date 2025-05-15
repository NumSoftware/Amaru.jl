# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechCohesive

MechCohesive_params = [
    FunInfo(:MechCohesive, "An isoparametric cohesive/cohesive element"),
]
@doc docstring(MechCohesive_params) MechCohesive
struct MechCohesiveProps<:ElemProperties
    function MechCohesiveProps(; props...)
        return new()
    end    
end

mutable struct MechCohesive<:Mech
    id    ::Int
    ctx::Context
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::MechCohesiveProps
    active::Bool
    linked_elems::Array{Element,1}

    function MechCohesive()
        return new()
    end
end

# Return the shape family that works with this element
compat_shape_family(::Type{MechCohesive}) = JOINTCELL
compat_elem_props(::Type{MechCohesive}) = MechCohesiveProps


function elem_init(elem::MechCohesive)

    # Avg volume of linked elements
    V = 0.0
    for elem in elem.linked_elems
        V += cell_extent(elem)
    end
    V /= length(elem.linked_elems)

    # Area of cohesive element
    A = 0.0
    C = getcoords(elem)
    n = div(length(elem.nodes), 2)
    C = C[1:n, :]
    fshape = elem.shape.facet_shape
    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)
        J    = C'*dNdR
        detJ = norm2(J)
        detJ <= 0 && error("Invalid Jacobian norm for cohesive element")
        A += detJ*ip.w
    end

    # Calculate and save h at cohesive element's integration points
    h = V/A
    for ip in elem.ips
        ip.state.h = h
    end

end

# function matrixT(J::Matrix{Float64})
#     if size(J,2)==2
#         L2 = vec(J[:,1])
#         L3 = vec(J[:,2])
#         L1 = cross(L2, L3)  # L1 is normal to the first cohesive face
#         L2 = cross(L3, L1)
#         normalize!(L1)
#         normalize!(L2)
#         normalize!(L3)
#         return collect([L1 L2 L3]') # collect is used to avoid Adjoint type
#     else
#         L2 = normalize(vec(J))
#         L1 = [ L2[2], -L2[1] ] # It follows the anti-clockwise numbering of 2D elements: L1 should be normal to the first cohesive face
#         return collect([L1 L2]')
#     end
# end


function elem_stiffness(elem::MechCohesive)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    fshape = elem.shape.facet_shape

    C = getcoords(elem)[1:hnodes,:]
    B = zeros(ndim, nnodes*ndim)
    K = zeros(nnodes*ndim, nnodes*ndim)

    DB = zeros(ndim, nnodes*ndim)
    J  = zeros(ndim, ndim-1)
    NN = zeros(ndim, nnodes*ndim)

    for ip in elem.ips
    	if elem.ctx.stressmodel==:axisymmetric
            th = 2*pi*ip.coord.x
        end
        
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        T   = matrixT(J)
        NN .= 0.0  # NN = [ -N[]  N[] ]
        for i in 1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @mul B = T*NN

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.mat, ip.state)
        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function update_elem!(elem::MechCohesive, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    fshape = elem.shape.facet_shape
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    C = getcoords(elem)[1:hnodes,:]
    B = zeros(ndim, nnodes*ndim)

    J  = zeros(ndim, ndim-1)
    NN = zeros(ndim, nnodes*ndim)
    Δω = zeros(ndim)

    for ip in elem.ips
    	if elem.ctx.stressmodel==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        T = matrixT(J)
        NN[:,:] .= 0.0  # NN = [ -N[]  N[] ]
        for i in 1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        @mul B = T*NN

        # internal force
        @mul Δω = B*dU
        Δσ, status = update_state!(elem.mat, ip.state, Δω)
        failed(status) && return dF, map, status
        coef = detJ*ip.w*th
        @mul dF += coef*B'*Δσ
    end

    return dF, map, success()
end


function elem_recover_nodal_values(elem::MechCohesive)
    nips = length(elem.ips)

    keys = output_keys(elem.mat)
    vals = zeros(nips, length(keys))
    for (i,ip) in enumerate(elem.ips)
        dict = ip_state_vals(elem.mat, ip.state)
        vals[i,:] = [ dict[key] for key in keys ]
    end
    
    node_vals = OrderedDict{Symbol, Array{Float64,1}}()
    E = extrapolator(elem.shape.facet_shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V ]
    end

    return node_vals
end
