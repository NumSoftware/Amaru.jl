# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru
mutable struct HydroJoint<:Hydromech
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    active::Bool
    linked_elems::Array{Element,1}
    env ::ModelEnv

    function HydroJoint()
        return new()
    end
end

# Return the shape family that works with this element
matching_shape_family(::Type{HydroJoint}) = JOINTCELL


function elem_config_dofs(elem::HydroJoint)
    nnodes = length(elem.nodes)
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :uw, :fw)
    end
end


function elem_init(elem::HydroJoint)
    # Get linked elements
    e1 = elem.linked_elems[1]
    e2 = elem.linked_elems[2]

    # Volume from first linked element
    V1 = 0.0
    C1 = getcoords(e1)
    for ip in e1.ips
        dNdR = e1.shape.deriv(ip.R)
        J    = C1'*dNdR
        detJ = det(J)
        V1  += detJ*ip.w
    end

    # Volume from second linked element
    V2 = 0.0
    C2 = getcoords(e2)
    for ip in e2.ips
        dNdR = e2.shape.deriv(ip.R)
        J    = C2'*dNdR
        detJ = det(J)
        V2  += detJ*ip.w
    end

    # Area of joint element
    A = 0.0
    C = getcoords(elem)
    n = div(length(elem.nodes), 3)
    C = C[1:n, :]
    fshape = elem.shape.facet_shape

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)
        J    = C'*dNdR
        detJ = norm2(J)
        A += detJ*ip.w
    end

    # Calculate and save h at joint element's integration points
    h = (V1+V2)/(2.0*A)
    for ip in elem.ips
        ip.state.h = h
    end
end


function elem_conductivity_matrix(elem::HydroJoint)
    ndim     = elem.env.ndim
    th       = elem.env.ana.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape

    C        = getcoords(elem)[1:nbsnodes,:]
    Cl       = zeros(nbsnodes, ndim-1)
    J        = Array{Float64}(undef, ndim-1, ndim)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)

    H        = zeros(3*nbsnodes, 3*nbsnodes)

    nodes_p  = []
    
    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Bp matrix  
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # new coordinate nodes

        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 

        # compute NN matrix
        Np = elem.shape.basic_shape.func(ip.R)  
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']

        # compute H
        coef  = detJ*ip.w*th*elem.mat.kt
        H -= coef*Nb'*Nb
        H -= coef*Nt'*Nt

         # compute crack aperture
        coef = detJ*ip.w*th*(elem.mat.w^3)/(12*elem.mat.η) 
        H -= coef*Bf'*Bf
    end
    
    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return H, map, map, nodes_p
end


function elem_compressibility_matrix(elem::HydroJoint)
    ndim     = elem.env.ndim
    th       = elem.env.ana.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = getcoords(elem)[1:nbsnodes,:]

    J   = Array{Float64}(undef, ndim-1, ndim)
    Cpp = zeros(3*nbsnodes, 3*nbsnodes)

    nodes_p  = []
    
    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 

    for ip in elem.ips
        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J = C'*dR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np matrix
        Np = elem.shape.basic_shape.func(ip.R)  
        N0 = 0*Np
        Nf = [N0' N0' Np']

        # compute Cpp
        coef = detJ*ip.w*elem.mat.β*elem.mat.w*th
        Cpp -= coef*Nf'*Nf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return Cpp, map, map
end


function elem_RHS_vector(elem::HydroJoint)
    ndim     = elem.env.ndim
    th       = elem.env.ana.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = getcoords(elem)[1:nbsnodes,:]

    J        = Array{Float64}(undef, ndim-1, ndim)
    Cl       = zeros(nbsnodes, ndim-1)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)
    Z        = zeros(ndim) 
    Z[end]   = 1.0
    bf       = zeros(ndim-1) 
    Q        = zeros(3*nbsnodes)

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 
 
    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) #rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  #coordinate of new nodes

        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 
        
        # compute Q

        # compute crack aperture
        coef = detJ*ip.w*th*(elem.mat.w^3)/(12*elem.mat.η)   
        bf = T[(2:end), (1:end)]*Z*elem.env.ana.γw
        
        @gemm Q += coef*Bf'*bf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return Q, map
end
#=
function elem_internal_forces(elem::HydroJoint, F::Array{Float64,1})
    ndim     = elem.env.ndim
    th       = elem.env.ana.thickness
    nnodes   = length(elem.nodes)
    bsnodes  = elem.shape.basic_shape.npoints
    nlnodes  = div(nnodes, 3)
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape

    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dFw    = zeros(nnodes)
    Bp     = zeros(ndim-1, nlnodes)

    J      = Array{Float64}(undef, ndim-1, ndim)

    C      = getcoords(elem)[1:nlnodes,:]
    Cl     = zeros(nlnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    Bpuwf  = zeros(ndim-1)
    Z      = zeros(ndim)
    Z[end] = 1.0


    for ip in elem.ips
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = C'*dNdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np =  N'

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes
        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR #dNdX
        B0 = 0*Bp
        Bf = [B0 B0 Bp]

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.env.ana.γw

        # compute Bu matrix
        Np = elem.shape.basic_shape.func(ip.R)
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']
        Nf = [N0' N0'  Np']

        # internal volumes dFw
        uwf  = ip.state.uw[3]
        coef = detJ*ip.w*elem.mat.β*th
        dFw -= coef*Nf'*uwf

        # longitudinal flow
        coef = detJ*ip.w*th  
        S = ip.state.S
        dFw -= coef*Bf'*S

        # transverse flow
        D = ip.state.D
        dFw -= coef*Nt'*D[1]
        dFw -= coef*Nb'*D[2]
    end

    F[map_w] += dFw
end
=#

function update_elem!(elem::HydroJoint, U::Array{Float64,1}, Δt::Float64)
    ndim     = elem.env.ndim
    th       = elem.env.ana.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.basic_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape


    map_w  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dUw    = U[map_w] # nodal pore-pressure increments
    Uw     = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ] # nodal pore-pressure at step n
    Uw    += dUw # nodal pore-pressure at step n+1

    dFw    = zeros(nnodes)

    J      = Array{Float64}(undef, ndim-1, ndim)
    Δω     = zeros(ndim)
    Δuw    = zeros(3)
    Bu     = zeros(ndim, dnlnodes*ndim)
    C      = getcoords(elem)[1:nbsnodes,:]
    Cl     = zeros(nlnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    Bp     = zeros(ndim-1, nbsnodes)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    BpUwf  = zeros(ndim-1)
    Z      = zeros(ndim)
    Z[end] = 1.0

    for ip in elem.ips

        # compute shape Jacobian
        dNdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = C'*dNdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np = elem.shape.basic_shape.func(ip.R)
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']
        Nf = [N0' N0'  Np']

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes

        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR #dNdX
        B0 = 0*Bp
        Bf = [B0 B0 Bp]

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.env.ana.γw

        # interpolation to the integ. point
        Δuw  = [Np'*dUw[1:nbsnodes]; Np'*dUw[nbsnodes+1:2*nbsnodes]; Np'*dUw[2*nbsnodes+1:end]]
        G    = [ dot(Nt,Uw); dot(Nb,Uw)]
        BfUw = Bf*Uw + bf

        Vt, L = update_state!(elem.mat, ip.state, Δuw, G, BfUw, Δt)

        # internal volumes dFw

        # compute fluid compressibility
        coef = detJ*ip.w*elem.mat.β*elem.mat.w*th
        dFw -= coef*Nf'*Δuw[3]

        # longitudinal flow
        coef = Δt*detJ*ip.w*th
        dFw -= coef*Bf'*L

        # transverse flow
        dFw += coef*Nt'*Vt[1]  
        dFw += coef*Nb'*Vt[2] 
    end

    F[map_w] .+= dFw
    return success()
end


function elem_extrapolated_node_vals(elem::HydroJoint)
    nips = length(elem.ips)

    E  = extrapolator(elem.shape.facet_shape, nips)
    Uwn = E*[ ip.state.uw[3] for ip in elem.ips ]

    node_vals = OrderedDict{Symbol, Array{Float64,1}}()
    node_vals[:uwn] = [ Uwn; Uwn; Uwn]

    return node_vals
end