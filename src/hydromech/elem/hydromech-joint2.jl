# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru
mutable struct HydroMechJoint2<:Hydromechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env ::ModelEnv

    function HydroMechJoint2()
        return new()
    end
end

# Return the shape family that works with this element
matching_shape_family(::Type{HydroMechJoint2}) = JOINT_SHAPE

function elem_config_dofs(elem::HydroMechJoint2)
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            add_dof(node, :uw, :fw)
        end
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)
    end
end

function elem_init(elem::HydroMechJoint2)
    # Get linked elements
    e1 = elem.linked_elems[1]
    e2 = elem.linked_elems[2]

    # Volume from first linked element
    V1 = 0.0
    C1 = elem_coords(e1)

    for ip in e1.ips
        dNdR = e1.shape.deriv(ip.R)
        J    = dNdR*C1
        detJ = det(J)
        V1  += detJ*ip.w
    end

    # Volume from second linked element
    V2 = 0.0
    C2 = elem_coords(e2)
    for ip in e2.ips
        dNdR = e2.shape.deriv(ip.R)
        J    = dNdR*C2
        detJ = det(J)
        V2  += detJ*ip.w
    end

    # Area of joint element
    A = 0.0
    C = elem_coords(elem)
    n = div(length(elem.nodes), 2)
    C = C[1:n, :]
    fshape = elem.shape.facet_shape
    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)
        J    = dNdR*C
        detJ = norm2(J)
        A += detJ*ip.w
    end

    # Calculate and save h at joint element's integration points
    h = (V1+V2)/(2.0*A)
    for ip in elem.ips
        ip.data.h = h
    end
end

function elem_stiffness(elem::HydroMechJoint2)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nlnodes,:]

    J        = Array{Float64}(undef, ndim-1, ndim)
    NN       = zeros(ndim, nnodes*ndim)
    Bu       = zeros(ndim, nnodes*ndim)
    DBu      = zeros(ndim, nnodes*ndim)
    K        = zeros(nnodes*ndim, nnodes*ndim)

    for ip in elem.ips
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C
        detJ = norm2(J)

        # compute Bu matrix
        T   = matrixT(J)
        NN .= 0.0  # NN = [ -N[]  N[] ]

        for i=1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] = N[i]
            end
        end

        @gemm Bu = T*NN

        # compute K
        coef = detJ*ip.w*th
        D    = mountD(elem.mat, ip.data)
        @gemm DBu = D*Bu
        @gemm K  += coef*Bu'*DBu
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


function elem_coupling_matrix(elem::HydroMechJoint2)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nlnodes,:]

    J        = Array{Float64}(undef, ndim-1, ndim)
    NN       = zeros(ndim, nnodes*ndim)
    Bu       = zeros(ndim, nnodes*ndim)
    mf       = [1.0, 0.0, 0.0][1:ndim]
    mfNf     = zeros(ndim, 2*nbsnodes)
    Cup      = zeros(nnodes*ndim, 2*nbsnodes) # u-p coupling matrix

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np = elem.shape.basic_shape.func(ip.R)
        Nf = [0.5*Np' 0.5*Np']

        # compute Bu matrix
        N  = fshape.func(ip.R)
        T   = matrixT(J)
        NN .= 0.0

        for i=1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @gemm Bu = T*NN
        # compute Cup
        coef = detJ*ip.w*th
        mfNf = mf*Nf
        Cup -= coef*Bu'*mfNf
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    map_p = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    return Cup, map_u, map_p
end


function elem_conductivity_matrix(elem::HydroMechJoint2)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape

    C        = elem_coords(elem)[1:nbsnodes,:]
    Cl       = zeros(nbsnodes, ndim-1)
    J        = Array{Float64}(undef, ndim-1, ndim)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)

    H        = zeros(2*nbsnodes, 2*nbsnodes)

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = dNpdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # new coordinate nodes

        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR
        Bf = [0.5*Bp 0.5*Bp]

        # compute NN matrix
        Np = elem.shape.basic_shape.func(ip.R)
        Nb = [ Np' -Np']
        Nt = [-Np'  Np']

        # compute H
        coef  = detJ*ip.w*th*elem.mat.kt
        H += coef*Nb'*Nb
        H += coef*Nt'*Nt

         # compute crack aperture
        if elem.mat.kl == 0.0
            if ip.data.upa == 0.0 || ip.data.w[1] <= 0.0
                kl = 0.0
            else
                kl = ip.data.w[1]
            end
        else
            if elem.mat.kl >= ip.data.w[1]
                kl = elem.mat.kl
            else
                kl = ip.data.w[1]
            end
        end

        coef = detJ*ip.w*th*(kl^3)/(12*elem.mat.η)
        H -= coef*Bf'*Bf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return H, map, map, nodes_p
end

#=
function elem_compressibility_matrix(elem::HydroMechJoint2)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nbsnodes,:]

    J   = Array{Float64}(undef, ndim-1, ndim)
    Cpp = zeros(2*nbsnodes, 2*nbsnodes)

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips
        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J = dNpdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np matrix
        Np = elem.shape.basic_shape.func(ip.R)
        Nf = [ 0.5*Np' 0.5*Np']

        # compute Cpp
        coef = detJ*ip.w*elem.mat.β*th
        Cpp -= coef*Nf'*Nf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p ]

    return Cpp, map, map
end
=#

function elem_RHS_vector(elem::HydroMechJoint2)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nbsnodes,:]

    J        = Array{Float64}(undef, ndim-1, ndim)
    Cl       = zeros(nbsnodes, ndim-1)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)
    Z        = zeros(ndim)
    Z[end]   = 1.0
    bf       = zeros(ndim-1)
    Q        = zeros(2*nbsnodes)

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = dNpdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) #rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  #coordinate of new nodes

        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR
        Bf = [0.5*Bp 0.5*Bp]

        # compute Q

        # compute crack aperture
        if elem.mat.kl == 0.0
            if ip.data.upa == 0.0 || ip.data.w[1] <= 0.0
                kl = 0.0
            else
                kl = ip.data.w[1]
            end
        else
            if elem.mat.kl >= ip.data.w[1]
                kl = elem.mat.kl
            else
                kl = ip.data.w[1]
            end
        end

        coef = detJ*ip.w*th*(kl^3)/(12*elem.mat.η)
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw

        @gemm Q += coef*Bf'*bf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return Q, map
end

function elem_internal_forces(elem::HydroMechJoint2, F::Array{Float64,1})
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes[1:nnodes] for key in keys ]
    map_p  = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    dF     = zeros(nnodes*ndim)
    Bu     = zeros(ndim, nnodes*ndim)
    dFw    = zeros(2*nbsnodes)
    Bp     = zeros(ndim-1, nbsnodes)

    J      = Array{Float64}(undef, ndim-1, ndim)
    NN     = zeros(ndim, nnodes*ndim)

    C      = elem_coords(elem)[1:nbsnodes,:]
    Cl     = zeros(nbsnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    Bpuwf  = zeros(ndim-1)
    Z      = zeros(ndim)
    Z[end] = 1.0


    for ip in elem.ips
        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = dNpdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes
        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR #dNdX
        Bf = [0.5*Bp 0.5*Bp]

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw

        # compute Np vector
        Np = elem.shape.basic_shape.func(ip.R)
        Nb = [ Np' -Np']
        Nt = [-Np'  Np']
        Nf = [ 0.5*Np' 0.5*Np']

        # compute NN matrix
        N    = fshape.func(ip.R)
        NN .= 0.0
        for i=1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @gemm Bu = T*NN

        # internal force
        uwf  = (ip.data.uw[1] + ip.data.uw[2])/2
        σ    = ip.data.σ[1:ndim] - mf*uwf # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ

        # internal volumes dFw
        w  = ip.data.w[1:ndim]
        coef = detJ*ip.w*th
        mfw = mf'*w
        dFw-= coef*Nf'*mfw
#=
        coef = detJ*ip.w*elem.mat.β*th
        dFw -= coef*Nf'*uwf
=#
        # longitudinal flow
        coef = detJ*ip.w*th
        S = ip.data.S
        dFw -= coef*Bf'*S

        # transverse flow
        D = ip.data.D
        dFw += coef*Nt'*D[1]
        dFw += coef*Nb'*D[2]
    end

    F[map_u] += dF
    F[map_p] += dFw
end


function elem_update!(elem::HydroMechJoint2, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 2) # half the number of total nodes
    nbsnodes = elem.shape.basic_shape.npoints
    fshape   = elem.shape.facet_shape

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes[1:nnodes] for key in keys ]
    map_p  = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    dU     = U[map_u] # nodal displacement increments
    dUw    = U[map_p] # nodal pore-pressure increments
    Uw     = [ node.dofdict[:uw].vals[:uw] for node in nodes_p ] # nodal pore-pressure at step n
    Uw    += dUw # nodal pore-pressure at step n+1

    dF     = zeros(nnodes*ndim)
    dFw    = zeros(2*nbsnodes)

    J      = Array{Float64}(undef, ndim-1, ndim)
    NN     = zeros(ndim, nnodes*ndim)
    Δω     = zeros(ndim)
    Δuw    = zeros(3)
    Bu     = zeros(ndim, nnodes*ndim)
    C      = elem_coords(elem)[1:nbsnodes,:]
    Cl     = zeros(nbsnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    Bp     = zeros(ndim-1, nbsnodes)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    BpUwf  = zeros(ndim-1)
    Z      = zeros(ndim)
    Z[end] = 1.0

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.basic_shape.deriv(ip.R)

        @gemm J = dNpdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np = elem.shape.basic_shape.func(ip.R)
        Nb = [ Np' -Np']
        Nt = [-Np'  Np']
        Nf = [ 0.5*Np' 0.5*Np']

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes

        @gemm Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR #dNdX
        Bf = [0.5*Bp 0.5*Bp]

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw

        # compute NN matrix
        N    = fshape.func(ip.R)
        NN .= 0.0

        for i=1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @gemm Bu = T*NN

        # interpolation to the integ. point
        uwb  = Np'*dUw[1:nbsnodes]
        uwt  = Np'*dUw[nbsnodes+1:2*nbsnodes]
        uwf  = (uwb + uwt)/2
        Δuw  = [uwb; uwt; uwf]
        G    = [ dot(Nt,Uw); dot(Nb,Uw)]
        BfUw = Bf*Uw + bf

        @gemv Δω = Bu*dU

        # internal force dF
        Δσ, Vt, L = stress_update(elem.mat, ip.data, Δω, Δuw, G, BfUw, Δt)
        Δσ -= mf*uwf # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFw
        coef = detJ*ip.w*th
        mfΔω = mf'*Δω
        dFw -= coef*Nf'*mfΔω
#=
        coef = detJ*ip.w*elem.mat.β*th
        dFw -= coef*Nf'*uwf
=#
        # longitudinal flow
        coef = Δt*detJ*ip.w*th
        dFw -= coef*Bf'*L

        # transverse flow
        dFw -= coef*Nt'*Vt[1]
        dFw -= coef*Nb'*Vt[2]
    end

    F[map_u] .+= dF
    F[map_p] .+= dFw
end
