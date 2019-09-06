# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru
mutable struct HydroJoint<:Hydromechanical
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

    function HydroJoint()
        return new()
    end
end

# Return the shape family that works with this element
matching_shape_family(::Type{HydroJoint}) = JOINT_SHAPE

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
    n = div(length(elem.nodes), 3)
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


function elem_conductivity_matrix(elem::HydroJoint)
    ndim     = elem.env.ndim
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 3) 
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape

    C        = elem_coords(elem)[1:nlnodes,:]
    Cl       = zeros(nlnodes, ndim-1)
    J        = Array{Float64}(undef, ndim-1, ndim)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nlnodes)

    H        = zeros(nnodes, nnodes)

    for ip in elem.ips

        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C 
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Bp matrix  
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # new coordinate nodes

        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 
 
        Np = N'
        N0 = 0*N'
        Nb = [Np N0 -Np]
        Nt = [N0 Np -Np]

        # compute H
        coef  = detJ*ip.w*elem.mat.kt/elem.mat.γw
        H += coef*Nb'*Nb
        H += coef*Nt'*Nt

        if elem.mat.kl == 0.0
            # compute W crack aperture
            if ip.data.t == 0.0 || ip.data.w[1] <= 0.0
                W = 0.0
            else
                W = ip.data.w[1]
            end
            coef = detJ*ip.w*(W^3)/(12*elem.mat.η) 
        else
            coef = detJ*ip.w*elem.mat.kl
        end    

        H -= coef*Bf'*Bf
    end
    
    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return H, map, map
end


function elem_compressibility_matrix(elem::HydroJoint)
    ndim     = elem.env.ndim
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 3) 
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nlnodes,:]

    J   = Array{Float64}(undef, ndim-1, ndim)
    Cpp = zeros(nnodes, nnodes) 

    for ip in elem.ips
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)
        
        @gemm J = dNdR*C
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np matrix
        Np = N'
        N0 = 0*N'
        Nf = [N0 N0 Np]

        # compute Cpp
        coef = detJ*ip.w*elem.mat.β
        Cpp -= coef*Nf'*Nf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Cpp, map, map
end


function elem_RHS_vector(elem::HydroJoint)
    ndim     = elem.env.ndim
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 3) 
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape
    C        = elem_coords(elem)[1:nlnodes,:]

    J        = Array{Float64}(undef, ndim-1, ndim)
    Cl       = zeros(nlnodes, ndim-1)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nlnodes)
    Z        = zeros(ndim) 
    Z[end]   = 1.0
    bf       = zeros(ndim-1) 
    Q        = zeros(nnodes)
 
    for ip in elem.ips

        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C 
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) #rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  #coordinate of new nodes

        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 
        
        # compute Q
        if elem.mat.kl == 0.0
            # compute W crack aperture
            if ip.data.t == 0.0 ||  ip.data.w[1] <= 0.0
                W = 0.0
            else
                W = ip.data.w[1]
            end
            coef = detJ*ip.w*(W^3)/(12*elem.mat.η)
        else
            coef = detJ*ip.w*elem.mat.kl            
        end  
        
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw
        @gemm Q += coef*Bf'*bf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in elem.nodes  ]

    return Q, map
end

function elem_internal_forces(elem::HydroJoint, F::Array{Float64,1})
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 3) 
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape

    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dFw    = zeros(nnodes)
    Bp     = zeros(ndim-1, nlnodes)     

    J      = Array{Float64}(undef, ndim-1, ndim)

    C      = elem_coords(elem)[1:nlnodes,:]
    Cl     = zeros(nlnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    mf     = [1.0, 0.0, 0.0][1:ndim]    
    Bpuwf  = zeros(ndim-1)
    Z      = zeros(ndim) 
    Z[end] = 1.0


    for ip in elem.ips
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C 
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np =  N'

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes
        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR #dNdX

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw
        
        # compute Bu matrix
        N0 = 0*N'
        Nb = [Np N0 -Np]
        Nt = [N0 Np -Np]
        Nf = [N0 N0 Np]

        # internal volumes dFw

        coef = detJ*ip.w*elem.mat.β
        dFw -= coef*Nf'*ip.data.uw[3] 

        # longitudinal flow
        if elem.mat.kl == 0.0
            # compute W crack aperture
            if ip.data.t == 0.0 || ip.data.w[1] <= 0.0
                W = 0.0
            else
                W = ip.data.w[1]
            end
            coef = detJ*ip.w*(W^3)/(12*elem.mat.η)
        else
            coef = detJ*ip.w*elem.mat.kl            
        end  

        #= NÃO RESOLVIDO
        # longitudinal flow
        Bfuw = Bf*uw + bf
        dFw -= coef*Bf'*BfUw 
        =# 

        # transverse flow
        coef  = detJ*ip.w*elem.mat.kt/elem.mat.γw

        Nbuw = ip.data.uw[1] - ip.data.uw[3] 
        dFw += coef*Nb'*Nbuw 

        Ntuw = ip.data.uw[2] - ip.data.uw[3] 
        dFw += coef*Nt'*Ntuw
    end

    F[map_p] += dFw
end


function elem_update!(elem::HydroJoint, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim     = elem.env.ndim
    th       = elem.env.thickness
    nnodes   = length(elem.nodes)
    nlnodes  = div(nnodes, 3) 
    dnlnodes = 2*nlnodes
    fshape   = elem.shape.facet_shape


    map_p  = [ node.dofdict[:uw].eq_id for node in elem.nodes ]

    dUw    = U[map_p] # nodal pore-pressure increments
    Uw     = [ node.dofdict[:uw].vals[:uw] for node in elem.nodes ] # nodal pore-pressure at step n
    Uw    += dUw # nodal pore-pressure at step n+1

    dFw    = zeros(nnodes)  
    dFw2    = zeros(nnodes)    

    J      = Array{Float64}(undef, ndim-1, ndim)
    Δω     = zeros(ndim)
    Δuw    = zeros(3)
    Bu     = zeros(ndim, dnlnodes*ndim)
    C      = elem_coords(elem)[1:nlnodes,:]
    Cl     = zeros(nlnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    Bp     = zeros(ndim-1, nlnodes)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    BpUwf  = zeros(ndim-1)
    Z      = zeros(ndim) 
    Z[end] = 1.0

    for ip in elem.ips

        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @gemm J = dNdR*C 
        detJ = norm2(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np =  N'
        N0 = 0*N'
        Nb = [Np N0 -Np]
        Nt = [N0 Np -Np]
        Nf = [N0 N0 Np]

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes

        @gemm Jl = dNdR*Cl
        Bp = inv(Jl)*dNdR #dNdX
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.mat.γw
        
        # interpolation to the integ. point 
        Δuw  = [Np*dUw[1:nlnodes]; Np*dUw[nlnodes+1:dnlnodes]; Np*dUw[dnlnodes+1:end]] 
        Gt = dot(Nt,Uw)/(elem.mat.γw*1)
        Gb = dot(Nb,Uw)/(elem.mat.γw*1)

        Vt, Vb  = update_state!(elem.mat, ip.data, Δuw, Gt, Gb)

        # internal volumes dFw
        Δuwf = Nf*dUw
        coef = detJ*ip.w*elem.mat.β
        dFw -= coef*Nf'*Δuwf

        # compute kt and kb
        if elem.mat.kl == 0.0
            # compute W crack aperture
            if ip.data.t == 0.0 || ip.data.w[1] <= 0.0
                W = 0.0
            else
                W = ip.data.w[1]
            end
            coef = Δt*detJ*ip.w*(W^3)/(12*elem.mat.η)
        else
            coef = Δt*detJ*ip.w*elem.mat.kl            
        end  

        # longitudinal flow
        BfUw = Bf*Uw + bf
        dFw -= coef*Bf'*BfUw

        # transverse flow
        coef = Δt*detJ*ip.w
        dFw -= coef*Nt'*Vt 
        dFw -= coef*Nb'*Vb 
    end
       

    F[map_p] .+= dFw
end
