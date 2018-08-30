# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MechRod<:Mechanical
    id    ::Int
    shape ::ShapeType
    cell  ::Cell
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::TagType
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    shared_data ::SharedAnalysisData

    function MechRod()
        return new() 
    end
end

matching_shape_family(::Type{MechRod}) = LINE_SHAPE

function elem_stiffness(elem::MechRod)
    ndim   = elem.shared_data.ndim
    nnodes = length(elem.nodes)
    A = elem.mat.A
    C = elem_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, 1, ndim)

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

        E    = calcD(elem.mat,ip.data)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end
            

function elem_mass(elem::MechRod) 
                
    ndim   = elem.shared_data.ndim
    nnodes = length(elem.nodes)
    ρ = elem.mat.ρ
    A = elem.mat.A
    
    C = elem_coords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    J  = Array{Float64}(undef, 1, ndim)
    N = zeros(ndim, ndim*nnodes)
    #N = Array{Float64}(undef, ndim, ndim*nnodes)
    
    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        Ni = elem.shape.func(ip.R) #encontrei em shape.jl FemMesh
        setNt(ndim,Ni,N)    

        @gemm J = dNdR*C
        detJ = norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*A*detJ*ip.w
        @gemm M += coef*N'*N 

    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return M, map, map
end
                        
function setNt(ndim::Int,Ni::Vect, N::Matx)
    nnodes = length(Ni) 
    N .= 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
        end
    elseif ndim==3
        for i in 1:nnodes
            j    = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
            N[3,3+j*ndim] = Ni[i]
       end
    else
        for i in 1:nodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
        end    
    end
    
end

function distributed_bc(elem::MechRod, facet::Union{Facet, Nothing}, key::Symbol, fun::Functor)
    ndim  = elem.shared_data.ndim

    # Check bcs
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")
    !(key in (:tx, :ty, :tz, :tn)) && error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.shared_data.t
    A      = elem.mat.A

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = nodes_coords(nodes, ndim)

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
        nJ = norm2(J)
        X = C'*N
        if ndim==2
            x, y = X
            val = fun(t,x,y,0.0)
            if key == :tx
                Q = [val, 0.0]
            elseif key == :ty
                Q = [0.0, val]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = val*n/norm(n)
            end
        else
            x, y, z = X
            val = fun(t,x,y,z)
            if key == :tx
                Q = [val, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, val, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, val]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = val*n/norm(n)
            end
        end
        F += N*Q'*(A*nJ*w) # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end
                        
                        
function elem_update!(elem::MechRod, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim   = elem.shared_data.ndim
    nnodes = length(elem.nodes)
    A      = elem.mat.A
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    C  = elem_coords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array{Float64}(undef, 1, ndim)
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

        deps = (B*dU)[1]
        dsig = stress_update(elem.mat, ip.data, deps)
        coef = A*detJ*ip.w
        dF  += coef*B'*dsig
    end

    F[map] += dF
end

function elem_vals(elem::MechRod)
    #get area and average axial force
    vals = OrderedDict(:A => elem.mat.A )
    mean_sa = mean( ip_state_vals(elem.mat, ip.data)[:sa] for ip in elem.ips )
    vals[:sa] = mean_sa
    vals[:fa] = elem.mat.A*mean_sa
    return vals
end
