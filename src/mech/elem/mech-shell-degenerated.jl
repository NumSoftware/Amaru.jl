# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    ShellDegenerated
A bulk finite element for mechanical equilibrium analyses.
"""
mutable struct ShellDegenerated<:Mechanical
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv

    function ShellDegenerated();
        return new()
    end
end

matching_shape_family(::Type{ShellDegenerated}) = SOLID_CELL

function elem_init(elem::ShellDegenerated)
    ipdata_ty = typeof(elem.ips[1].state)
    if :h in fieldnames(ipdata_ty)
        # Element volume/area
        V = 0.0
        C = getcoords(elem)
        for ip in elem.ips
            dNdR = elem.shape.deriv(ip.R)
            J    = dNdR*C
            detJ = det(J)
            @assert detJ>0
            V   += detJ*ip.w
        end

        # Representative length size for the element
        nips = length(elem.ips)
        ndim = elem.env.ndim
        h = V^(1/ndim)

        for ip in elem.ips
            ip.state.h = h
        end
    end

    return nothing
end


function distributed_bc(elem::ShellDegenerated, facet::Union{Facet, Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.thickness
    suitable_keys = (:tx, :ty, :tz, :tn)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet!==nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = getcoords(nodes, ndim)

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
        J = C'*D
        #J = D*C
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
            if elem.env.modeltype=="axisymmetric"
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            if key == :tx
                Q = [vip, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, vip, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, vip]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = vip*normalize(n)
            end
        end
        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end



# Rotation Matrix
function Rot_Matrix(elem::ShellDegenerated, J::Matx, T::Matx)

    t = elem.mat.t
    # J (3x2)
    V1 = vec(J[:,1])
    V2 = vec(J[:,2])

    V3 = t/2*cross(V1, V2)  # multiplicar ou não por t/2?

    V = [ normalize(V1)
          normalize(V2)
          normalize(V3) ]

          l1, m1, n1 = normalize(V1)
          l2, m2, n2 = normalize(V2)
          l3, m3, n3 = normalize(V3)

#=
    l1, m1, n1 = V[1,:]
    l2, m2, n2 = V[2,:]
    l3, m3, n3 = V[3,:]
=#


    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;   T[1,4] =   SR2*m1*n1;  T[1,5] =   SR2*n1*l1;  T[1,6] =   SR2*l1*m1;
    T[2,1] =     l2*l2;  T[2,2] =     m2*m2;  T[2,3] =     n2*n2;   T[2,4] =   SR2*m2*n2;  T[2,5] =   SR2*n2*l2;  T[2,6] =   SR2*l2*m2;
    T[3,1] =     l3*l3;  T[3,2] =     m3*m3;  T[3,3] =     n3*n3;   T[3,4] =   SR2*m3*n3;  T[3,5] =   SR2*n3*l3;  T[3,6] =   SR2*l3*m3;
    T[4,1] = SR2*l2*l3;  T[4,2] = SR2*m2*m3;  T[4,3] = SR2*n2*n3;   T[4,4] = m2*n3+m3*n2;  T[4,5] = l2*n3+l3*n2;  T[4,6] = l2*m3+l3*m2;
    T[5,1] = SR2*l3*l1;  T[5,2] = SR2*m3*m1;  T[5,3] = SR2*n3*n1;   T[5,4] = m3*n1+m1*n3;  T[5,5] = l3*n1+l1*n3;  T[5,6] = l3*m1+l1*m3;
    T[6,1] = SR2*l1*l2;  T[6,2] = SR2*m1*m2;  T[6,3] = SR2*n1*n2;   T[6,4] = m1*n2+m2*n1;  T[6,5] = l1*n2+l2*n1;  T[6,6] = l1*m2+l2*m1;
end


function setB(elem::ShellDegenerated, J::Matrix{Float64}, ip::Ip, dNdX::Matx, N::Vect, B::Matx)
    nnodes, ndim = size(dNdX)
    B .= 0.0
    t = elem.mat.t
    ζ = ip.R[3] #não sei se está certo
  
    C = getcoords(elem)

    # J (3x2)
   
    t = elem.mat.t
    # J (3x2)
    V1 = vec(J[:,1])
    V2 = vec(J[:,2])

    l1, m1, n1 = normalize(V1)
    l2, m2, n2 = normalize(V2)


        for i in 1:nnodes

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        dNdz = dNdX[i,3]
        
        j    = i-1

        B[1,1+j*ndim] = dNdx;                                               B[1,4+j*ndim] = -ζ*dNdx*t/2*l2;              B[1,5+j*ndim] = -ζ*dNdx*t/2*l1

                              B[2,2+j*ndim] = dNdy;                         B[2,4+j*ndim] = -ζ*dNdy*t/2*m2;              B[2,5+j*ndim] = -ζ*dNdy*t/2*m1

                                                     B[3,3+j*ndim] = dNdz;  B[3,4+j*ndim] = -ζ*dNdz*t/2*n2;              B[3,5+j*ndim] = -ζ*dNdz*t/2*n1

        B[4,1+j*ndim] = dNdx;  B[4,2+j*ndim] = dNdy;                        B[4,4+j*ndim] = -ζ*t/2*(dNdy*l2+dNdx*m2);    B[4,5+j*ndim] = -ζ*t/2*(dNdy*l1+dNdx*m1)

                               B[5,2+j*ndim] = dNdy; B[5,3+j*ndim] = dNdx;  B[5,4+j*ndim] = -ζ*t/2*(dNdz*m2+dNdy*n2);    B[5,5+j*ndim] = -ζ*t/2*(dNdz*m1+dNdy*n1)

        B[6,1+j*ndim] = dNdx;                        B[6,3+j*ndim] = dNdz;  B[6,4+j*ndim] = -ζ*t/2*(dNdz*l2+dNdx*n2);    B[6,5+j*ndim] = -ζ*t/2*(dNdz*l1+dNdx*n1)
 
 

    end

end


function Dmatrix(elem::ShellDegenerated)

    nu = elem.mat.nu
    E1 = elem.mat.E/(1-elem.mat.nu^2)
    G  = elem.mat.E/(2*(1+elem.mat.nu))
    G1 = 5/6*G

        D =   [E1   nu*E1 0 0 0  0
             nu*E1    E1  0 0 0  0
              0         0 0 0 0  0
              0         0 0 G 0  0
              0         0 0 0 G1 0
              0         0 0 0  0 G1 ]
#=
              D =   [E1   nu*E1 0 0  0
                    nu*E1    E1 0 0  0
                     0          0 G 0  0
                      0         0 0 G1 0
                      0         0 0  0 G1 ]
=#
    return D
end

function elem_config_dofs(elem::ShellDegenerated)
    ndim = elem.env.ndim
    ndim == 1 && error("ShellDegenerated: Shell elements do not work in 1d analyses")
    #if ndim==2
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
    #else
        #error("ShellDegenerated: Shell elements do not work in this analyses")
        #=
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
        =#
    #end
end


function elem_map(elem::ShellDegenerated)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) 
    #end

    #dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    #dof_keys = (:ux, :uy, :uz)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function elem_stiffness(elem::ShellDegenerated)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    #@show nnodes

    C = getcoords(elem)
    K = zeros(6*nnodes, 6*nnodes)
    B = zeros(6, 5*nnodes)
    JJ  = zeros(9,9)
    Rot_K  = zeros(5*nnodes,6*nnodes)
    Klocal = zeros(5*nnodes, 5*nnodes)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    V  = Array{Float64}(undef, ndim, ndim)
    T = zeros(6,6)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    D = Dmatrix(elem)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute B matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        #@gemm J = C'*dNdR
        #@gemm dNdX = dNdR*inv(J)
        J = C'*dNdR
        #dNdX = dNdR*inv(J)
        dNdX = dNdR*pinv(J) #! correto
        Rot_Matrix(elem, J, T)
        #@show T
        #@show V
        
        # aux = (elem.mat.t/2)*cross(J[:,1],J[:,2])
        #@show aux

        # normalize!(aux)
        #@show aux
        #falta normalizar o vetor     


        detJ = norm2(J)
        #detJ = det(J)
        #@show detJ

        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

           
        setB(elem, J, ip, dNdX, N, B)
        #@show B

       # coef = detJ*ip.w*th
        coef = detJ*ip.w


        Klocal += (B'*T'*D*T*B)*coef
        @show Klocal
    end

    for i in 1:1
        Rot_K[(i-1)*5+1:i*5, (i-1)*6+1:i*6] = [T[1:2,:]
                                               T[4:6,:] ]
    end

    K = Rot_K'*Klocal*Rot_K

    #keys = (:ux, :uy, :uz)[1:ndim]
    #keys =(:ux, :uy, :uz, :rx, :ry, :rz)[1:ndim]
    #map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    map = elem_map(elem)
    return K, map, map
end


function elem_update!(elem::ShellDegenerated, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
    return success()
end


#=
function elem_update!(elem::ShellDegenerated, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    dU = U[map]
    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)
    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    Δε = zeros(6)
    C = getcoords(elem)
    for ip in elem.ips
        if elem.env.modeltype=="axisymmetric"
            th = 2*pi*ip.coord.x
        end
        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        setB(elem, ip, dNdX, B)
        @gemv Δε = B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return failure("ShellDegenerated: Error at integration point $(ip.id)")
        #if failed(status)
            #status.message = "ShellDegenerated: Error at integration point $(ip.id)\n" * status.message
            #return status
        #end
        coef = detJ*ip.w*th
        @gemv dF += coef*B'*Δσ
    end
    F[map] += dF
    return success()
end
function elem_vals(elem::ShellDegenerated)
    vals = OrderedDict{Symbol,Float64}()
    if haskey(ip_state_vals(elem.mat, elem.ips[1].state), :damt)
        mean_dt = mean( ip_state_vals(elem.mat, ip.state)[:damt] for ip in elem.ips )
        vals[:damt] = mean_dt
        mean_dc = mean( ip_state_vals(elem.mat, ip.state)[:damc] for ip in elem.ips )
        vals[:damc] = mean_dc
    end
    #vals = OrderedDict{String, Float64}()
    #keys = elem_vals_keys(elem)
#
    #dicts = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
    #nips = length(elem.ips)
#
    #for key in keys
        #s = 0.0
        #for dict in dicts
            #s += dict[key]
        #end
        #vals[key] = s/nips
    #end
    return vals
end
=#
