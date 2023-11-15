# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechBeam

struct MechBeamProps<:ElemProperties
    ρ::Float64
    γ::Float64
    αs::Float64
    thy::Float64
    thz::Float64

    function MechBeamProps(; args...)
        args = checkargs(args, func_params(MechBeamProps))

        if haskey(args, :A)
            thy = thz = √args.A # assuming circular section
        else
            thy, thz = args.thy, args.thz
        end

        return new(args.rho, args.gamma, args.alpha_s, thy, thz)
    end    
end


func_params(::Type{MechBeamProps}) = [
    FunInfo( :MechBeamProps, "Creates a `MechBeamProps` instance.", ()),
    ArgInfo( :thy, "y' thickness", condition=:(thy>0) ),
    ArgInfo( :thz, "z' thickness", condition=:(thz>0.0) ),
    ArgInfo( :A, "Section area",  condition=:(A>0.0)  ),
    ArgInfo( :gamma, "Specific weight", 0, condition=:(gamma>=0.0) ),
    ArgInfo( :rho, "Density", 0, condition=:(rho>=0.0)  ),
    ArgInfo( :alpha_s, "Shear correction coef.", 5/6, condition=:(alpha_s>0) ),
    ArgOpt( :A, (:thy, :thz) ),
]


"""
    MechBeam
A beam finite element for mechanical equilibrium analyses.
"""
mutable struct MechBeam<:Mech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::MechBeamProps
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv
    Dlmn  ::Array{ Array{Float64,2}, 1}

    function MechBeam()
        return new()
    end
end


compat_shape_family(::Type{MechBeam}) = LINECELL
compat_elem_props(::Type{MechBeam}) = MechBeamProps


function elem_init(elem::MechBeam)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    Dlmn = Array{Float64,2}[]
    C = getcoords(elem)

    for i in 1:nnodes
        Ri = elem.shape.nat_coords[i,:]
        dNdR = elem.shape.deriv(Ri)
        J = C'*dNdR
        L = zeros(ndim,ndim)
        set_rot_x_xp(elem, J, L)
        push!(Dlmn, L)
    end
    elem.Dlmn = Dlmn

    return nothing
end


function setquadrature!(elem::MechBeam, n::Int=0)
    ndim = elem.env.ndim

    if ndim==3
        if n in (0,8)
            nl, nj, nk = 2, 2, 2
        elseif n==12
            nl, nj, nk = 3, 2, 2
        elseif n==2
            nl, nj, nk = 2, 1, 1
        elseif n==3
            nl, nj, nk = 3, 1, 1
        elseif n==32
            nl, nj, nk = 2, 4, 4
        elseif n==48
            nl, nj, nk = 3, 4, 4
        elseif n==16
            nl, nj, nk = 4, 2, 2
        elseif n==18
            nl, nj, nk = 2, 3, 3
        elseif n==27
            nl, nj, nk = 3, 3, 3
        end
    else
        if n in (0,4)
            nl, nj, nk = 2, 2, 1
        # elseif n==10
            # nl, nj, nk = 5, 2, 1
        elseif n==12
            nl, nj, nk = 4, 3, 1
        elseif n==9
            nl, nj, nk = 3, 3, 1
        elseif n==6
            nl, nj, nk = 3, 2, 1
        elseif n==8
            nl, nj, nk = 4, 2, 1
        elseif n==2
            nl, nj, nk = 2, 1, 1
        elseif n==3
            nl, nj, nk = 3, 1, 1
        end
    end

    ipL = get_ip_coords(LIN2, nl) # longitudinal
    ipT = get_ip_coords(LIN2, nj) # transversal
    
    resize!(elem.ips, nl*nj*nk)
    for i in 1:nl
        for j in 1:nj
            for k in 1:nk
                if ndim==2
                    R = [ ipL[i].coord[1], ipT[j].coord[1], 0.0 ]
                    w = ipL[i].w*ipT[j].w
                else
                    R = [ ipL[i].coord[1], ipT[j].coord[1], ipT[k].coord[1] ]
                    w = ipL[i].w*ipT[j].w*ipT[k].w
                end
                m = (i-1)*nj*nk + (j-1)*nk + k
                elem.ips[m] = Ip(R, w)
                elem.ips[m].id = m
                elem.ips[m].state = compat_state_type(typeof(elem.mat), typeof(elem), elem.env)(elem.env)
                elem.ips[m].owner = elem
            end
        end
    end

    # finding ips global coordinates
    C     = getcoords(elem)
    shape = elem.shape

    s = 0.0
    for ip in elem.ips
        # @show ip.w ip.R
        s+= ip.w
    end
    # @show s

    for ip in elem.ips
        R = [ ip.R[1], 0.0, 0.0 ]
        N = shape.func(R)
        ip.coord = C'*N
    end

end


function distributed_bc(elem::MechBeam, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function body_c(elem::MechBeam, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function elem_config_dofs(elem::MechBeam)
    ndim = elem.env.ndim
    ndim in (2,3) || error("MechBeam: Beam elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if ndim==3
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
        end
        add_dof(node, :rz, :mz)
    end
end


function elem_map(elem::MechBeam)
    if elem.env.ndim==2
        keys =(:ux, :uy, :rz)
    else
        keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    end
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
end


# Rotation Matrix
function set_rot_x_xp(elem::MechBeam, J::Matx, R::Matx)
    ndim = elem.env.ndim
    V1 = normalize(vec(J))
    V1 = round.(V1, digits=14)

    if ndim==2
        V2 = [ -V1[2], V1[1] ]
        R[1,:] .= V1
        R[2,:] .= V2
    else
        if V1[1]==0.0
            V2 = [ 1.0, 0.0, 0.0 ]
        elseif V1[2]==0.0
            V2 = [ 0.0, 1.0, 0.0 ]
        else
            V2 = [ 0.0, 0.0, 1.0 ]
        end

        V3 = cross(V1, V2)
        V2 = cross(V3, V1)

        normalize!(V1)
        normalize!(V2)
        normalize!(V3)

        R[1,:] .= V1
        R[2,:] .= V2
        R[3,:] .= V3
    end
end

function calcS(elem::MechBeam, αs::Float64)
    return @SMatrix [ 
        1.  0.  0.
        0.  αs   0.
        0.  0.  αs 
    ]
end

function setB(elem::MechBeam, ip::Ip, L::Matx, N::Vect, dNdX::Matx, Rθ::Matx, Bil::Matx, Bi::Matx, B::Matx)
    ndim = elem.env.ndim
    ndof = ndim==2 ? 3 : 6
    nnodes = size(dNdX,1)
    thz = elem.props.thz
    thy = elem.props.thy
    # Note that matrix B is designed to work with tensors in Mandel's notation
    if ndim==2
        Rθ[3,3] = 1.0
        for i in 1:nnodes
            η = ip.R[2]
            Rθ[1:2,1:2] .= L

            Ni = N[i]
            dNdx = dNdX[i]

            Bil[1,1] = dNdx;                              Bil[1,3] = -dNdx*η*thy/2
                                 Bil[3,2] = dNdx/SR2;     Bil[3,3] = -1/SR2*Ni

            c = (i-1)*ndof
            @mul Bi = Bil*Rθ
            B[:, c+1:c+ndof] .= Bi
            # @showm L
            # @showm Rθ
            # @showm Bil
            # error()
        end
    else
        for i in 1:nnodes
            η = ip.R[2]
            ζ = ip.R[3]
            Rθ[1:3,1:3] .= L
            # Rθ[1:3,1:3] .= elem.Dlmn[i]
            Rθ[4:6,4:6] .= elem.Dlmn[i]

            Ni = N[i]
            dNdx = dNdX[i]

            Bil[1,1] = dNdx;                                                                             Bil[1,5] = dNdx*ζ*thz/2;  Bil[1,6] = -dNdx*η*thy/2
                                                   Bil[2,3] = dNdx/SR2;  Bil[2,4] = 1/SR2*dNdx*η*thy/2;  Bil[2,5] = 1/SR2*Ni
                             Bil[3,2] = dNdx/SR2;                        Bil[3,4] = -1/SR2*dNdx*ζ*thz/2;                           Bil[3,6] = -1/SR2*Ni

            c = (i-1)*ndof
            @mul Bi = Bil*Rθ
            B[:, c+1:c+ndof] .= Bi

            # @showm L
            # @showm Rθ
            # @showm Bil
            # error()
        end
    end
end

function elem_stiffness(elem::MechBeam)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    thz = elem.props.thz
    thy = elem.props.thy
    ndof = ndim==2 ? 3 : 6
    # nstr = ndim==2 ? 2 : 3
    nstr = 3

    C   = getcoords(elem)
    K   = zeros(ndof*nnodes, ndof*nnodes)
    B   = zeros(nstr, ndof*nnodes)
    L   = zeros(ndim, ndim)
    Rθ  = zeros(ndof, ndof)
    Bil = zeros(nstr, ndof)
    Bi  = zeros(nstr, ndof)
    S   = calcS(elem, elem.props.αs)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D  = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        D     = calcD(elem.mat, ip.state)

        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)
        
        if ndim==2
            detJ′ = dx′dξ*thz*thy/2
        else
            detJ′ = dx′dξ*thz/2*thy/2
        end
        coef = detJ′*ip.w
        K += coef*B'*S*D*B
    end

    map = elem_map(elem)
    return K, map, map
end


function update_elem!(elem::MechBeam, U::Array{Float64,1}, dt::Float64)
    ndim = elem.env.ndim
    nnodes = length(elem.nodes)
    thz = elem.props.thz
    thy = elem.props.thy
    ndof = ndim==2 ? 3 : 6
    # nstr = ndim==2 ? 2 : 3
    nstr = 3

    C   = getcoords(elem)
    B   = zeros(nstr, ndof*nnodes)
    L   = zeros(ndim,ndim)
    Rθ  = zeros(ndof,ndof)
    Bil = zeros(nstr,ndof)
    Bi  = zeros(nstr,ndof)

    map = elem_map(elem)
    dU  = U[map]
    dF  = zeros(length(dU))
    Δε  = zeros(nstr)
    S   = calcS(elem, elem.props.αs)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D  = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)
        Δε = B*dU
        Δσ, status = update_state!(elem.mat, ip.state, Δε)
        failed(status) && return dF, map, status
        
        if ndim==2
            detJ′ = dx′dξ*thz*thy/2
        else
            detJ′ = dx′dξ*thz/2*thy/2
        end
        coef = detJ′*ip.w
        dF += coef*B'*S*Δσ
    end

    return dF, map, success()
end


function elem_vals(elem::MechBeam)
    # get ip average values
    ipvals = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
    sum  = merge(+, ipvals... )
    nips = length(elem.ips)
    vals = OrderedDict( k=>v/nips for (k,v) in sum)
    return vals
end


function elem_extrapolated_node_vals(elem::MechBeam)
    ndim = elem.env.ndim
    thz = elem.props.thz
    thy = elem.props.thy
    ndof = ndim==2 ? 3 : 6

    nnodes = length(elem.nodes)
    C = getcoords(elem) # global coordinates
    Ξ = elem.shape.nat_coords # natural coordinates

    # get local displacementhy from global
    if ndim==2
        keys =(:ux, :uy, :rz)
    else
        keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    end
    U = [ node.dofdict[key].vals[key] for node in elem.nodes for key in keys ]
    U′ = similar(U)
    Rθ = zeros(ndof,ndof)

    # coefficient matrices to compute moment and shear
    Am = zeros(nnodes, nnodes)
    Av = zeros(nnodes, nnodes)

    for i in 1:nnodes
        if ndim==2
            Rθ[1:2,1:2] .= elem.Dlmn[i]
            Rθ[3,3] = 1.0 
        else
            Rθ[1:3,1:3] .= elem.Dlmn[i]
            Rθ[4:6,4:6] .= elem.Dlmn[i]
        end

        U′[(i-1)*ndof+1:i*ndof] = Rθ*U[(i-1)*ndof+1:i*ndof]
        ξ = Ξ[i]

        # first derivatives and jacobian
        dNdξ = elem.shape.deriv([ξ])
        jac  = norm(C'*dNdξ)

        # second derivatives
        d2Ndξ2 = elem.shape.deriv2([ξ])
        d2xdξ2 = dot(C'*d2Ndξ2, elem.Dlmn[i][1,:])
        d2ξNdx2 = -d2xdξ2/jac^3
        
        for j in 1:nnodes
            dNdx = dNdξ[j]/jac
            Am[i,j] = dNdx
            d2Ndx2 = d2Ndξ2[j]*(1/jac)^2 + dNdξ[j]*d2ξNdx2
            Av[i,j] = d2Ndx2
        end
        
    end

    E = elem.mat.E
    if ndim==2
        θZ = U′[3:ndof:ndof*nnodes]
        Izy = thz*thy^3/12

        # Bending moment
        Mxy = (Am*θZ).*(E*Izy)
        # Shear (negative to match convention)
        Vxy = -(Av*θZ).*(E*Izy)
        Vxy = ones(nnodes)*mean(Vxy)

        return OrderedDict(
            :MXY => Mxy,
            :VXY => Vxy
        )
    else
        θY = U′[5:ndof:ndof*nnodes-1]
        θZ = U′[6:ndof:ndof*nnodes]

        Iyz = thy*thz^3/12
        Izy = thz*thy^3/12

        # Bending moment
        Mxz = (Am*θY).*(E*Iyz)
        Mxy = (Am*θZ).*(E*Izy)
        # Shear (negative to match convention)
        Vxz = -(Av*θY).*(E*Iyz) 
        Vxy = -(Av*θZ).*(E*Izy)
        Vxz = ones(nnodes)*mean(Vxz)
        Vxy = ones(nnodes)*mean(Vxy)

        return OrderedDict(
            :MXZ => Mxz,
            :MXY => Mxy,
            :VXZ => Vxz,
            :VXY => Vxy
        )
    end
end
 