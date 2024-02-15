# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechBar, MechRod, MechTruss

struct MechBarProps<:ElemProperties
    ρ::Float64
    γ::Float64
    A::Float64

    function MechBarProps(; props...)
        names = (rho="Density", gamma="Specific weight", A="Section area")
        default = (rho=0.0, gamma=0.0)
        props   = merge(default, props)
        required = (:A,)
        @checkmissing props required names

        @check props.rho>=0
        @check props.gamma>=0
        @check props.A>0.0

        return new(props.rho, props.gamma, props.A)
    end    
end

# 

"""
    MechBar

A line finite element for mechanical equilibrium analyses.
"""
mutable struct MechBar<:Mech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat::Material
    props ::MechBarProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function MechBar()
        return new()
    end
end

const MechRod = MechBar
const MechTruss = MechBar

compat_shape_family(::Type{MechBar}) = LINECELL
# compat_elem_types(::Type{MechBarProps}) = MechBar
compat_elem_props(::Type{MechBar}) = MechBarProps

embedded_type(::Type{MechBar}) = MechEmbBar
# embedded_type(::Type{MechEmbBar}) = MechEmbBar


# function check_props(::Type{MechSolid}; props...)
#     names = (rho="Density", gamma="Specific weight", A="Section area", diameter="diameter")
#     required = (:A,)
#     @checkmissing props required names

#     default = (rho=0.0, gamma=0.0)
#     props   = merge(default, props)
#     rho     = props.rho
#     gamma   = props.gamma

#     @check rho>=0
#     @check gamma>=0
#     @check A>0

#     return props
# end


function elem_stiffness(elem::MechBar)
    local E::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)

    A = elem.props.A
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        E    = calcD(elem.mat, ip.state)
        coef = E*A*detJ*ip.w
        @mul K += coef*B'*B
    end
    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_mass(elem::MechBar)

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    ρ = elem.props.ρ
    A = elem.props.A

    C = getcoords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, 1)
    N = zeros(ndim, ndim*nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        Ni = elem.shape.func(ip.R)
        setNt(ndim,Ni,N)

        @mul J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*A*detJ*ip.w
        @mul M += coef*N'*N

    end

    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
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


function distributed_bc(elem::MechBar, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function body_c(elem::MechBar, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function elem_internal_forces(elem::MechBar, F::Array{Float64,1})
    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    A      = elem.props.A
    keys   = [:ux, :uy, :uz][1:ndim]
    map    = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dF = zeros(nnodes*ndim)
    C = getcoords(elem)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        σ = ip.state.σ
        coef = A*detJ*ip.w
        dF .+= coef*σ*vec(B')
    end

    F[map] += dF
end


function elem_activate(elem::MechBar, F::Array{Float64,1})
    elem_internal_forces(elem, F)
end


function update_elem!(elem::MechBar, U::Array{Float64,1}, Δt::Float64)

    ndim   = elem.env.ndim 
    nnodes = length(elem.nodes)
    A      = elem.props.A
    keys   = [:ux, :uy, :uz][1:ndim]
    map    = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    C  = getcoords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, 1)


    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        deps = (B*dU)[1]
        dsig, _ = update_state!(elem.mat, ip.state, deps)
        coef = A*detJ*ip.w
        dF  .+= coef*vec(B')*dsig
    end

    return dF, map, success()
end


function elem_vals(elem::MechBar)
    # get ip average values
    ipvals = [ ip_state_vals(elem.mat, ip.state) for ip in elem.ips ]
    sum  = merge(+, ipvals... )
    nips = length(elem.ips)
    vals = OrderedDict( k=>v/nips for (k,v) in sum)
    return vals
end
