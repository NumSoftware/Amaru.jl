# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechRod

struct MechRodProps<:ElemProperties
    ρ::Float64
    γ::Float64
    A::Float64

    function MechRodProps(;rho=0.0, gamma=0.0, A=0.0, diameter=0.0)
        @check rho>=0
        @check gamma>=0
        @check A>0.0 || diameter>0.0
        A>0 || (A=π*diameter^2/4)

        return new(rho, gamma, A)
    end    
end

MechRod = MechRodProps

"""
    MechRod

A line finite element for mechanical equilibrium analyses.
"""
mutable struct MechRodElem<:MechElem
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    props ::MechRodProps
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function MechRodElem(props=MechRodProps())
        this = new()
        this.props = props
        return this
    end
end

matching_shape_family(::Type{MechRodElem}) = LINECELL
matching_elem_type(::Type{MechRodProps}) = MechRodElem
matching_props_type(::Type{MechRodElem}) = MechRodProps

function elem_stiffness(elem::MechRodElem)
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
        @gemm J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        E    = calcD(elem.matparams, ip.state)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end
    keys = [:ux, :uy, :uz][1:ndim]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_mass(elem::MechRodElem)

    ndim   = elem.env.ndim
    nnodes = length(elem.nodes)
    ρ = elem.matparams.ρ
    A = elem.props.A

    C = getcoords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, 1)
    N = zeros(ndim, ndim*nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        Ni = elem.shape.func(ip.R)
        setNt(ndim,Ni,N)

        @gemm J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*A*detJ*ip.w
        @gemm M += coef*N'*N

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



function distributed_bc(elem::MechRodElem, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end

function body_c(elem::MechRodElem, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end

function elem_internal_forces(elem::MechRodElem, F::Array{Float64,1})
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
        @gemm J = C'*dNdR
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

function elem_activate(elem::MechRodElem, F::Array{Float64,1})
    elem_internal_forces(elem, F)
end



function update_elem!(elem::MechRodElem, U::Array{Float64,1}, Δt::Float64)

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
        @gemm J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        deps = (B*dU)[1]
        dsig, _ = update_state(elem.matparams, ip.state, deps)
        coef = A*detJ*ip.w
        dF  .+= coef*vec(B')*dsig
    end

    return dF, map, success()
end

function elem_vals(elem::MechRodElem)
    # get ip average values
    ipvals = [ ip_state_vals(elem.matparams, ip.state) for ip in elem.ips ]
    sum  = merge(+, ipvals... )
    nips = length(elem.ips)
    vals = OrderedDict( k=>v/nips for (k,v) in sum)
    return vals
end
