# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechFluid

struct MechFluidProps<:ElemProperties
    ρ::Float64
    γ::Float64

    function MechFluidProps(; props...)
        default = (rho=0.0, gamma=0.0)
        props   = merge(default, props)
        rho     = props.rho
        gamma   = props.gamma

        @check rho>=0
        @check gamma>=0

        return new(rho, gamma)
    end    
end


"""
    MechFluid

A bulk finite element for fluid analyses.
"""
mutable struct MechFluid<:Mech
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::MechFluidProps
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv

    function MechFluid()
        return new()
    end
end

compat_shape_family(::Type{MechFluid}) = BULKCELL
compat_elem_props(::Type{MechFluid}) = MechFluidProps


function distributed_bc(elem::MechFluid, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_boundary_forces(elem, facet, key, val)
end


function body_c(elem::MechFluid, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_solid_body_forces(elem, key, val)
end


function setB(::MechFluid, dNdX::Matx, B::Matx)
    nnodes, ndim = size(dNdX)
    # Note that matrix B is designed to work with tensors in Mandel's notation

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[i,1]
            B[1,2+j*ndim] = dNdX[i,2]
        end
    else
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[i,1]
            B[1,2+j*ndim] = dNdX[i,2]
            B[1,3+j*ndim] = dNdX[i,3]
        end
    end

end


function elem_stiffness(elem::MechFluid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)

    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    
    B    = zeros(1, nnodes*ndim)
    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        setB(elem, dNdX, B)

        # compute K
        Kb = elem.mat.K # bulk modulus
        coef = Kb*detJ*ip.w*th
        
        @mul K += coef*B'*B
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


function elem_mass(elem::MechFluid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    ρ = elem.props.ρ
    C = getcoords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    N = zeros(ndim, nnodes*ndim)
    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.ana.stressmodel==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute N matrix
        Ni   = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        for i in 1:nnodes
            for j in 1:ndim
                N[j, (i-1)*ndim+j] = Ni[i]
            end
        end

        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*detJ*ip.w*th
        @mul M += coef*N'*N
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return M, map, map
end


function elem_internal_forces(elem::MechFluid)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    C = getcoords(elem)
    for ip in elem.ips
        if elem.env.ana.stressmodel==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        setB(elem, ip, dNdX, B)

        p    = ip.state.p
        coef = detJ*ip.w*th
        @mul dF += coef*B'*p
    end
    
    return dF, map, success()
end


function update_elem!(elem::MechFluid, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.ana.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    B  = zeros(1, nnodes*ndim)

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    C = getcoords(elem)
    for ip in elem.ips

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        setB(elem, dNdX, B)

        Δεv = dot(B,dU)
        Δp, status = update_state!(elem.mat, ip.state, Δεv)
        failed(status) && return dF, map, status
        coef = detJ*ip.w*th
        @mul dF += coef*B'*Δp
    end

    return dF, map, success()
end


function elem_vals(elem::MechFluid)
    vals = OrderedDict{Symbol,Float64}()
    return vals
end


