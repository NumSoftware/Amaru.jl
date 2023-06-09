# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    MechSolid

A bulk finite element for mechanical equilibrium analyses.
"""
mutable struct MechSolid<:Mechanical
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env   ::ModelEnv

    function MechSolid();
        return new()
    end
end

matching_shape_family(::Type{MechSolid}) = BULKCELL

# function elem_init(elem::MechSolid)
#     ipdata_ty = typeof(elem.ips[1].state)
#     if :h in fieldnames(ipdata_ty)
#         # Element volume/area
#         V = 0.0
#         C = getcoords(elem)
#         for ip in elem.ips
#             dNdR = elem.shape.deriv(ip.R)
#             J    = dNdR*C
#             detJ = det(J)
#             @assert detJ>0
#             V   += detJ*ip.w
#         end

#         # Representative length size for the element
#         nips = length(elem.ips)
#         ndim = elem.env.ndim
#         h = V^(1/ndim)

#         for ip in elem.ips
#             ip.state.h = h
#         end
#     end

#     return nothing
# end


function distributed_bc(elem::MechSolid, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_solid_boundary_forces(elem, facet, key, val)
end


function body_c(elem::MechSolid, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_solid_body_forces(elem, key, val)
end


function setB(elem::Element, ip::Ip, dNdX::Matx, B::Matx)
    nnodes, ndim = size(dNdX)
    # Note that matrix B is designed to work with tensors in Mandel's notation

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[i,1]
            B[2,2+j*ndim] = dNdX[i,2]
            B[6,1+j*ndim] = dNdX[i,2]/SR2
            B[6,2+j*ndim] = dNdX[i,1]/SR2
        end
        if elem.env.modeltype=="axisymmetric"
            N = elem.shape.func(ip.R)
            for i in 1:nnodes
                j = i-1
                r = ip.coord.x
                B[1,1+j*ndim] = dNdX[i,1]
                B[2,2+j*ndim] = dNdX[i,2]
                B[3,1+j*ndim] =    N[i]/r
                B[6,1+j*ndim] = dNdX[i,2]/SR2
                B[6,2+j*ndim] = dNdX[i,1]/SR2
            end
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[i,1]
            dNdy = dNdX[i,2]
            dNdz = dNdX[i,3]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,2+j*ndim] = dNdz/SR2;   B[4,3+j*ndim] = dNdy/SR2
            B[5,1+j*ndim] = dNdz/SR2;   B[5,3+j*ndim] = dNdx/SR2
            B[6,1+j*ndim] = dNdy/SR2;   B[6,2+j*ndim] = dNdx/SR2
        end
    end

end


function elem_stiffness(elem::MechSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = C'*dNdR
        @gemm dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        setB(elem, ip, dNdX, B)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.mat, ip.state)
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


function elem_mass(elem::MechSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    ρ = elem.mat.ρ
    C = getcoords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    N = zeros(ndim, nnodes*ndim)
    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute N matrix
        Ni   = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        for i in 1:nnodes
            for j in 1:ndim
                N[j, (i-1)*ndim+j] = Ni[i]
            end
        end

        @gemm J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*detJ*ip.w*th
        @gemm M += coef*N'*N
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return M, map, map
end


function elem_internal_forces(elem::MechSolid)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

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
        setB(elem, ip, dNdX, B)

        σ    = ip.state.σ
        coef = detJ*ip.w*th
        @gemv dF += coef*B'*σ
    end
    
    return dF, map, success()
end


function elem_update!(elem::MechSolid, U::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    dU = U[map]
    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)

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
        setB(elem, ip, dNdX, B)

        @gemv Δε = B*dU
        Δσ, status = stress_update(elem.mat, ip.state, Δε)
        failed(status) && return dF, map, failure("MechSolid: Error at integration point $(ip.id)")
        coef = detJ*ip.w*th
        @gemv dF += coef*B'*Δσ
    end

    return dF, map, success()
end


function elem_vals(elem::MechSolid)
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


