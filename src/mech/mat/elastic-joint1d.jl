# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint1D

mutable struct Joint1DIpState<:IpState
    shared_data::SharedAnalysisData
    ndim::Int
    sig ::Array{Float64,1}
    eps ::Array{Float64,1}
    function Joint1DIpState(shared_data::SharedAnalysisData=SharedAnalysisData())
        this = new(shared_data)
        ndim = shared_data.ndim
        this.sig = zeros(ndim)
        this.eps = zeros(ndim)
        return this
    end
end

mutable struct ElasticJoint1D<:Material
    ks::Float64
    kn::Float64
    h ::Float64    # section perimeter

    function ElasticJoint1D(prms::Dict{Symbol,Float64})
        return  ElasticJoint1D(;prms...)
    end

    function ElasticJoint1D(;ks=NaN, kn=NaN, h=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # h : section perimeter
        @assert ks>=0
        @assert kn>=0
        @assert (h>0 || A>0 || dm>0)

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(ks, kn, h)
        return this
    end
end

# Returns the element type that works with this material
matching_elem_type(::ElasticJoint1D) = MechJoint1D

# Create a new instance of Ip data
new_ip_state(mat::ElasticJoint1D, shared_data::SharedAnalysisData) = Joint1DIpState(shared_data)

function set_state(ipd::Joint1DIpState, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("MecElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("MecElasticSolid: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::ElasticJoint1D, ipd::Joint1DIpState)
    ks = mat.ks
    kn = mat.kn
    if ipd.shared_data.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function stress_update(mat::ElasticJoint1D, ipd::Joint1DIpState, deps)
    D = calcD(mat, ipd)
    dsig = D*deps

    ipd.eps += deps
    ipd.sig += dsig
    return dsig
end

function ip_state_vals(mat::ElasticJoint1D, ipd::Joint1DIpState)
    return Dict(
      :ur   => ipd.eps[1] ,
      :tau  => ipd.sig[1] )
end
