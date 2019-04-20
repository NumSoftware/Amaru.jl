# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointIpState<:IpState
    analysis_data::AnalysisData
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointIpState(analysis_data::AnalysisData=AnalysisData())
        this = new(analysis_data)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

mutable struct ElasticJoint<:Material
    E::Float64 # Young modulus from bulk material
    ν::Float64 # Poisson ration from bulk material
    α::Float64 # elastic displacement scale factor

    function ElasticJoint(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJoint(;E=NaN, nu=NaN, alpha=1.0)
        @assert E>=0
        @assert nu>=0

        this = new(E, nu, alpha)
        return this
    end
end

# Returns the element type that works with this material model
@static if @isdefined MechJoint
    matching_elem_type(::ElasticJoint) = MechJoint
end

# Create a new instance of Ip data
new_ip_state(mat::ElasticJoint, analysis_data::AnalysisData) = JointIpState(analysis_data)


function mountD(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.analysis_data.ndim
    G  = mat.E/(1.0+mat.ν)/2.0
    kn = mat.E*mat.α/ipd.h
    ks =     G*mat.α/ipd.h
    if ndim==2
        return [  kn  0.0 
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJoint, ipd::JointIpState, Δu)
    ndim = ipd.analysis_data.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ
    return Δσ
end

function ip_state_vals(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.analysis_data.ndim
    if ndim == 2
        return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] )
    else
        return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] )
    end
end
