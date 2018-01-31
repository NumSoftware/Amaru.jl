# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticBeam

mutable struct BeamIpState<:IpState
    shared_data::SharedAnalysisData
    #σ::Float64
    #ε::Float64
    function BeamIpState(shared_data::SharedAnalysisData=SharedAnalysisData())
        return new(shared_data)
    end
end

mutable struct ElasticBeam<:Material
    E::Float64
    A::Float64
    I::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, A=NaN, I=NaN, Ix=NaN, Iy=NaN, Iz=NaN)
        @assert E>0.0
        @assert A>0.0
        @assert I>0.0
        this = new(E,A,I)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeam

# Create a new instance of Ip data
new_ip_state(mat::ElasticBeam, shared_data::SharedAnalysisData) = BeamIpState(shared_data)

#function set_state(ipd::BeamIpState, σ=NaN, ε=NaN)
    #if !isnan(σ); ipd.σ = σ end
    #if !isnan(ε); ipd.ε = ε end
#end

function ip_state_vals(mat::ElasticBeam, ipd::BeamIpState)
    return Dict()
end
