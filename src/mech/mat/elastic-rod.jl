# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticRod

mutable struct ElasticRodIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Float64
    ε::Float64
    function ElasticRodIpState(shared_data::SharedAnalysisData=SharedAnalysisData())
        this = new(shared_data)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end

mutable struct ElasticRod<:Material
    E::Float64
    A::Float64
    ρ::Float64

    function ElasticRod(prms::Dict{Symbol,Float64})
        return  ElasticRod(;prms...)
    end

    function ElasticRod(;E::Number=NaN, A::Number=NaN, rho::Number=0.0)
        E<=0.0 && error("Invalid value for E: $E")
        A<=0.0 && error("Invalid value for A: $A")
        rho<0.0 && error("Invalid value for rho: $rho")
        return new(E, A, rho)
    end
end

matching_elem_type(::ElasticRod) = MechRod

# Create a new instance of Ip data
new_ip_state(mat::ElasticRod, shared_data::SharedAnalysisData) = ElasticRodIpState(shared_data)

function set_state(ipd::ElasticRodIpState, σ=NaN, ε=NaN)
    if !isnan(σ); ipd.σ = σ end
    if !isnan(ε); ipd.ε = ε end
end

function stress_update(mat::ElasticRod, ipd::ElasticRodIpState, Δε::Float64)
    E  = mat.E
    Δσ = mat.E*Δε
    ipd.ε += Δε
    ipd.σ += Δσ
    return Δσ
end

function ip_state_vals(mat::ElasticRod, ipd::ElasticRodIpState)
    return OrderedDict(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A )
end

function calcD(mat::ElasticRod, ips::ElasticRodIpState)
    return mat.E
end

