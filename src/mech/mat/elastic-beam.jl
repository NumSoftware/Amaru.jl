# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticBeam

mutable struct ElasticBeamState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticBeamState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(3)
        this.ε = zeros(3)
        return this
    end
end

mutable struct ElasticBeam<:Material
    E::Float64
    nu::Float64
    th::Float64
    ts::Float64
    A::Float64
    ρ::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, nu=0.0, th=NaN, ts=NaN, rho=0.0)
        @check E>0.0
        @check 0.0<=nu<0.5
        @check th>0.0
        @check ts>0.0
        @check rho>=0.0
        this = new(E, nu, ts, th, ts*th, rho)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeam

# Type of corresponding state structure
ip_state_type(mat::ElasticBeam) = ElasticBeamState


function calcD(mat::ElasticBeam, ipd::ElasticBeamState)
    E = mat.E
    ν = mat.nu
    c = E/(1.0-ν^2)
    g = E/(1+ν)

    return [ c    0.0    0.0  
             0.0  5/6*g  0.0  
             0.0  0.0    5/6*g ]
end


function stress_update(mat::ElasticBeam, ipd::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(mat, ipd)
    dσ = D*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::ElasticBeam, ipd::ElasticBeamState)
    return OrderedDict(
      :(sx') => ipd.σ[1],
      :(ex') => ipd.ε[1],
      :A  => mat.A )
end
