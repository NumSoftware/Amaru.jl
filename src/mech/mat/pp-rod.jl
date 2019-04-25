# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PPRod

mutable struct PPRodIpState<:IpState
    env::ModelEnv
    σ::Float64
    ε::Float64
    εpa::Float64
    Δγ ::Float64

    function PPRodIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        this.εpa = 0.0
        this.Δγ  = 0.0
        return this
    end
end

mutable struct PPRod<:Material
    E::Float64
    A::Float64
    σy0::Float64
    H::Float64
    ρ::Float64

    function PPRod(prms::Dict{Symbol,Float64})
        return  PPRod(;prms...)
    end

    function PPRod(;E=NaN, A=NaN, sig_y=NaN, H=0.0, rho::Number=0.0)
        E<=0.0 && error("Invalid value for E: $E")
        A<=0.0 && error("Invalid value for A: $A")
        sig_y<0.0 && error("Invalid value for sig_y: $sig_y")
        rho<0.0 && error("Invalid value for rho: $rho")

        return new(E, A, sig_y, H, rho)
    end
end

matching_elem_type(::PPRod) = MechRod

# Create a new instance of Ip data
new_ip_state(mat::PPRod, env::ModelEnv) = PPRodIpState(env)

function yield_func(mat::PPRod, ipd::PPRodIpState, σ::Float64)
    σya = mat.σy0 + mat.H*ipd.εpa
    return abs(σ) - σya
end

function calcD(mat::PPRod, ipd::PPRodIpState)
    if ipd.Δγ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end

function stress_update(mat::PPRod, ipd::PPRodIpState, Δε::Float64)
    E, H = mat.E, mat.H
    σini = ipd.σ
    σtr    = σini + E*Δε
    ftr    = yield_func(mat, ipd, σtr)
    ipd.Δγ = ftr>0.0 ? ftr/(E+H) : 0.0
    Δεp    = ipd.Δγ*sign(σtr)
    ipd.εpa += ipd.Δγ
    ipd.σ  = σtr - E*Δεp
    Δσ     = ipd.σ - σini
    ipd.ε += Δε
    return Δσ
end

function ip_state_vals(mat::PPRod, ipd::PPRodIpState)
    return OrderedDict{Symbol,Float64}(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :epa_ppt => ipd.εpa,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A )
end

