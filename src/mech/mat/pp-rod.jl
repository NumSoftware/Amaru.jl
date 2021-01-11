# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PPRod

mutable struct PPRodIpState<:IpState
    env::ModelEnv
    σ::Float64
    ε::Float64
    εp::Float64
    Δγ ::Float64

    function PPRodIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        this.εp = 0.0
        this.Δγ = 0.0
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

    function PPRod(;E=NaN, A=NaN, fy=NaN, sig_y=NaN, H=0.0, rho=0.0)
        @check E>0.0     
        @check A>0.0     
        @check sig_y>0.0 || fy>0
        @check rho>=0.0

        isnan(fy) && (fy=sig_y)

        return new(E, A, fy, H, rho)
    end
end

matching_elem_type(::PPRod) = MechRod
matching_elem_type_if_embedded(::PPRod) = MechEmbRod

# Type of corresponding state structure
ip_state_type(mat::PPRod) = PPRodIpState

function yield_func(mat::PPRod, ipd::PPRodIpState, σ::Float64)
    σya = mat.σy0 + mat.H*ipd.εp
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
    E, H    = mat.E, mat.H
    σini    = ipd.σ
    σtr     = σini + E*Δε
    ftr     = yield_func(mat, ipd, σtr)
    ipd.Δγ  = ftr>0.0 ? ftr/(E+H) : 0.0
    Δεp     = ipd.Δγ*sign(σtr)
    ipd.εp += ipd.Δγ
    ipd.σ   = σtr - E*Δεp
    Δσ      = ipd.σ - σini
    ipd.ε  += Δε
    return Δσ, ReturnStatus(true)
end

function ip_state_vals(mat::PPRod, ipd::PPRodIpState)
    return OrderedDict{Symbol,Float64}(
        :sa  => ipd.σ,
        :ea  => ipd.ε,
        :eap => ipd.εp,
        :Fa  => ipd.σ*mat.A,
        :A   => mat.A
    )
end

