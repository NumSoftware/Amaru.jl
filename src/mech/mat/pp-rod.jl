# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PPRod


"""
    PPRod

A type for linear elastic perfectly plastic materials in rods.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PPRod<:Material
    "Young modulus"
    E::Float64
    "Section area"
    A::Float64
    "Yielding stress"
    σy0::Float64
    "Hardening parameter"
    H::Float64
    "Density"
    ρ::Float64

    function PPRod(prms::Dict{Symbol,Float64})
        return  PPRod(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `PPRod` material type

    # Arguments
    - `E`: Young modulus
    - `A`: Section area
    - `dm`: Diameter (only if `A` is not provided)
    - `fy`: Yielding stress
    - `H`: Hardening parameter
    - `rho`: Density
    """
    function PPRod(;E=NaN, A=NaN, fy=NaN, sig_y=NaN, H=0.0, rho=0.0, dm=NaN)
        !isnan(dm) && (A=pi*dm^2/4)
        isnan(fy) && (fy=sig_y)
        @check E>0.0     
        @check A>0.0     
        @check fy>0
        @check rho>=0.0


        return new(E, A, fy, H, rho)
    end
end

"""
    ElasticRodState

A type for the state data of a PPRod type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PPRodState<:IpState
    env::ModelEnv
    σ::Float64
    ε::Float64
    εp::Float64
    Δγ ::Float64

    function PPRodState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        this.εp = 0.0
        this.Δγ = 0.0
        return this
    end
end

matching_elem_type(::PPRod) = MechRod
matching_elem_type_if_embedded(::PPRod) = MechEmbRod

# Type of corresponding state structure
ip_state_type(mat::PPRod) = PPRodState

function yield_func(mat::PPRod, ipd::PPRodState, σ::Float64)
    σya = mat.σy0 + mat.H*ipd.εp
    return abs(σ) - σya
end

function calcD(mat::PPRod, ipd::PPRodState)
    if ipd.Δγ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end

function stress_update(mat::PPRod, ipd::PPRodState, Δε::Float64)
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

function ip_state_vals(mat::PPRod, ipd::PPRodState)
    return OrderedDict{Symbol,Float64}(
        :sa  => ipd.σ,
        :ea  => ipd.ε,
        :eap => ipd.εp,
        :fa  => ipd.σ*mat.A,
        :A   => mat.A
    )
end

