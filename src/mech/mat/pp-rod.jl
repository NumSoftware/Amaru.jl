# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PPRod


"""
    PPRod

A type for linear elastic perfectly plastic materials in rods.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PPRod<:MatParams
    "Young modulus"
    E::Float64
    "Yielding stress"
    fy::Float64
    "Hardening parameter"
    H::Float64

    function PPRod(prms::Dict{Symbol,Float64})
        return  PPRod(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `PPRod` material type

    # Arguments
    - `E`: Young modulus
    - `fy`: Yielding stress
    - `H`: Hardening parameter
    """
    function PPRod(;E=NaN, fy=NaN, H=0.0)
        @check E>0.0
        @check fy>0
        return new(E, fy, H)
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

matching_elem_type(::PPRod) = MechRodElem
matching_elem_type_if_embedded(::PPRod) = MechEmbRodElem

# Type of corresponding state structure
ip_state_type(::MechRodElem, ::PPRod) = PPRodState
ip_state_type(::MechEmbRodElem, ::PPRod) = PPRodState


function yield_func(matparams::PPRod, state::PPRodState, σ::Float64)
    σya = matparams.fy + matparams.H*state.εp
    return abs(σ) - σya
end


function calcD(matparams::PPRod, state::PPRodState)
    if state.Δγ == 0.0
        return matparams.E
    else
        E, H = matparams.E, matparams.H
        return E*H/(E+H)
    end
end


function update_state(matparams::PPRod, state::PPRodState, Δε::Float64)
    E, H    = matparams.E, matparams.H
    σini    = state.σ
    σtr     = σini + E*Δε
    ftr     = yield_func(matparams, state, σtr)
    state.Δγ  = ftr>0.0 ? ftr/(E+H) : 0.0
    Δεp     = state.Δγ*sign(σtr)
    state.εp += state.Δγ
    state.σ   = σtr - E*Δεp
    Δσ      = state.σ - σini
    state.ε  += Δε
    return Δσ, ReturnStatus(true)
end


function ip_state_vals(matparams::PPRod, state::PPRodState)
    return OrderedDict{Symbol,Float64}(
        :sa  => state.σ,
        :ea  => state.ε,
        :eap => state.εp,
    )
end

