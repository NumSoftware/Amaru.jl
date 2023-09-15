# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PPBar


"""
    PPBar

A type for linear elastic perfectly plastic materials in rods.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PPBar<:Material
    "Young modulus"
    E::Float64
    "Yielding stress"
    fy::Float64
    "Hardening parameter"
    H::Float64

    @doc """
        $(SIGNATURES)

    Creates an `PPBar` material type

    # Arguments
    - `E`: Young modulus
    - `fy`: Yielding stress
    - `H`: Hardening parameter
    """
    function PPBar(; params...)
        names = (E="Young modulus", fy="Yielding stress", H="Hardening modulus")
        required = (:E, :fy)    
        @checkmissing params required names

        default = (H=0.0,)
        params  = merge(default, params)
        E       = params.E
        fy      = params.fy
        H       = params.H

        @check E>=0.0
        @check H>=0.0
        @check fy>0.0
        this = new(E, fy, H)
        return this
    end
end

"""
    PPBarState

A type for the state data of a PPBar type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PPBarState<:IpState
    env::ModelEnv
    σ::Float64
    ε::Float64
    εp::Float64
    Δγ ::Float64

    function PPBarState(env::ModelEnv)
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        this.εp = 0.0
        this.Δγ = 0.0
        return this
    end
end

compat_state_type(::Type{PPBar}, ::Type{MechBar}, env::ModelEnv) = PPBarState
compat_state_type(::Type{PPBar}, ::Type{MechEmbBar}, env::ModelEnv) = PPBarState


# Type of corresponding state structure
compat_state_type(::Type{PPBar}) = PPBarState

# Element types that work with this material
compat_elem_types(::Type{PPBar}) = (MechBar,)
compat_elem_type_if_embedded(::Type{PPBar}) = MechEmbBar


function yield_func(mat::PPBar, state::PPBarState, σ::Float64)
    σya = mat.fy + mat.H*state.εp
    return abs(σ) - σya
end


function calcD(mat::PPBar, state::PPBarState)
    if state.Δγ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end


function update_state!(mat::PPBar, state::PPBarState, Δε::Float64)
    E, H    = mat.E, mat.H
    σini    = state.σ
    σtr     = σini + E*Δε
    ftr     = yield_func(mat, state, σtr)
    state.Δγ  = ftr>0.0 ? ftr/(E+H) : 0.0
    Δεp     = state.Δγ*sign(σtr)
    state.εp += state.Δγ
    state.σ   = σtr - E*Δεp
    Δσ      = state.σ - σini
    state.ε  += Δε
    return Δσ, ReturnStatus(true)
end


function ip_state_vals(mat::PPBar, state::PPBarState)
    return OrderedDict{Symbol,Float64}(
        :sX  => state.σ,
        :eX  => state.ε,
        :eXp => state.εp,
    )
end

