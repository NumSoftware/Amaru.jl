export VonMises


"""
    VonMises

A type for linear elastic materials with Von Mises failure criterion.

# Fields

$(TYPEDFIELDS)
"""
mutable struct VonMises<:Material
    "Young modulus"
    E ::Float64
    "Poisson ratio"
    ν ::Float64
    "Yielding stress"
    σy::Float64
    "Hardening parameter"
    H ::Float64
    "Density"
    ρ::Float64

    function VonMises(prms::Dict{Symbol,Float64})
        return VonMises(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `VonMises` material type

    # Arguments
    - `E`: Young modulus
    - `nu`: Poisson ratio
    - `fy`: Yielding stress
    - `H`: Hardening parameter
    """
    function VonMises(; params...)
        names = (E="Young modulus", nu="Poisson ratio", fy="Yield stress", H="Hardening modulus")
        required = (:E, :nu, :fy, :H)
        @checkmissing params required names

        default = (nu=0.0,)
        params  = merge(default, params)
        E       = params.E
        nu      = params.nu
        fy      = params.fy
        H       = params.H

        @check E>=0.0
        @check 0<=nu<0.5
        @assert fy>0.0
        @assert H>=0.0

        return new(E, nu, fy, H)
    end
end

"""
    DruckerPragerState

A type for the state data of a `DruckerPrager` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct VonMisesState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Vec6
    "Strain tensor"
    ε::Vec6
    "Accumulated plastic strain"
    εpa::Float64
    "Plastic multiplier"
    Δγ::Float64
    function VonMisesState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end


# Type of corresponding state structure
ip_state_type(::Type{VonMises}) = VonMisesState

# Element types that work with this material
matching_elem_types(::Type{VonMises}) = (MechSolid, MechShell)

function yield_func(mat::VonMises, state::VonMisesState, σ::AbstractArray)
    j1  = J1(σ)
    j2d = J2D(σ)
    σy  = mat.σy
    H   = mat.H
    εpa = state.εpa
    return √(3*j2d) - σy - H*εpa
end


function calcD(mat::VonMises, state::VonMisesState, stressmodel::String=state.env.ana.stressmodel)
    σy = mat.σy
    H  = mat.H
    #De = mat.De
    De  = calcDe(mat.E, mat.ν, stressmodel)
    @show De

    if state.Δγ==0.0
        return De
    end

    j2d = J2D(state.σ)
    @assert j2d>0
    s  = dev(state.σ)
    su = s/norm(s)
    V  = √(3/2)*su # df/dσ
    Nu = su

    # return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
    return De - De*Nu*V'*De / (V'*De*Nu + H)

end


function update_state!(mat::VonMises, state::VonMisesState, Δε::Array{Float64,1}, stressmodel::String=state.env.ana.stressmodel)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, stressmodel)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, state, σtr)

    if ftr < 1.e-8
        # elastic
        state.Δγ = 0.0
        state.σ  = σtr
    else
        # plastic
        K, G  = mat.E/(3.0*(1.0-2.0*mat.ν)), mat.E/(2.0*(1.0+mat.ν))
        H     = mat.H
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        # @show √j2dtr
        # @show state.Δγ

        √j2dtr - state.Δγ*√2*G <= 0.0 && return state.σ, failure("VonMisses: √j2dtr - state.Δγ*√2*G>0.0")
        state.Δγ = ftr/(√6*G + H)
        j1     = j1tr
        m      = 1.0 - state.Δγ*√2*G/√j2dtr
        state.σ  = m*dev(σtr) + j1/3.0*tI

        state.εpa += state.Δγ

    end

    state.ε += Δε
    Δσ     = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::VonMises, state::VonMisesState, stressmodel::String=state.env.ana.stressmodel)
    ndim  = state.env.ndim
    σ, ε  = state.σ, state.ε
    j1    = tr(σ)
    srj2d = √J2D(σ)

    D = stress_strain_dict(σ, ε, stressmodel)
    D[:epa]   = state.εpa
    D[:j1]    = j1
    D[:srj2d] = srj2d

    return D
end
