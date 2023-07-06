export VonMises


"""
    VonMises

A type for linear elastic materials with Von Mises failure criterion.

# Fields

$(TYPEDFIELDS)
"""
mutable struct VonMises<:MatParams
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
    - `rho`: Density
    """
    function VonMises(;E=NaN, nu=0.0, fy=0.0, H=0.0, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert fy>0.0
        @assert H>=0.0

        this    = new(E, nu, fy, H)
        #this.De = calcDe(E, nu)
        return this
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
    σ::Tensor2
    "Strain tensor"
    ε::Tensor2
    "Accumulated plastic strain"
    εpa::Float64
    "Plastic multiplier"
    Δγ::Float64
    function VonMisesState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end


matching_elem_type(::VonMises) = MechSolidElem
ip_state_type(::MechSolidElem, ::VonMises) = VonMisesState


function yield_func(matparams::VonMises, state::VonMisesState, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    σy  = matparams.σy
    H   = matparams.H
    εpa = state.εpa
    return √(3*j2d) - σy - H*εpa
end

function calcD(matparams::VonMises, state::VonMisesState)
    σy = matparams.σy
    H  = matparams.H
    #De = matparams.De
    De  = calcDe(matparams.E, matparams.ν, state.env.anaprops.stressmodel)

    if state.Δγ==0.0
        return De
    end

    j2d = J2D(state.σ)
    @assert j2d>0
    s  = dev(state.σ)
    su = s/norm(s)
    V  = √(3/2)*su # df/dσ
    Nu = su

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
end

function update_state(matparams::VonMises, state::VonMisesState, Δε::Array{Float64,1})
    σini = state.σ
    De   = calcDe(matparams.E, matparams.ν, state.env.anaprops.stressmodel)
    σtr  = state.σ + inner(De, Δε)
    ftr  = yield_func(matparams, state, σtr)

    if ftr < 1.e-8
        # elastic
        state.Δγ = 0.0
        state.σ  = σtr
    else
        # plastic
        K, G  = matparams.E/(3.0*(1.0-2.0*matparams.ν)), matparams.E/(2.0*(1.0+matparams.ν))
        H     = matparams.H
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        @assert √j2dtr - state.Δγ*√2*G > 0.0
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

function ip_state_vals(matparams::VonMises, state::VonMisesState)
    ndim  = state.env.ndim
    σ, ε  = state.σ, state.ε
    j1    = tr(σ)
    srj2d = √J2D(σ)

    D = stress_strain_dict(σ, ε, state.env.anaprops.stressmodel)
    D[:epa]   = state.εpa
    D[:j1]    = j1
    D[:srj2d] = srj2d

    return D
end
