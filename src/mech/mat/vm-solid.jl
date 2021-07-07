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
    DruckerPragerIpState

A type for the state data of a `DruckerPrager` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct VonMisesIpState<:IpState
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
    function VonMisesIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end

matching_elem_type(::VonMises) = MechSolid
ip_state_type(mat::VonMises) = VonMisesIpState

function yield_func(mat::VonMises, ipd::VonMisesIpState, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    σy  = mat.σy
    H   = mat.H
    εpa = ipd.εpa
    return √(3*j2d) - σy - H*εpa
end

function calcD(mat::VonMises, ipd::VonMisesIpState)
    σy = mat.σy
    H  = mat.H
    #De = mat.De
    De  = calcDe(mat.E, mat.ν, ipd.env.modeltype)

    if ipd.Δγ==0.0
        return De
    end

    j2d = J2D(ipd.σ)
    @assert j2d>0
    s  = dev(ipd.σ)
    su = s/norm(s)
    V  = √(3/2)*su # df/dσ
    Nu = su

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
end

function stress_update(mat::VonMises, ipd::VonMisesIpState, Δε::Array{Float64,1})
    σini = ipd.σ
    De   = calcDe(mat.E, mat.ν, ipd.env.modeltype)
    σtr  = ipd.σ + inner(De, Δε)
    ftr  = yield_func(mat, ipd, σtr)

    if ftr < 1.e-8
        # elastic
        ipd.Δγ = 0.0
        ipd.σ  = σtr
    else
        # plastic
        K, G  = mat.E/(3.0*(1.0-2.0*mat.ν)), mat.E/(2.0*(1.0+mat.ν))
        H     = mat.H
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        @assert √j2dtr - ipd.Δγ*√2*G > 0.0
        ipd.Δγ = ftr/(√6*G + H)
        j1     = j1tr
        m      = 1.0 - ipd.Δγ*√2*G/√j2dtr
        ipd.σ  = m*dev(σtr) + j1/3.0*tI

        ipd.εpa += ipd.Δγ

    end

    ipd.ε += Δε
    Δσ     = ipd.σ - σini
    return Δσ, success()
end

function ip_state_vals(mat::VonMises, ipd::VonMisesIpState)
    ndim  = ipd.env.ndim
    σ, ε  = ipd.σ, ipd.ε
    j1    = tr(σ)
    srj2d = √J2D(σ)

    D = stress_strain_dict(σ, ε, ipd.env.modeltype)
    D[:epa]   = ipd.εpa
    D[:j1]    = j1
    D[:srj2d] = srj2d

    return D
end
