# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidLinSeep

mutable struct ElasticSolidLinSeepState<:IpState
    env::ModelEnv
    σ::Array{Float64,1} # stress
    ε::Array{Float64,1} # strain
    V::Array{Float64,1} # fluid velocity
    D::Array{Float64,1} # distance traveled by the fluid
    uw::Float64         # pore pressure
    function ElasticSolidLinSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.V = zeros(env.ndim)
        this.D = zeros(env.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct ElasticSolidLinSeep<:Material
    E ::Float64 # Young's modulus
    nu::Float64 # Poisson ratio
    k ::Float64 # specific permeability
    γw::Float64 # specific weight of the fluid
    α ::Float64 # Biot's coefficient
    S ::Float64 # Storativity coefficient

    function ElasticSolidLinSeep(prms::Dict{Symbol,Float64})
        return  ElasticSolidLinSeep(;prms...)
    end

    function ElasticSolidLinSeep(;E=NaN, nu=NaN, k=NaN, kappa=NaN, gammaw=NaN, alpha=1.0, S=0.0, n=NaN, Ks=NaN, Kw=NaN, eta=NaN)

        if isnan(k)
            k = (kappa*gammaw)/eta # specific permeability = (intrinsic permeability * fluid specific weight)/viscosity
        end

        if isnan(S)
            S = (alpha - n)/Ks + n/Kw # S = (alpha - porosity)/(bulk module of the solid) + (porosity)/(bulk module of the fluid)
        end

        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu")
        k>0         || error("Invalid value for k: $k")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        0<alpha<=1.0|| error("Invalid value for alpha: $alpha")
        S>=0.0      || error("Invalid value for S: $S")

        this = new(E, nu, k, gammaw, alpha, S)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidLinSeep) = HMSolid

# Type of corresponding state structure
ip_state_type(mat::ElasticSolidLinSeep) = ElasticSolidLinSeepState

function set_state(state::ElasticSolidLinSeepState; sig=zeros(0), eps=zeros(0))
    sq2 = √2.0
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        state.σ .= sig.*mdl
    else
        if length(sig)!=0; error("ElasticSolidLinSeep: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        state.ε .= eps.*mdl
    else
        if length(eps)!=0; error("ElasticSolidLinSeep: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::ElasticSolidLinSeep, state::ElasticSolidLinSeepState)
    return calcDe(mat.E, mat.nu, state.env.modeltype) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidLinSeep, state::ElasticSolidLinSeepState) # Hydraulic conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidLinSeep, state::ElasticSolidLinSeepState, Δε::Array{Float64,1}, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcD(mat, state)
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ
    K = calcK(mat, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return Δσ, state.V
end

function ip_state_vals(mat::ElasticSolidLinSeep, state::ElasticSolidLinSeepState)
    D = stress_strain_dict(state.σ, state.ε, state.env.modeltype)

    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.env.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
