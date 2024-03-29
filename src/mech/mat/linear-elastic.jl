# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

LinearElastic_params = [
    FunInfo(:LinearElastic, "Linear-elastic material model"),
    KwArgInfo(:E, "Young modulus", cond=:(E>0.0)),
    KwArgInfo(:nu, "Poisson ratio", 0.0, cond=:(0.0<=nu<0.5)) # default value to easy the use of truss elements
]
@doc docstring(LinearElastic_params) LinearElastic

mutable struct LinearElastic<:Material
    E ::Float64
    ν::Float64

    function LinearElastic(; args...)
        args = checkargs(args, LinearElastic_params)
        return new(args.E, args.nu)
    end
end


mutable struct ElasticSolidState<:IpState
    env::ModelEnv
    σ::Vec6
    ε::Vec6

    function ElasticSolidState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


mutable struct ElasticPlaneStressState<:IpState
    env::ModelEnv
    σ::Vec6
    ε::Vec6

    function ElasticPlaneStressState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


mutable struct ElasticBarState<:IpState
    env::ModelEnv
    σ::Float64
    ε::Float64
    function ElasticBarState(env::ModelEnv)
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end


compat_state_type(::Type{LinearElastic}, ::Type{MechSolid}, env::ModelEnv)  = env.ana.stressmodel==:planestress ? ElasticPlaneStressState : ElasticSolidState
compat_state_type(::Type{LinearElastic}, ::Type{MechShell}, env::ModelEnv)  = ElasticPlaneStressState
compat_state_type(::Type{LinearElastic}, ::Type{MechBeam}, env::ModelEnv)   = ElasticBeamState
compat_state_type(::Type{LinearElastic}, ::Type{MechBar}, env::ModelEnv)    = ElasticBarState
compat_state_type(::Type{LinearElastic}, ::Type{MechEmbBar}, env::ModelEnv) = ElasticBarState


function calcDe(E::Real, ν::Real, stressmodel::Symbol=:d3)
    if stressmodel==:planestress
        c = E/(1-ν^2)
        return @SArray [
            c     c*ν   0.0   0.0        0.0        0.0
            c*ν   c     0.0   0.0        0.0        0.0
            0.0   0.0   0.0   0.0        0.0        0.0
            0.0   0.0   0.0   c*(1.0-ν)  0.0        0.0
            0.0   0.0   0.0   0.0        c*(1.0-ν)  0.0
            0.0   0.0   0.0   0.0        0.0        c*(1.0-ν) ]
        ezz = -ν/E*(sxx+syy)
    else
        c = E/((1+ν)*(1-2*ν))
        return @SArray [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ]
    end
end


# LinearElastic model for 3D and 2D bulk elements under plain-strain state


function calcD(mat::LinearElastic, state::ElasticSolidState)
    return calcDe(mat.E, mat.ν)
end


function update_state!(mat::LinearElastic, state::ElasticSolidState, dε::AbstractArray)
    De = calcDe(mat.E, mat.ν)

    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticSolidState)
    return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
end


# LinearElastic model for 2D bulk elements under plane-stress state and shell elements


function calcD(mat::LinearElastic, state::ElasticPlaneStressState)
    return calcDe(mat.E, mat.ν, :planestress)
end


function update_state!(mat::LinearElastic, state::ElasticPlaneStressState, dε::AbstractArray)
    De = calcDe(mat.E, mat.ν, :planestress)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticPlaneStressState)
    return stress_strain_dict(state.σ, state.ε, :planestress)
end


# LinearElastic for beam elements

mutable struct ElasticBeamState<:IpState
    env::ModelEnv
    σ::Vec3
    ε::Vec3

    function ElasticBeamState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(Vec3)
        this.ε = zeros(Vec3)
        return this
    end
end


function calcD(mat::LinearElastic, state::ElasticBeamState)
    E, ν = mat.E, mat.ν
    
    return @SMatrix [ 
        E  0.0  0.0  
        0.0  E/(1+ν)  0.0  
        0.0  0.0  E/(1+ν)
    ]
end


function update_state!(mat::LinearElastic, state::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticBeamState)
    vals = OrderedDict{Symbol,Float64}(
      :sX  => state.σ[1],
      :eX  => state.ε[1],
      :sXY => state.σ[2]/SR2
    )
    if state.env.ndim==3
        vals[:sXZ] = state.σ[2]/SR2
        vals[:sXY] = state.σ[3]/SR2
    end
    return vals
end


# LinearElastic model for bar elements

function calcD(mat::LinearElastic, ips::ElasticBarState)
    return mat.E
end


function update_state!(mat::LinearElastic, state::ElasticBarState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticBarState)
    return OrderedDict(
      :sX => state.σ,
      :eX => state.ε,
      )
end