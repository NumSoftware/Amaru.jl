# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

"""
    LinearElastic

A type for linear elastic materials.

# Fields

$(TYPEDFIELDS)
"""
mutable struct LinearElastic<:Material
    "Young Modulus"
    E ::Float64
    "Poisson ratio"
    ν::Float64

    function LinearElastic(; args...)
        args = checkargs(args, arg_rules(LinearElastic))       
        return new(args.E, args.nu)
    end
end

arg_rules(::Type{LinearElastic}) = 
[
    @arginfo E E>0.0 "Young modulus"
    @arginfo nu=0.0 0.0<=nu<0.5 "Poisson ratio"
]


"""
    ElasticSolidState

A type for the state data of a `LinearElastic` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticSolidState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Vec6
    "Strain tensor"
    ε::Vec6

    function ElasticSolidState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


"""
    ElasticSolidState

A type for the state data of a `LinearElastic` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticPlaneStressState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::SArray
    "Strain tensor"
    ε::SArray

    function ElasticPlaneStressState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end

# compat_state_type(::Type{LinearElastic}, ::Type{MechSolid}, env::ModelEnv) = env.ana.stressmodel=="plane-stress" ? ElasticSolidPlaneStressState : ElasticSolidState

compat_state_type(::Type{LinearElastic}, ::Type{MechSolid}, env::ModelEnv)  = ElasticSolidState
compat_state_type(::Type{LinearElastic}, ::Type{MechShell}, env::ModelEnv)  = ElasticSolidState


# Element types that work with this material
# compat_elem_types(::Type{LinearElastic}) = (MechSolid, MechShell)


function calcDe(E::Number, ν::Number, stressmodel::String)
    if stressmodel=="plane-stress"
        c = E/(1-ν^2)
        return @SArray [
            c     c*ν   0.0   0.0        0.0        0.0
            c*ν   c     0.0   0.0        0.0        0.0
            0.0   0.0   1.0   0.0        0.0        0.0
            0.0   0.0   0.0   c*(1.0-ν)  0.0        0.0
            0.0   0.0   0.0   0.0        c*(1.0-ν)  0.0
            0.0   0.0   0.0   0.0        0.0        c*(1.0-ν) ]
        ezz = -ν/E*(sxx+syy)
    # elseif stressmodel=="shell"
    #     c = E/(1.0-ν^2)
    #     g = E/(1+ν)
    #     return @SArray [
    #         c    c*ν   0.0  0.0    0.0    0.0
    #         c*ν  c     0.0  0.0    0.0    0.0
    #         0.0  0.0   0.0  0.0    0.0    0.0
    #         0.0  0.0   0.0  5/6*g  0.0    0.0
    #         0.0  0.0   0.0  0.0    5/6*g  0.0
    #         0.0  0.0   0.0  0.0    0.0    g ]
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


function calcD(mat::LinearElastic, state::ElasticSolidState, stressmodel::String=state.env.ana.stressmodel)
    return calcDe(mat.E, mat.ν, stressmodel)
end


function update_state!(mat::LinearElastic, state::ElasticSolidState, dε::AbstractArray, stressmodel::String=state.env.ana.stressmodel)
    De = calcDe(mat.E, mat.ν, stressmodel)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticSolidState, stressmodel::String=state.env.ana.stressmodel)
    return stress_strain_dict(state.σ, state.ε, stressmodel)
end



"""
    ElasticRodState

A type for the state data of a `ElasticRod` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticRodState<:IpState
    "environment information"
    env::ModelEnv
    "Axial stress"
    σ::Float64
    "Axial strain"
    ε::Float64
    function ElasticRodState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end


compat_state_type(::Type{LinearElastic}, ::Type{MechRod}, env::ModelEnv) = ElasticRodState
compat_state_type(::Type{LinearElastic}, ::Type{MechEmbRod}, env::ModelEnv) = ElasticRodState


function calcD(mat::LinearElastic, ips::ElasticRodState)
    return mat.E
end


function update_state!(mat::LinearElastic, state::ElasticRodState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticRodState)
    return OrderedDict(
      :sX => state.σ,
      :eX => state.ε,
      )
end




mutable struct ElasticBeamState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticBeamState(env::ModelEnv=ModelEnv())
        this = new(env)
        nstr = env.ndim==2 ? 2 : 3
        this.σ = zeros(nstr)
        this.ε = zeros(nstr)
        return this
    end
end


compat_state_type(::Type{LinearElastic}, ::Type{MechBeam}, env::ModelEnv) = ElasticBeamState
# compat_state_type(::Type{LinearElastic}, ::Type{MechShell}, env::ModelEnv) = ElasticShellState
# compat_state_type(::Type{LinearElastic}, ::Type{MechBeam}, env::ModelEnv)  = ElasticBeamState

# Type of corresponding state structure
# compat_state_type(::Type{LinearElastic}) = ElasticBeamState

# Element types that work with this material
# compat_elem_types(::Type{LinearElastic}) = (MechBeam,)


function calcD(mat::LinearElastic, state::ElasticBeamState)
    E = mat.E
    ν = mat.ν
    G = E/2/(1+ν)

    if state.env.ndim==2
        return [ E      0.0
                 0.0  5/6*2*G ]
    else
        return [ E   0.0    0.0  
                0.0  5/6*2*G  0.0  
                0.0  0.0    5/6*2*G ]
    end
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



# mutable struct ElasticShellState<:IpState
#     "Environment information"
#     env::ModelEnv
#     "Stress tensor"
#     σ::Vec6 # in the local system
#     "Strain tensor"
#     ε::Vec6 # in the local system

#     function ElasticShellState(env::ModelEnv=ModelEnv())
#         this = new(env)
#         this.σ = zeros(6)
#         this.ε = zeros(6)
#         return this
#     end
# end

# function calcD(mat::LinearElastic, state::ElasticShellState)
#     E = mat.E
#     ν = mat.ν
#     c = E/(1-ν^2)
#     g = E/(1+ν)
#     return [
#         c    c*ν   0.0  0.0    0.0    0.0
#         c*ν  c     0.0  0.0    0.0    0.0
#         0.0  0.0   0.0  0.0    0.0    0.0
#         0.0  0.0   0.0  5/6*g  0.0    0.0
#         0.0  0.0   0.0  0.0    5/6*g  0.0
#         0.0  0.0   0.0  0.0    0.0    g ]
#     # ezz = -ν/E*(sxx+syy)
# end


# function update_state!(mat::LinearElastic, state::ElasticShellState, dε::Array{Float64,1})
#     D = calcD(mat, state)
#     dσ = D*dε
#     state.ε += dε
#     state.σ += dσ
#     return dσ, success()
# end


# function ip_state_vals(mat::LinearElastic, state::ElasticShellState)
#     return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
# end

