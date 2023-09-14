# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Mazars

mutable struct MazarsState<:IpState
    env::ModelEnv
    σ::Vec6
    ε::Vec6
    φt::Float64
    φc::Float64
    φ::Float64  # damage
    ε̅max::Float64
    function MazarsState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.φ = 0.0
        this.φc = 0.0
        this.φt = 0.0
        this.ε̅max = 0.0
        return this
    end
end

mutable struct Mazars<:Material
    E ::Float64
    ν::Float64
    ε̅0::Float64
    At::Float64
    Bt::Float64
    Ac::Float64
    Bc::Float64
    ρ::Float64
    De::Mat6x6
    invDe::Mat6x6

    function Mazars(prms::Dict{Symbol,Float64})
        return Mazars(;prms...)
    end

    function Mazars(;E=NaN, nu=0.0, eps0=NaN, At=NaN, Bt=NaN, Ac=NaN, Bc=NaN, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert At>0.0
        @assert Ac>0.0
        @assert Bt>0.0
        @assert Bc>0.0
        @assert rho>=0.0
        @assert eps0>0

        this     = new(E, nu, eps0, At, Bt, Ac, Bc, rho)
        this.De  = calcDe(E, nu, "3d")
        this.invDe = inv(this.De)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{Mazars}, ::Type{MechSolid}, env::ModelEnv) = MazarsState


# Element types that work with this material
compat_elem_types(::Type{Mazars}) = (MechSolid,)


function calcD(mat::Mazars, state::MazarsState)
    # There is something wrong with the derivatives here

    # Equivalent strain scalar
    εp = eigvals(state.ε)
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    ε̅max = max(state.ε̅max, mat.ε̅0)

    if ε̅<ε̅max
        #@show "elastic"
        return (1.0 - state.φ)*mat.De
    else
        #@show "plastic"
        #return (1.0 - state.φ)*mat.De
        return mat.De

        # Principal stresses and principal directions
        σp, V = eigen(state.σ)
        σp = [ σp; zeros(3) ]
        # Eigen vectors
        p1 = V[:, 1]
        p2 = V[:, 2]
        p3 = V[:, 3]
        P  = [ matrix2Mandel(p1*p1'), matrix2Mandel(p2*p2'), matrix2Mandel(p3*p3') ]

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        Dei= inv(mat.De)
        εt = Dei*σt
        εc = Dei*σc
        εv = sum(pos.(εt)) + sum(neg.(εc))

        εv +=  1e-15 # avoid division by zero

        # Tensile and compression damage weights
        αt = clamp(sum(pos.(εt))/εv, 0.0, 1.0)
        αc = clamp(sum(neg.(εc))/εv, 0.0, 1.0)

        # Constitutive matrix calculation
        dφtdε̅ = (1.0-mat.At)*mat.ε̅0*ε̅^-2 + mat.At*mat.Bt*exp( -mat.Bt*(ε̅-mat.ε̅0) )
        dφcdε̅ = (1.0-mat.Ac)*mat.ε̅0*ε̅^-2 + mat.Ac*mat.Bc*exp( -mat.Bc*(ε̅-mat.ε̅0) )
        dε̅dε  = sum( pos(εp[i])*(1+sign(εp[i]))*P[i] for i in 1:3 ) / (2*ε̅)
        dφdε  = (αt*dφtdε̅ + αc*dφcdε̅) * dε̅dε

        #@show dφdε'*mat.De
        #D     = (1.0 - state.φ)*mat.De - (dφdε'*mat.De)'*state.ε'
        D     = (1.0 - state.φ)*mat.De - dφdε*(mat.De*state.ε)'
        return D
    end
end


function update_state!(mat::Mazars, state::MazarsState, Δε::AbstractArray)
    σini  = state.σ
    state.ε = state.ε + Δε

    E  = mat.E
    nu = mat.ν

    # Principal stresses tensor
    εp = eigvals(state.ε)

    # Equivalent strain scalar
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    state.ε̅max = max(state.ε̅max, mat.ε̅0)

    if ε̅ < state.ε̅max  # linear-elastic increment
        #@show "Elast"
        state.σ = (1.0 - state.φ)*mat.De*state.ε
    else # increment with damage: A previous elastic step is desired
        #@show "Plast"
        state.ε̅max = ε̅

        # Principal stresses and principal directions
        σp = eigvals(state.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        state.φt = 1.0 - (1-mat.At)*mat.ε̅0/ε̅ - mat.At/exp(mat.Bt*(ε̅-mat.ε̅0))
        state.φc = 1.0 - (1-mat.Ac)*mat.ε̅0/ε̅ - mat.Ac/exp(mat.Bc*(ε̅-mat.ε̅0))

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        εt = mat.invDe*σt
        εc = mat.invDe*σc

        #εv = sum(pos(εt[i]) + pos(εc[i]) for i in 1:3)
        εv = sum(pos.(εt)) + sum(neg.(εc))
        εv == 0.0 && (εv += 1e-15)

        # Tensile and compressive damage weights
        αt = clamp(sum(pos.(εt))/εv, 0.0, 1.0)
        αc = clamp(sum(neg.(εc))/εv, 0.0, 1.0)

        # Damage variable
        φ = αt*state.φt + αc*state.φc
        state.φ = clamp(φ, state.φ, 0.999)

        # Total stress and stress increment
        state.σ = (1.0 - state.φ)*mat.De*state.ε
    end

    Δσ    = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::Mazars, state::MazarsState)
    ndim  = state.env.ndim
    σ, ε  = state.σ, state.ε

    D = stress_strain_dict(σ, ε, state.env.ana.stressmodel)
    D[:dam]  = state.φ
    D[:damt] = state.φt
    D[:damc] = state.φc
    D[:eq]   = state.ε̅max

    return D
end
