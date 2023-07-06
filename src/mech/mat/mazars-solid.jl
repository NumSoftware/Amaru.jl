# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Mazars

mutable struct MazarsState<:IpState
    env::ModelEnv
    σ::Tensor2
    ε::Tensor2
    φt::Float64
    φc::Float64
    φ::Float64  # damage
    ε̅max::Float64
    function MazarsState(env::ModelEnv=ModelEnv())
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

mutable struct Mazars<:MatParams
    E ::Float64
    nu::Float64
    ε̅0::Float64
    At::Float64
    Bt::Float64
    Ac::Float64
    Bc::Float64
    ρ::Float64
    De::Tensor4
    invDe::Tensor4

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
        this.invDe  = inv(this.De)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::Mazars) = MechSolidElem

# Type of corresponding state structure
ip_state_type(::MechSolidElem, ::Mazars) = MazarsState


function calcD(matparams::Mazars, state::MazarsState)
    # There is something wrong with the derivatives here

    # Equivalent strain scalar
    εp = eigvals(state.ε)
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    ε̅max = max(state.ε̅max, matparams.ε̅0)

    if ε̅<ε̅max
        #@show "elastic"
        return (1.0 - state.φ)*matparams.De
    else
        #@show "plastic"
        #return (1.0 - state.φ)*matparams.De
        return matparams.De

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
        Dei= inv(matparams.De)
        εt = Dei*σt
        εc = Dei*σc
        εv = sum(pos.(εt)) + sum(neg.(εc))

        εv +=  1e-15 # avoid division by zero

        # Tensile and compression damage weights
        αt = clamp(sum(pos.(εt))/εv, 0.0, 1.0)
        αc = clamp(sum(neg.(εc))/εv, 0.0, 1.0)

        # Constitutive matrix calculation
        dφtdε̅ = (1.0-matparams.At)*matparams.ε̅0*ε̅^-2 + matparams.At*matparams.Bt*exp( -matparams.Bt*(ε̅-matparams.ε̅0) )
        dφcdε̅ = (1.0-matparams.Ac)*matparams.ε̅0*ε̅^-2 + matparams.Ac*matparams.Bc*exp( -matparams.Bc*(ε̅-matparams.ε̅0) )
        dε̅dε  = sum( pos(εp[i])*(1+sign(εp[i]))*P[i] for i in 1:3 ) / (2*ε̅)
        dφdε  = (αt*dφtdε̅ + αc*dφcdε̅) * dε̅dε

        #@show dφdε'*matparams.De
        #D     = (1.0 - state.φ)*matparams.De - (dφdε'*matparams.De)'*state.ε'
        D     = (1.0 - state.φ)*matparams.De - dφdε*(matparams.De*state.ε)'
        return D
    end
end


function update_state(matparams::Mazars, state::MazarsState, Δε::Array{Float64,1})
    σini  = state.σ
    state.ε = state.ε + Δε

    E  = matparams.E
    nu = matparams.nu

    # Principal stresses tensor
    εp = eigvals(state.ε)

    # Equivalent strain scalar
    ε̅ = norm(pos.(εp))
    ε̅ == 0.0 && (ε̅ += 1e-15)
    state.ε̅max = max(state.ε̅max, matparams.ε̅0)

    if ε̅ < state.ε̅max  # linear-elastic increment
        #@show "Elast"
        state.σ = (1.0 - state.φ)*matparams.De*state.ε
    else # increment with damage: A previous elastic step is desired
        #@show "Plast"
        state.ε̅max = ε̅

        # Principal stresses and principal directions
        σp = eigvals(state.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        state.φt = 1.0 - (1-matparams.At)*matparams.ε̅0/ε̅ - matparams.At/exp(matparams.Bt*(ε̅-matparams.ε̅0))
        state.φc = 1.0 - (1-matparams.Ac)*matparams.ε̅0/ε̅ - matparams.Ac/exp(matparams.Bc*(ε̅-matparams.ε̅0))

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        εt = matparams.invDe*σt
        εc = matparams.invDe*σc

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
        state.σ = (1.0 - state.φ)*matparams.De*state.ε
    end

    Δσ    = state.σ - σini
    return Δσ, success()
end

function ip_state_vals(matparams::Mazars, state::MazarsState)
    ndim  = state.env.ndim
    σ, ε  = state.σ, state.ε

    D = stress_strain_dict(σ, ε, state.env.anaprops.stressmodel)
    D[:dam]  = state.φ
    D[:damt] = state.φt
    D[:damc] = state.φc
    D[:eq]   = state.ε̅max

    return D
end
