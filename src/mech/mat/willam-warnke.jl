# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export WillamWarnke

WillamWarnke_params = [
    FunInfo(:WillamWarnke, "Linear-elastic model with Willam-Warnke yield criterion and linear/non-linear hardening"),
    KwArgInfo(:E, "Young's modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson's ratio", cond=:(nu>=0)),
    KwArgInfo(:ft, "Tensile strength", cond=:(ft>0)),
    KwArgInfo(:fc, "Compressive strength", cond=:(fc<0)),
    KwArgInfo(:fb, "Biaxial compressive strength", cond=:(fb<0)),
    KwArgInfo(:psi, "Dilatancy angle in degrees", nothing, cond=:(0<psi<90)),
    KwArgInfo(:H, "Hardening modulus", nothing, cond=:(H>=0)),
    KwArgInfo((:xi_fun, :ξfun), "Hardening curve in terms of ξ0", nothing),
]
@doc docstring(WillamWarnke_params) WillamWarnke

mutable struct WillamWarnke<:Material
    E::Float64
    ν::Float64
    e::Float64  # eccentricity of the yield surface
    ξ0::Float64 # distance from the origin to the apex
    α::Float64  # slope of the yield surface
    ᾱ::Float64  # slope of the plastic potential surface
    H::Float64  # hardening modulus
    ξfun::Union{Nothing,PathFunction}

    function WillamWarnke(; kwargs...)
        args = checkargs(kwargs, WillamWarnke_params)

        fc, fb, ft = abs(args.fc), abs(args.fb), args.ft
        xi_fun, H = args.xi_fun, args.H

        if isnothing(xi_fun)
            isnothing(args.H) && error("WillamWarnke: H or xi_fun must be provided")
        else
            H = 0.0
        end
        
        e  = (fc*fb - fc*ft + 3*fb*ft)/(2*fc*fb + ft*fc)
        @assert 0<e<=1

        ξ0 = √3*fb*ft/(fb-ft) # ok
        ρ0 = √2*fc*ξ0/(fc + √3*ξ0) # ok
        α = ρ0/ξ0

        ψ = args.psi
        if isnothing(ψ)
            ᾱ = e*α
        else
            ᾱ = tand(ψ)
        end

        this = new(args.E, args.nu, e, ξ0, α, ᾱ, H, args.xi_fun)
        return this
    end
end


mutable struct WillamWarnkeState<:IpState
    ctx::Context
    σ::Vec6
    ε::Vec6
    εpa::Float64
    Δλ::Float64
    function WillamWarnkeState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εpa = 0.0
        this.Δλ  = 0.0
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{WillamWarnke}, ::Type{MechSolid}, ctx::Context) = ctx.stressmodel!=:planestress ? WillamWarnkeState : error("WillamWarnke: This model is not compatible with planestress")


function calc_r(mat::WillamWarnke, σ::Vec6)
    e = mat.e
    j2, j3 = J2(σ), J3(σ)

    norm_s = √(2*j2)
    det_s  = j3
    θ      = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden
    return r
end


@inline function calc_ξmax(mat::WillamWarnke, εpa::Float64)
    if mat.ξfun !== nothing
        ξmax = mat.ξfun(εpa)
    else
        ξmax = mat.ξ0 + mat.H*εpa
    end
end


@inline function calc_deriv_ξmax_εpa(mat::WillamWarnke, εpa::Float64)
    if mat.ξfun !== nothing
        dξdεpa = derive(mat.ξfun, εpa)
    else
        dξdεpa = mat.H
    end
    return dξdεpa
end


function yield_func(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray, εpa::Float64)
    α = mat.α
    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)
    r = calc_r(mat, σ)
    ξmax = calc_ξmax(mat, εpa)

    return ρ - r*α*(ξmax - ξ)
end


function yield_derivs(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray, εpa::Float64)
    
    e, α = mat.e, mat.α
    i1, j2 = tr(σ), J2(σ)

    ρ = √(2*j2)
    ξ = i1/√3

    # deviatoric derivatives
    s      = dev(σ)
    det_s  = J3(σ)
    adj_s  = det_s*inv(s)
    norm_s = ρ

    # θ and derivatives
    θ       = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    rnum    = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden    = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r       = rnum/rden
    drnumdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ)
    drdendθ = 4*sin(2*θ)*(e^2-1)
    drdθ    = (drnumdθ*rden - rnum*drdendθ)/rden^2

    if 1-abs(cos(3*θ)) > 1e-6 # condition to avoid division by zero
        dθds = -√6*(adj_s/ρ^3 - 3*s*det_s/ρ^5)/√abs(1 - 54*det_s^2/ρ^6)
    else
        dθds = 0.0*I2
    end

    # f derivatives
    dfdρ = 1.0
    dfdξ = α*r
    ξmax = calc_ξmax(mat, εpa)
    dfdr = α*(ξmax - ξ)
    dfdθ = dfdr*drdθ

    dρdσ = s/norm(s)
    dξdσ = √3/3*I2
    dsdσ = Psd
    dθdσ = dsdσ*dθds
    
    dfdσ      = dfdρ*dρdσ + dfdξ*dξdσ + dfdθ*dθdσ
    dfdξmax   = -α*r
    dξmaxdεpa = calc_deriv_ξmax_εpa(mat, εpa)
    dfdεpa    = dfdξmax*dξmaxdεpa

    return dfdσ, dfdεpa
end


function potential_derivs(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray)
    # g(σ) = ρ - ᾱ*(ξ0-ξ)  Drucker-Prager like surface

    ᾱ = mat.ᾱ
    s = dev(σ)

    # g derivatives
    # dgdρ = 1.0
    # dgdξ = ᾱ
    # dρdσ = s/norm(s)
    # dξdσ = √3/3*I2

    dgdσ = s/norm(s) + ᾱ*√3/3*I2
    return dgdσ

end


function calcD(mat::WillamWarnke, state::WillamWarnkeState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stressmodel)

    if state.Δλ==0.0
        return De
    end

    j2 = J2(state.σ)

    if j2 != 0.0
        dfdσ, dfdεpa = yield_derivs(mat, state, state.σ, state.εpa)
        dgdσ = potential_derivs(mat, state, state.σ)
    else # apex
        dfdσ = √3/3*I2
        dgdσ = dfdσ
    end

    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεpa*norm(dgdσ))
    return Dep

end


function calc_σ_εpa_Δλ(mat::WillamWarnke, state::WillamWarnkeState, σtr::Vec6)

    α  = mat.α
    ᾱ  = mat.ᾱ

    E, ν = mat.E, mat.ν
    K, G  = E/(3.0*(1.0-2.0*ν)), E/(2.0*(1.0+ν))

    maxits = 20
    tol = 1e-3

    # ftr = yield_func(mat, state, σtr, state.εpa)

    # @show f
    str   = dev(σtr)
    n_str = norm(str)
    ξtr   = tr(σtr)/√3
    ρtr   = n_str
    r     = calc_r(mat, σtr)
    ξmax  = calc_ξmax(mat, state.εpa)
    
    # iterative process since θ is unknown at step n+1
    for i in 1:maxits

        numΔλ = ρtr - α*r*(ξmax - ξtr)
        denΔλ = 2*G + 3*α*ᾱ*K*r

        Δλ = numΔλ/denΔλ

        # @assert Δλ >= 0.0
        # todo: check if Δλ is negative and return to the apex
        
        # estimative at n+1
        s   = (1 - 2*G*Δλ/n_str)*str
        ξ   = ξtr - 3*ᾱ*K*Δλ
        σ   = ξ/√3*I2 + s
        εpa = state.εpa + Δλ*√(1 + ᾱ^2)
        ρ   = ρtr - 2*G*Δλ
        r   = calc_r(mat, σ)
        ξmax = calc_ξmax(mat, εpa)

        # yield function at n+1
        f = ρ - r*α*(ξmax - ξ)

        if abs(f) < tol
            @assert Δλ >= 0.0

            return σ, εpa, Δλ, success()
        end

    end


    return state.σ, 0.0, 0.0, failure("WillamWarnke: maximum iterations reached")

end


function calc_σ_εpa_Δλ_bis(mat::WillamWarnke, state::WillamWarnkeState, σtr::Vec6)
    α  = mat.α
    ᾱ  = mat.ᾱ

    E, ν = mat.E, mat.ν
    K, G  = E/(3.0*(1.0-2.0*ν)), E/(2.0*(1.0+ν))

    # De = calcDe(mat.E, mat.ν, state.ctx.stressmodel)

    ξtr = tr(σtr)/√3
    str = dev(σtr)
    n_str = norm(str)
    r     = calc_r(mat, σtr)
    ξmax  = calc_ξmax(mat, state.εpa)

    # estimative of Δλ
    # dgdσ = potential_derivs(mat, state, σtr)
    # Δλ0  = norm(σtr-state.σ)/norm(De*dgdσ)

    # str   = dev(σtr)
    ρtr   = n_str
    numΔλ = ρtr - α*r*(ξmax - ξtr)
    denΔλ = 2*G + 3*α*ᾱ*K*r
    Δλ0 = numΔλ/denΔλ

    # function of Δλ
    ff(Δλ)  = begin
        # quantities at n+1
        ρ = ρtr - 2*G*Δλ
        if ρ>0
            ξ = ξtr - 3*ᾱ*K*Δλ
            s = (1 - 2*G*Δλ/n_str)*str  # todo: avoid to compute s
            σ = ξ/√3*I2 + s
        else
            ξ = calc_ξmax(mat, state.εpa)
            σ = ξ/√3*I2
        end
        
        εpa = state.εpa + Δλ*√(1 + ᾱ^2)
        return yield_func(mat, state, σ, εpa)
    end

    a, b, status = findrootinterval(ff, 0.0, Δλ0)
    failed(status) && return state.σ, 0.0, 0.0, status

    Δλ, status = findroot(ff, a, b, ftol=1e-3, method=:bisection)
    failed(status) && return state.σ, 0.0, 0.0, status
    @assert Δλ >= 0.0

    σ   = σtr - 2*G*Δλ*str/n_str - √3*K*Δλ*ᾱ*I2
    εpa = state.εpa + Δλ*√(1 + ᾱ^2)
    
    return σ, εpa, Δλ, success()

end


function update_state!(mat::WillamWarnke, state::WillamWarnkeState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stressmodel)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, state, σtr, state.εpa)

    if ftr < 1e-8
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        σ, εpa, Δλ, status = calc_σ_εpa_Δλ(mat, state, σtr)
        # σ, εpa, Δλ, status = calc_σ_εpa_Δλ_bis(mat, state, σtr)
        failed(status) && return state.σ, status
        
        state.σ, state.εpa, state.Δλ = σ, εpa, Δλ
    end

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::WillamWarnke, state::WillamWarnkeState)
    σ, ε  = state.σ, state.ε
    ρ = √(2*J2(σ))
    ξ = tr(σ)/√3

    D       = stress_strain_dict(σ, ε, state.ctx.stressmodel)
    D[:ep]  = state.εpa
    D[:xi]  = ξ
    D[:rho] = ρ

    return D
end
