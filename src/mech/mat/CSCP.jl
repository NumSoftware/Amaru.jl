# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CSCP

CSCP_params = [
    FunInfo(:CSCP, "Closed surface concrete plasticity"),
    KwArgInfo(:E, "Young's modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson's ratio", cond=:(nu>=0)),
    # KwArgInfo(:e, "Excentricity", nothing, cond=:(1>=e>0.50)),
    KwArgInfo(:alpha, "Curvature coefficient", 1.5, cond=:(0.4<alpha<=1.0)),
    KwArgInfo(:beta, "Factor to get fb from fc_max", 1.15, cond=:(1<=beta<=1.5)),
    KwArgInfo(:fc, "Compressive strength OR plastic curve for fc as a function of εc", nothing),
    KwArgInfo(:ft, "Tensile strength OR plastic curve for ft as a function of w", nothing),
    KwArgInfo(:pmin, "Elastic limit for p under compression OR curve for p as a function of εv", nothing),
    KwArgInfo(:GF, "Tensile fracture energy (only used if ft is a number)", nothing, cond=:(GF>0)),
]


mutable struct CSCP<:Material
    E::Float64
    ν::Float64
    e::Float64
    α::Float64
    ft_fun::Union{Nothing,PathFunction}
    fc_fun::Union{Nothing,PathFunction}
    p_fun::Union{Nothing,PathFunction}
    fc_peak::Float64
    fb::Float64

    function CSCP(; kwargs...)
        args = checkargs(kwargs, CSCP_params)

        # Elastic constants
        E = args.E
        ν = args.nu
        K = E/(3*(1-2*ν))

        # Input functions
        ft_fun = args.ft isa PathFunction ? args.ft : PathFunction(:M, 0.0, args.ft)
        fc_fun = args.fc isa PathFunction ? args.fc : PathFunction(:M, 0.0, args.fc)
        p_fun  = args.pmin isa PathFunction ? args.pmin : PathFunction(:M, 0.0, args.pmin)

        # Fix functions to be in terms of plastic strains
        fc_fun = map((x, y) -> (abs(x-y/E), y), fc_fun)
        p_fun  = map((x, y) -> (abs(x-y/K), y), p_fun)

        GF = args.GF
        ft = args.ft
        if !isnothing(GF) && ft isa Number
            ft_fun = PathFunction(:M, 0.0, ft)
            wc     = GF/(0.1947*ft)
            W      = range(wc/10, wc, 10)
            for w in W
                z = (1 + 27*(w/wc)^3)*exp(-6.93*w/wc) - 28*(w/wc)*exp(-6.93)
                append!(ft_fun, :L, w, ft*max(z,0.0))
            end
        end

        # fc peak
        εcp_vals   = range(limits(fc_fun)..., 100)
        fc_vals    = fc_fun.(εcp_vals)
        fc_peak, _ = findmin(fc_vals) # fc peak value

        α = args.alpha
        β = args.beta

        # value of exentricity to match fb in a biaxial trajectory, assuming the state when ξb=0 
        e = β/(2*β)^α # forgot how to demonstate this

        fc_peak = abs(fc_peak)
        fb = β*fc_peak

        this = new(args.E, args.nu, e, α, ft_fun, fc_fun, p_fun, fc_peak, fb)
        return this
    end
end


mutable struct CSCPState<:IpState
    ctx::Context
    σ  ::Vec6
    ε  ::Vec6
    εcp::Float64
    εtp::Float64
    Δλ ::Float64
    h  ::Float64
    function CSCPState(ctx::Context)
        this     = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εcp = 0.0 # plastic strain in compression
        this.εtp = 0.0 # plastic strain in tension
        this.Δλ  = 0.0 # increment of plastic multiplier
        this.h   = 0.0 
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{CSCP}, ::Type{MechSolid}, ctx::Context) = ctx.stressmodel!=:planestress ? CSCPState : error("CSCP: This model is not compatible with planestress")


function calc_θ(::CSCP, σ::Vec6)
    j2 = J2(σ)
    if j2==0.0
        θ = 0.0
    else
        norm_s = √(2*j2)
        det_s  = J3(σ)
        θ      = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    end
    return θ
end


function calc_rθ(mat::CSCP, σ::Vec6)
    e = mat.e
    θ = calc_θ(mat, σ)

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end

function calc_rξ(mat::CSCP, ξb::Float64, ξ::Float64)
    α  = mat.α
    fc_peak = mat.fc_peak

    return spow((ξb-ξ)/fc_peak, α)
end

function calc_rc(mat::CSCP, ξa::Float64, ξ::Float64)
    ξc = 2*mat.fb/√3
    ξ>=ξc && return 1.0
    ξ<ξa  && return 0.0
    return √(1 - ((ξc-ξ)/(ξc-ξa))^2)
end


function calc_ξa_ξb_κ(mat::CSCP, state::CSCPState, εcp::Float64, εtp::Float64)
    α  = mat.α
    e  = mat.e
    fc_peak = mat.fc_peak

    w = εtp*state.h
    ft = mat.ft_fun(w)
    fc = mat.fc_fun(εcp)

    ξa = √3*mat.p_fun(√3*εcp)    # ξ = √3p ; plastic volumetric strain εvp = √3*εcp in isotropic compression
    @assert ξa<0
    @assert ξa<fc/√3

    Ω  = (-ft/(fc*e))^(1/α)
    ξb = 1/√3*(fc*Ω - ft)/(Ω-1)
    ξb<0 && @show ξb
    @assert ξb>=0

    κ  = -√(2/3)*fc*((ξb-fc/√3)/fc_peak)^-α
    @assert κ>0

    return return ξa, ξb, κ
end


function yield_func(mat::CSCP, state::CSCPState, σ::AbstractArray, εcp::Float64, εtp::Float64)
    # f(σ) = ρ - rθ⋅rc⋅rξ⋅κ

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εcp, εtp)
    rθ = calc_rθ(mat, σ)
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)

    return ρ - rθ*rc*rξ*κ
end


function yield_derivs(mat::CSCP, state::CSCPState, σ::AbstractArray, εcp::Float64, εtp::Float64)
    e = mat.e
    α = mat.α
    fc_peak = mat.fc_peak

    i1, j2 = tr(σ), J2(σ)

    ρ = √(2*j2)
    ξ = i1/√3

    # deviatoric derivatives
    s      = dev(σ) 
    det_s  = J3(σ)
    adj_s  = det_s*inv(s)
    norm_s = ρ

    # θ and derivatives
    θ        = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    rnum     = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden     = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    rθ       = rnum/rden
    drθnumdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ)
    drθdendθ = 4*sin(2*θ)*(e^2-1)
    drθdθ    = (drθnumdθ*rden - rnum*drθdendθ)/rden^2

    if 1-abs(cos(3*θ)) > 1e-6 # condition to avoid division by zero
        dθds = -√6*(adj_s/ρ^3 - 3*s*det_s/ρ^5)/√abs(1 - 54*det_s^2/ρ^6)
    else
        dθds = 0.0*I2
    end

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εcp, εtp)

    ξc = 2*mat.fb/√3
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)

    # f derivative w.r.t. σ:
    dfdρ  = 1
    dfdrc = -rθ*rξ*κ
    dfdrξ = -rθ*rc*κ
    drcdξ = ξa<ξ<ξc ? (ξc-ξ)/(ξc-ξa)^2/√(1-((ξc-ξ)/(ξc-ξa))^2) : 0.0
    drξdξ = ξb-ξ!=0.0 ? -α/fc_peak * abs((ξb-ξ)/fc_peak)^(α-1) : 0.0
    dfdξ  = dfdrc*drcdξ + dfdrξ*drξdξ
    dfdrθ = -rc*rξ*κ
    dfdθ  = dfdrθ*drθdθ

    dρdσ = s/norm(s)
    dξdσ = √3/3*I2
    dsdσ = Psd
    dθdσ = dsdσ*dθds

    if ρ==0 # apex
        dfdσ = √3/3*I2
    else
        dfdσ = dfdρ*dρdσ + dfdξ*dξdσ + dfdθ*dθdσ
    end
    
    f_εcp  = εcp -> yield_func(mat, state, σ, εcp, εtp)
    dfdεcp = derive(f_εcp, εcp)
    
    f_εtp  = εtp -> yield_func(mat, state, σ, εcp, εtp)
    dfdεtp = derive(f_εtp, εtp)

    return dfdσ, dfdεcp, dfdεtp
end


function potential_derivs(mat::CSCP, state::CSCPState, σ::AbstractArray, εcp::Float64, εtp::Float64)
    # f(σ) = ρ - e⋅rc⋅rξ⋅κ

    e  = mat.e
    α  = mat.α
    fc_peak = mat.fc_peak

    i1 = tr(σ)
    ξ  = i1/√3
    s  = dev(σ)
    ρ  = norm(s)
    ρ == 0 && return  √3/3*I2

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εcp, εtp)

    ξc = 2*mat.fb/√3
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)
    
    dgdrc = -e*rξ*κ
    dgdrξ = -e*rc*κ
    drcdξ = ξa<ξ<ξc ? (ξc-ξ)/(ξc-ξa)^2/√(1-((ξc-ξ)/(ξc-ξa))^2) : 0.0
    drξdξ = ξb-ξ!=0.0 ? -α/fc_peak * abs((ξb-ξ)/fc_peak)^(α-1) : 0.0
    dgdξ  = dgdrc*drcdξ + dgdrξ*drξdξ
    
    dξdσ = √3/3*I2
    dgdρ = 1.0

   
    # if ξ>ξb && (ξ-ξb)*dgdξ > ρ*dgdρ  # apex
    #     dgdσ = ξb*I2 - σ
    # else
    #     dgdσ = s/ρ + dgdξ*dξdσ
    # end
    dgdσ = s/ρ + dgdξ*dξdσ
    return dgdσ

end


function calcD(mat::CSCP, state::CSCPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stressmodel)

    state.Δλ==0.0 && return De

    dfdσ, dfdεcp, dfdεtp = yield_derivs(mat, state, state.σ, state.εcp, state.εtp)
    dgdσ = potential_derivs(mat, state, state.σ, state.εcp, state.εtp)

    Λ = eigvals(dgdσ, sort=false)
    # Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*norm(max.(0.0, Λ)))
    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*maximum(max.(0.0, Λ)))
    return Dep
end


function calc_σ_εp_Δλ(mat::CSCP, state::CSCPState, σtr::Vec6)

    maxits = 40
    maxits = 60
    tol    = 1.0
    dgdσ   = potential_derivs(mat, state, state.σ, state.εcp, state.εtp)
    De     = calcDe(mat.E, mat.ν, state.ctx.stressmodel)
    # Δλ     = norm(σtr-state.σ)/norm(De*dgdσ)/10
    Δλ     = eps()

    σ  = σtr - Δλ*(De*dgdσ)

    εcp = state.εcp
    εtp = state.εtp
    f   = yield_func(mat, state, state.σ, state.εcp, state.εtp)
    η   = 1.0 # initial damping
    
    # iterative process
    for i in 1:maxits
        dfdσ, _ = yield_derivs(mat, state, σ, εcp, εtp)
        dgdσ    = potential_derivs(mat, state, σ, εcp, εtp)

        dfdΔλ = -dfdσ'*De*dgdσ

        Δλ = Δλ - η*f/dfdΔλ
        if Δλ<0
            # Δλ = abs(Δλ)
            # @show Δλ
        end

        if isnan(Δλ) 
            return state.σ, 0.0, 0.0, 0.0, failure("CSCP: Δλ is NaN")
        end

        σ  = σtr - Δλ*(De*dgdσ)
        # εp = state.εp + Δλ*norm(dgdσ)

        Λ   = eigvals(dgdσ, sort=false)
        εcp = state.εcp + Δλ*norm(min.(0.0, Λ))
        # εtp = state.εtp + Δλ*norm(max.(0.0, Λ))
        εtp = state.εtp + Δλ*maximum(max.(0.0, Λ))
        f   = yield_func(mat, state, σ, εcp, εtp)

        if abs(f) < tol
            Δλ < 0.0 && return σ, 0.0, 0.0, 0.0, failure("CSCP: negative Δλ")

            return σ, εcp, εtp, Δλ, success()
        end

        # dumping
        i>10 && (η = 0.6)
        i>15 && (η = 0.3)
        # i>20 && (η = 0.15)
    end


    return state.σ, 0.0, 0.0, 0.0, failure("CSCP: maximum iterations reached")
end


function update_state!(mat::CSCP, state::CSCPState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stressmodel)
    Δσtr = De*Δε
    σtr  = state.σ + Δσtr
    ftr  = yield_func(mat, state, σtr, state.εcp, state.εtp)

    Δλ  = 0.0
    tol = 1.0

    if ftr < tol
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        state.σ, state.εcp, state.εtp, state.Δλ, status = calc_σ_εp_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status
    end

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::CSCP, state::CSCPState)
    σ, ε  = state.σ, state.ε
    # ρ  = √(2*J2(σ))
    # ξ  = tr(σ)/√3
    # θ  = calc_θ(mat, σ)
    # r  = calc_rθ(mat, σ)
    
    # fc = mat.fc_fun(state.εcp)
    # ft = mat.ft_fun(state.εtp)
    
    # ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, state.εcp, state.εtp)
    # rc = calc_rc(mat, ξa, ξ)
    # rξ = calc_rξ(mat, ξb, ξ)

    vals_d = stress_strain_dict(σ, ε, state.ctx.stressmodel)

    vals_d[:εcp] = state.εcp
    vals_d[:εtp] = state.εtp
    # vals_d[:ξ]   = ξ
    # vals_d[:ρ]   = ρ
    # vals_d[:θ]   = θ
    # vals_d[:fc]  = fc
    # vals_d[:ft]  = ft
    # vals_d[:ξa]  = ξa
    # vals_d[:ξb]  = ξb
    # vals_d[:κ]   = κ
    # vals_d[:r]   = r
    # vals_d[:rξ]  = rξ
    # vals_d[:rc]  = rc
    # vals_d[:ξc]  = 2*mat.fb/√3
    # vals_d[:fcb] = mat.fc_peak

    return vals_d
end
