# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CSCP

CSCP_params = [
    FunInfo(:CSCP, "Closed surface concrete plasticity model"),
    KwArgInfo(:E, "Young's modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson's ratio", cond=:(nu>=0)),
    KwArgInfo(:ft, "Tensile strength", cond=:(ft>0)),
    KwArgInfo(:fc, "Compressive strength", cond=:(fc<0)),
    KwArgInfo(:fb, "Biaxial compressive strength", cond=:(fb<0)),
    KwArgInfo(:ft_fun, "Plastic function for ft", nothing),
    KwArgInfo(:fc_fun, "Plastic function for fc", nothing),
    KwArgInfo(:fic_fun, "Plastic function for fic", nothing),
]

mutable struct CSCP<:Material
    E::Float64
    ν::Float64
    e::Float64  # eccentricity of the yield surface
    ft_fun::Union{Nothing,PathFunction}
    fc_fun::Union{Nothing,PathFunction}
    fic_fun::Union{Nothing,PathFunction}

    function CSCP(; kwargs...)
        args = checkargs(kwargs, CSCP_params)
        
        fc_fun, ft_fun, fic_fun = args.fc_fun, args.ft_fun, args.fic_fun

        fc, fb, ft = abs(args.fc), abs(args.fb), args.ft
        e  = (fc*fb - fc*ft + 3*fb*ft)/(2*fc*fb + ft*fc)
        @assert 0<e<=1

        this = new(args.E, args.nu, e, ft_fun, fc_fun, fic_fun)
        return this
    end
end


mutable struct CSCPState<:IpState
    env::ModelEnv
    σ::Vec6
    ε::Vec6
    εp::Float64
    Δλ::Float64
    function CSCPState(env::ModelEnv)
        this = new(env)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εp = 0.0
        this.Δλ  = 0.0
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{CSCP}, ::Type{MechSolid}, env::ModelEnv) = env.ana.stressmodel!=:planestress ? CSCPState : error("CSCP: This model is not compatible with planestress")


function calc_r(mat::CSCP, σ::Vec6)
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


function yield_func(mat::CSCP, σ::AbstractArray, εp::Float64)
    e = mat.e
    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)
    r = calc_r(mat, σ)

    fc = mat.fc_fun(εp)
    ft = mat.ft_fun(εp)
    ξa = √3*mat.fic_fun(εp)
    Ω  = √((ft-√3*ξa)/(fc-√3*ξa))
    ξb = fc*ft/√3 * (1+e*Ω)/(ft+fc*e*Ω)
    c  = (fc/√3-ξa)*(fc/√3-ξb)^2/(2/3*fc^2)

    return ρ^2 - r/c*(ξ-ξa)*(ξ-ξb)^2
end


function yield_derivs(mat::CSCP, σ::AbstractArray, εp::Float64)
    e = mat.e
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

    # @show εp

    fc = mat.fc_fun(εp)
    ft = mat.ft_fun(εp)
    ξa = √3*mat.fic_fun(εp)
    Ω  = √((ft-√3*ξa)/(fc-√3*ξa))
    ξb = fc*ft/√3 * (1+e*Ω)/(ft+fc*e*Ω)
    c  = (fc/√3-ξa)*(fc/√3-ξb)^2/(2/3*fc^2)
    
    # f derivative w.r.t. σ:

    dfdρ = 2*ρ
    dfdξ = -r/c*((ξ-ξb)^2 + 2*(ξ-ξa)*(ξ-ξb))
    dfdr = -1/c*(ξ-ξa)*(ξ-ξb)^2
    dfdθ = dfdr*drdθ

    dρdσ = s/norm(s)
    dξdσ = √3/3*I2
    dsdσ = Psd
    dθdσ = dsdσ*dθds
    
    dfdσ = dfdρ*dρdσ + dfdξ*dξdσ + dfdθ*dθdσ

    # f derivative w.r.t. ξa and ξb:
    # f = ρ^2 - r*h/c
    
    # h = -r*(ξ-ξa)*(ξ-ξb)^2
    
    # dhdξa = r*(ξ-ξb)^2
    # dcdξa = -(fc/√3-ξb)^2/(2/3*fc^2)
    # dfdξa = (dhdξa*c - h*dcdξa)/c^2
    
    # dhdξb = 2*r*(ξ-ξa)*(ξ-ξb)
    # dcdξb = -2*(fc/√3-ξa)*(fc/√3-ξb)/(2/3*fc^2)
    # dfdξb = (dhdξb*c - h*dcdξb)/c^2

    # f_ξa = ξa -> begin
    #     c = (fc/√3-ξa)*(fc/√3-ξb)^2/(2/3*fc^2)
    #     return ρ^2 - r/c*(ξ-ξa)*(ξ-ξb)^2
    # end
    # f_ξb = ξb -> begin 
    #     c = (fc/√3-ξa)*(fc/√3-ξb)^2/(2/3*fc^2)
    #     ρ^2 - r/c*(ξ-ξa)*(ξ-ξb)^2
    # end
    # dfdξa = derive(f_ξa, ξa)
    # dfdξb = derive(f_ξb, ξb)

    # dξadεp = √3*derive(mat.fic_fun, εp)

    # ξb_εp = εp -> begin
    #     fc = mat.fc_fun(εp)
    #     ft = mat.ft_fun(εp)
        
    #     Ω  = √((ft-√3*ξa)/(fc-√3*ξa))
        
    #     ξb = fc*ft/√3 * (1-e*Ω)/(ft-fc*e*Ω)
    #     return ξb
    # end
    # dξbdεp = derive(ξb_εp, εp)


    f_εp = εp -> yield_func(mat, σ, εp)

    dfdεp = derive(f_εp, εp)

    # dfdεp = dfdξa*dξadεp + dfdξb*dξbdεp

    # @show dfdξa
    # @show dξadεp
    # @show dfdξb
    # @show dξbdεp
    # @show εp

    # @show dfdεp

    return dfdσ, dfdεp
end


function potential_derivs(mat::CSCP, σ::AbstractArray, εp::Float64)
    # g(σ) = ρ^2 - 1/c*(ξ-ξa)*(ξ-ξb)^2

    e = mat.e
    i1= tr(σ)
    ξ = i1/√3

    s  = dev(σ)
    fc = mat.fc_fun(εp)
    ft = mat.ft_fun(εp)

    ξa = √3*mat.fic_fun(εp)
    Ω  = √((ft-√3*ξa)/(fc-√3*ξa))
    ξb = fc*ft/√3 * (1+e*Ω)/(ft+fc*e*Ω)

    c  = (fc/√3-ξa)*(fc/√3-ξb)^2/(2/3*fc^2)

    dgdσ = s - 1/c*((ξ-ξb)^2 + 2(ξ-ξa)*(ξ-ξb))*√3/3*I2
    return dgdσ

end


function calcD(mat::CSCP, state::CSCPState)
    De  = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)

    state.Δλ==0.0 && return De

    j2 = J2(state.σ)

    if j2 != 0.0
        dfdσ, dfdεp = yield_derivs(mat, state.σ, state.εp)
        dgdσ = potential_derivs(mat, state.σ, state.εp)
    else # apex
        dfdσ = √3/3*I2
        dgdσ = dfdσ
    end

    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεp*norm(dgdσ))
    return Dep

end


function calc_σ_εp_Δλ(mat::CSCP, state::CSCPState, σtr::Vec6)

    maxits = 10
    tol    = 1e-3
    dgdσ   = potential_derivs(mat, state.σ, state.εp)
    De     = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)
    Δλ     = norm(σtr-state.σ)/norm(De*dgdσ)

    @show σtr


    # iterative process
    for i in 1:maxits
        println()
        @show i
        σ    = σtr
        dgdσ = I2
        εp   = state.εp

        f_Δλ = Δλ -> begin
            @show Δλ
            for j in 1:8
                @show j
                dgdσ = potential_derivs(mat, σ, εp)
                σ = σtr - Δλ*(De*dgdσ)
                @show dgdσ
                @show Δλ*(De*dgdσ)
                @show σ
            end
            println()
            @show εp
            @show Δλ
            @show σ
            εp = state.εp + Δλ*norm(dgdσ)
            yield_func(mat, σ, εp)
        end


        Δλ = Δλ - f_Δλ(Δλ)/derive(f_Δλ, Δλ)
        @show Δλ

        error()


        if f_Δλ(Δλ) < tol
            @assert Δλ >= 0.0
            return σ, εp, Δλ, success()
        end
    end

    return state.σ, 0.0, 0.0, failure("CSCP: maximum iterations reached")

end


function update_state!(mat::CSCP, state::CSCPState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, σtr, state.εp)
    Δλ = 0.0

    if ftr < 1e-8
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        # σ, εp, Δλ, status = calc_σ_εp_Δλ(mat, state, σtr)
        # failed(status) && return state.σ, status
        
        # state.σ, state.εp, state.Δλ = σ, εp, Δλ

        nsteps = 10
        ΔΔε = Δε/nsteps
        for i in 1:nsteps
            j2 = J2(state.σ)
            
            # Derivatives
            if j2 != 0.0
                dfdσ, dfdεp = yield_derivs(mat, state.σ, state.εp)
                dgdσ = potential_derivs(mat, state.σ, state.εp)
            else # apex
                dfdσ = √3/3*I2
                dgdσ = dfdσ
            end

            ΔΔσtr = De*ΔΔε
            ftr = yield_func(mat, state.σ + ΔΔσtr, state.εp)

            if ftr < 1e-8
                ΔΔσ = ΔΔσtr
                ΔΔεp = 0.0
            else
                Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεp*norm(dgdσ))
                ΔΔσ = Dep*ΔΔε
                if dfdεp!=0.0
                    # ΔΔεp = -dot(dfdσ, ΔΔσ)/dfdεp
                    Δλ = dot(dfdσ, ΔΔσ)/(-dfdεp*norm(dgdσ))
                    # @show Δλ
                    @assert Δλ >= 0.0
                else
                    Δλ = norm(ΔΔσtr-ΔΔσ)/norm(De*dgdσ)
                    # @show Δλ
                end
                ΔΔεp = Δλ*norm(dgdσ)

            end

            state.σ += ΔΔσ
            state.εp += ΔΔεp

        end

        state.Δλ = Δλ
    end

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::CSCP, state::CSCPState)
    σ, ε  = state.σ, state.ε
    ρ = √(2*J2(σ))
    ξ = tr(σ)/√3

    D       = stress_strain_dict(σ, ε, state.env.ana.stressmodel)
    D[:ep]  = state.εp
    D[:xi]  = ξ
    D[:rho] = ρ

    return D
end
