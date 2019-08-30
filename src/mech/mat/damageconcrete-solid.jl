# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DamageConcrete

mutable struct DamageConcrete<:Material
    E0::Float64  # initial Young modulus
    ν::Float64
    ft::Float64
    GF::Float64
    fc::Float64
    εc0::Float64
    ρ::Float64

    function DamageConcrete(; E=NaN, nu=0.2, ft=NaN, GF=NaN, fc=NaN, epsc=NaN, rho=0.0)
        E >0.0 || error("Invalid value for rho: $E")
        ft>0.0 || error("Invalid value for ft: $ft")
        GF>0.0 || error("Invalid value for GF: $GF")
        fc<0.0 || error("Invalid value for fc: $fc")
        nu>=0.0 || error("Invalid value for nu: $nu")
        epsc<0.0 || error("Invalid value for epsc: $epsc")
        rho>=0.0 || error("Invalid value for rho: $rho")
        epsc/(fc/E)>1.0 || error("epsc should be greater in magnitude than fc/E")

        this = new(E, nu, ft, GF, fc, epsc, rho)
        return this
    end
end


mutable struct DamageConcreteIpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}  # current stress
    ε::Array{Float64,1}  # current strain
    ε̅cmax::Float64 
    ε̅tmax::Float64 
    h::Float64        # characteristic length
    in_linear_range::Bool
    damt::Float64
    damc::Float64

    function DamageConcreteIpState(env::ModelEnv=ModelEnv()) 
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.ε̅cmax = 0.0
        this.ε̅tmax = 0.0
        this.h  = 0.0 # will be set by elem_init
        this.in_linear_range = false
        this.damt = 0.0
        this.damc = 0.0

        return this
    end
end


# Returns the element type that works with this material model
matching_elem_type(::DamageConcrete) = MechSolid

# Type of corresponding state structure
ip_state_type(mat::DamageConcrete) = DamageConcreteIpState


function uniaxial_σ(mat::DamageConcrete, ipd::DamageConcreteIpState, εi::Float64)
    σp = eigvals(ipd.σ)
    σ1c, σ2c, σ3c = neg.(σp)
    if εi>=0  # tension: Nilsson and Oldenburg 1982; Beshara and Virdi 1991; Wu and Yao 1998
        αt = 2.05
        #αt = 1.5 #*
        γ = 1.0 - αt/mat.fc*( √(σ1c^2 + σ2c^2+ σ3c^2)) # suggested value for coef: 2.0
        #γ = 1.0 - αt/mat.fc*( abs(σ1c + σ2c+ σ3c)) # suggested value for coef: 2.0
        #γ = 1.0
        ft = mat.ft*γ
        εt0 = ft/mat.E0
        if εi<εt0
            return εi*mat.E0
        else
            #γ = 1.0 - 0.5/mat.fc*( √(σ1c^2 + σ2c^2+ σ3c^2))
            σ̅t = norm(pos.(σp))
            w = (εi-εt0)*ipd.h # crack openning
            #w = (εi-εt0*σ̅t/ft)*ipd.h # crack openning

            #GF = mat.GF - 0.5*ft*εt0*ipd.h
            #@show ipd.h
            #@show mat.GF
            #@show GF
            #@show 0.5*ft*εt0*ipd.h
            #error()
            #@assert GF>0
            #return ft*exp(-ft/GF*w)
            return ft*exp(-ft/mat.GF*w)
        end
    else # compression: Popovics 1973; Carreira and Chu 1985
        αc = 0.3
        #αc = 0.2
        #γ = 1.0 + αc*((σ2c*σ3c+ σ1c*σ3c+ σ1c*σ2c)/mat.fc^2)^0.25 # suggested values for coef: 0.2
        #γ = 1.0 + αc*((σ2c*σ3c+ σ1c*σ3c+ σ1c*σ2c)/mat.fc^2) # suggested values for coef: 0.2
        γ = 1.0 - αc/mat.fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c)) # suggested values for coef: 0.2
        fc = mat.fc*γ

        β = 1/(1-fc/(mat.εc0*mat.E0))
        β = max(min(β,10),2) # limit the value of β
        εr = εi/mat.εc0
        return fc*(β*εr)/(β - 1.0 + εr^β)
    end
end


function uniaxial_E(mat::DamageConcrete, ipd::DamageConcreteIpState, εi::Float64)
    σp = eigvals(ipd.σ)
    σ1c, σ2c, σ3c = neg.(σp)
    if εi>=0 # tension
        αt = 2.05
        #αt = 1.5 #*
        γ = 1.0 - αt/mat.fc*( √(σ1c^2 + σ2c^2+ σ3c^2)) # suggested value for coef: 2.0
        #γ = 1.0 - αt/mat.fc*( abs(σ1c + σ2c+ σ3c)) # suggested value for coef: 2.0
        #γ = 1.0 - 2.05/mat.fc*( abs(σ1c + σ2c+ σ3c)) # suggested value for coef: 2.0
        #γ = 1.0
        ft = mat.ft*γ
        εt0 = ft/mat.E0
        if εi<εt0
            return mat.E0
        else
            #GF = mat.GF - 0.5*ft*εt0*ipd.h
            #Et = ft*exp(-ft/(GF/ipd.h)*(εi-εt0)) * (-ft/(GF/ipd.h))
            σ̅t = norm(pos.(σp))
            #w = (εi-εt0*σ̅t/ft)*ipd.h # crack openning
            #Et = ft*exp(-ft/mat.GF*w) * (-ft/mat.GF*ipd.h)
            Et = ft*exp(-ft/(mat.GF/ipd.h)*(εi-εt0)) * (-ft/(mat.GF/ipd.h))
            return Et
        end
    else # compression
        αc = 0.3
        #αc = 0.2
        #αc = 0.25
        #αc = 0.05 #*
        γ = 1.0 - αc/mat.fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c)) # suggested values for coef: 0.2
        #γ = 1.0 - 0.15/mat.fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c)) # suggested values for coef: 0.4
        #γ = 1.0 + 0.2*((σ2c*σ3c+ σ1c*σ3c+ σ1c*σ2c)/mat.fc^2)^0.25 # suggested values for coef: 0.3
        #γ = 1.0 + 0.45(√(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c))/abs(mat.fc) # suggested values for coef: 0.45
        #γ = 1.0
        fc = mat.fc*γ

        εc0 = mat.εc0*γ
        εr  = εi/εc0

        β = 1/(1-fc/(εc0*mat.E0))
        Ec = β*(fc/εc0)/(β-1+εr^β) - β^2*(fc/εc0)*εr^β/(β-1+εr^β)^2
        #Ec *= γ
        return Ec
    end
end


function calcD(mat::DamageConcrete, ipd::DamageConcreteIpState)
    # special functions
    pos(x) = (abs(x)+x)/2.0
    neg(x) = (-abs(x)+x)/2.0
    σfun(εi) = uniaxial_σ(mat, ipd, εi)
    Efun(εi) = uniaxial_E(mat, ipd, εi)

    # principal strains
    εp = eigvals(ipd.ε)
    ε̅t = norm(pos.(εp))
    ε̅c = norm(neg.(εp))
    # principal stresses
    σp = eigvals(ipd.σ)
    σ̅c = norm(neg.(σp))
    σ̅t = norm(pos.(σp))

    # check for tension or compression dominant state
    if σ̅c==0 
        in_tension = true
    else
        in_tension = σ̅t/σ̅c > 0.01
        #in_tension = σ̅t/σ̅c > 0.01 || ε̅t/ε̅c > 0.01
    end

    #in_tension = true

    # estimate tangent Young modulus
    if in_tension
        if ipd.in_linear_range
            E = ε̅t==0 ? mat.E0 : min(mat.E0, σ̅t/ε̅t)
        else
            E = Efun(ε̅t)
        end
        ν = mat.ν*(1-ipd.damt)
        #ν = mat.ν
    else
        if ipd.in_linear_range
            E = ε̅c==0 ? mat.E0 : min(mat.E0, σ̅c/ε̅c)
        else
            E = Efun(-ε̅c)
        end
        #ν = mat.ν*(1-ipd.damc)
        ν = mat.ν
        #@show E
    end

    #@show in_tension
    #@show ipd.in_linear_range

    Emin = mat.E0*1e-6
    abs(E)<Emin && (E=Emin)


    #D  = calcDe(E, mat.ν, :general)
    #ν = mat.ν
    D  = calcDe(E, ν, :general)
    return D
end


#function calcDsec(mat::DamageConcrete, ipd::DamageConcreteIpState, Δε::Array{Float64,1}, modeltype::Symbol)
function stress_update(mat::DamageConcrete, ipd::DamageConcreteIpState, Δε::Array{Float64,1})
    # special functions
    pos(x) = (abs(x)+x)/2.0
    neg(x) = (-abs(x)+x)/2.0
    σfun(εi) = uniaxial_σ(mat, ipd, εi)
    Efun(εi) = uniaxial_E(mat, ipd, εi)


    εp = eigvals(ipd.ε)
    #ecc = norm(neg.(εp))
    #ett = norm(pos.(εp))


    ipd.ε .+= Δε

    # principal strains
    εp = eigvals(ipd.ε)
    ε̅t = norm(pos.(εp))
    ε̅c = norm(neg.(εp))
    # principal stresses
    σp = eigvals(ipd.σ)
    σ̅c = norm(neg.(σp))
    σ̅t = norm(pos.(σp))
    Δε̅tmax = max(ε̅t-ipd.ε̅tmax, 0.0)
    Δε̅cmax = max(ε̅c-ipd.ε̅cmax, 0.0)

    if σ̅c==0 
        in_tension = true
    else
        in_tension = σ̅t/σ̅c > 0.01
        in_tension = σ̅t/σ̅c > 0.01 || ε̅t > ε̅c
        #@show ε̅t/ε̅c
    end
    #in_tension =  true

    ipd.ε̅tmax = max(ε̅t, ipd.ε̅tmax)
    ipd.ε̅cmax = max(ε̅c, ipd.ε̅cmax)

    # estimate tangent Young modulus
    if in_tension
        if ε̅t < ipd.ε̅tmax
            E = ε̅t==0 ? mat.E0 : min(mat.E0, σ̅t/ε̅t)

            ipd.in_linear_range = true
        else
            E = Efun(ε̅t)
            ipd.in_linear_range = false
        end
        ν = mat.ν*(1-ipd.damt)
        #ν = mat.ν
    else
        if ε̅c < ipd.ε̅cmax
            E = ε̅c==0 ? mat.E0 : min(mat.E0, σ̅c/ε̅c)
            ipd.in_linear_range = true
        else
            E = Efun(-ε̅c)
            ipd.in_linear_range = false
        end
        #ν = mat.ν*(1-ipd.damc)
        ν = mat.ν
    end

    # Fix negative E values to avoid stress signal change in compression
    if E<0
        if σ̅c+E*Δε̅cmax < 0 && Δε̅cmax > Δε̅tmax
            E = -σ̅c/Δε̅cmax
        end
        if σ̅t+E*Δε̅tmax < 0 && Δε̅cmax < Δε̅tmax
            E = -σ̅t/Δε̅tmax
        end
    end

    if E>0
        #@show σ̅t
        #@show σ̅c
        #@show E
        #@show mat.E0
        #@show in_tension
        #@show ipd.in_linear_range
    end


    # update maximum strains
    #if in_tension
        #ipd.ε̅tmax = max(ε̅t, ipd.ε̅tmax)
    #else
        #ipd.ε̅cmax = max(ε̅c, ipd.ε̅cmax)
    #end

    Emin = mat.E0*1e-6
    abs(E)<Emin && (E=Emin)

    D  = calcDe(E, ν, :general)
    Δσ = D*Δε
    ipd.σ .+= Δσ

    ipd.damt = clamp(1.0 - σfun(ipd.ε̅tmax)/ipd.ε̅tmax/mat.E0, 0.0, 1.0)
    ipd.damc = clamp(1.0 + σfun(-ipd.ε̅cmax)/ipd.ε̅cmax/mat.E0, 0.0, 1.0)
    #ipd.damt = 1.0 - σfun(ipd.ε̅tmax)/ipd.ε̅tmax/mat.E0
    #ipd.damc = 1.0 + σfun(-ipd.ε̅cmax)/ipd.ε̅cmax/mat.E0

    #if ipd.damt<0 || ipd.damt>1
        #@show ipd.damt
        #@show σfun(ipd.ε̅tmax)
        #@show ipd.ε̅tmax
        #error("damage error")
    #end

    #ipd.env.cstage==1 && ipd.env.cinc==2 && error()

    return Δσ
end

function ip_state_vals(mat::DamageConcrete, ipd::DamageConcreteIpState)
    dict = stress_strain_dict(ipd.σ, ipd.ε, ipd.env.ndim)
    dict[:damt] = ipd.damt
    dict[:damc] = ipd.damc

    εp = eigvals(ipd.ε)
    ε̅t = norm(pos.(εp))
    Efun(εi) = uniaxial_E(mat, ipd, εi)
    εt0 = mat.ft/mat.E0
    dict[:w] = abs(ε̅t-εt0)*ipd.h # crack openning
    return dict
end
