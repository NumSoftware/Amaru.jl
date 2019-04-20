# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export NIConcrete


mutable struct NIConcrete<:Material
    E0::Float64  # initial Young modulus
    ν::Float64
    ft::Float64
    GF::Float64
    fc::Float64
    εc0::Float64
    ρ::Float64

    function NIConcrete(prms::Dict{Symbol,Float64})
        return  NIConcrete(;prms...)
    end

    function NIConcrete(; E=NaN, nu=0.2, ft=NaN, GF=NaN, fc=NaN, epsc=NaN, rho=0.0)
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


mutable struct NIConcreteIpState<:IpState
    analysis_data::AnalysisData
    σ::Array{Float64,1}  # current stress
    ε::Array{Float64,1}  # current strain
    ε̅cmax::Float64 
    ε̅tmax::Float64 
    h::Float64        # characteristic length
    in_linear_range::Bool
    damt::Float64
    damc::Float64

    NIConcreteIpState() = new()

    function NIConcreteIpState(mat::NIConcrete, analysis_data::AnalysisData=AnalysisData()) 
        this = new(analysis_data)
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
matching_elem_type(::NIConcrete) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::NIConcrete, analysis_data::AnalysisData) = NIConcreteIpState(mat, analysis_data)


function uniaxial_σ(mat::NIConcrete, ipd::NIConcreteIpState, εi::Float64)
    if εi>=0  # tension: Nilsson and Oldenburg 1982; Beshara and Virdi 1991; Wu and Yao 1998
        εt0 = mat.ft/mat.E0
        if εi<εt0
            return εi*mat.E0
        else
            return mat.ft*exp(-mat.ft/(mat.GF/ipd.h)*(εi-εt0))
        end
    else # compression: Popovics 1973; Carreira and Chu 1985
        #fc = mat.fc
        #fcc = 1.5*fc
        #σ1c, σ2c, σ3c = neg.(eigvals(ipd.σ))
        #γ = 1.0 - (fcc-fc)/fc^2*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c))


        β = 1/(1-mat.fc/(mat.εc0*mat.E0))
        β = max(min(β,10),2) # limit the value of β

        εr = εi/mat.εc0
        return mat.fc*(β*εr)/(β - 1.0 + εr^β)
        #return γ*mat.fc*(β*εr)/(β - 1.0 + εr^β)
    end
end


function uniaxial_E(mat::NIConcrete, ipd::NIConcreteIpState, εi::Float64)
    if εi>=0 # tension
        εt0 = mat.ft/mat.E0
        if εi<εt0
            return mat.E0
        else
            Et = mat.ft*exp(-mat.ft/(mat.GF/ipd.h)*(εi-εt0)) * (-mat.ft/(mat.GF/ipd.h))
            return Et
        end
    else # compression

        fc = mat.fc
        #fcc = 1.41*fc
        σ1c, σ2c, σ3c = neg.(eigvals(ipd.σ))
        #γ = 1.0 - (fcc-fc)/fc^2*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c))
        α = 0.37
        #α = 0.00
        γ = 1.0 - α/fc*( √(σ2c*σ3c)+ √(σ1c*σ3c)+ √(σ1c*σ2c))
        #@show γ


        εc0 = mat.εc0
        εr  = εi/εc0

        #εr  = εi/εc0/γ

        β = 1/(1-mat.fc/(εc0*mat.E0))
        #β = max(min(β,10),2)
        Ec = β*(mat.fc/εc0)/(β-1+εr^β) - β^2*(mat.fc/εc0)*εr^β/(β-1+εr^β)^2
        return Ec
    end
end


function calcD(mat::NIConcrete, ipd::NIConcreteIpState)
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
        in_tension = σ̅t/σ̅c > abs(mat.ft/mat.fc)
        #in_tension = σ̅t/σ̅c > abs(mat.ft/mat.fc) || ε̅t>ε̅c
        #in_tension = ε̅t>ε̅c
        #in_tension = ε̅t/ε̅c>0.1
        #in_tension = σ̅t/σ̅c > 0.01
    end

    # estimate tangent Young modulus
    if in_tension
        if ipd.in_linear_range
            #@show "linear tension"
            #E = σ̅t/ε̅t
            #E = ε̅t==0 ? mat.E0 : σ̅t/ε̅t
            E = ε̅t==0 ? mat.E0 : min(mat.E0, σ̅t/ε̅t)
            #if E>2e7
                #@show "linear tension"
                #@show E
                #error()
            #end
        else
            #@show "tension"
            E = Efun(ε̅t)

            #if σ̅c>σ̅t
                #E = sign(E)*√abs(E*Efun(ε̅c))
                #E = Efun(ε̅c)
            #end
        end
    else
        if ipd.in_linear_range
            #@show "linear compression"
            #E = ε̅c==0 ? mat.E0 : σ̅c/ε̅c
            E = ε̅c==0 ? mat.E0 : min(mat.E0, σ̅c/ε̅c)
            #if E>2e7
                #@show "linear compression"
                #@show E
                #error()
            #end
        else
            #@show "compression"
            E = Efun(-ε̅c)
        end
    end

    Emin = mat.E0*1e-6
    abs(E)<Emin && (E=Emin)

    return calcDe(E, mat.ν, :general)
end


#function calcDsec(mat::NIConcrete, ipd::NIConcreteIpState, Δε::Array{Float64,1}, model_type::Symbol)
function stress_update(mat::NIConcrete, ipd::NIConcreteIpState, Δε::Array{Float64,1})
    # special functions
    pos(x) = (abs(x)+x)/2.0
    neg(x) = (-abs(x)+x)/2.0
    σfun(εi) = uniaxial_σ(mat, ipd, εi)
    Efun(εi) = uniaxial_E(mat, ipd, εi)


    εp = eigvals(ipd.ε)
    ecc = norm(neg.(εp))
    ett = norm(pos.(εp))


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
        in_tension = σ̅t/σ̅c > abs(mat.ft/mat.fc)
        #in_tension = σ̅t/σ̅c > abs(mat.ft/mat.fc) || ε̅t>ε̅c
        #in_tension = ε̅t>ε̅c
        #in_tension = ε̅t/ε̅c>0.1
        #in_tension = σ̅t/σ̅c > 0.01
        #@show σ̅t/σ̅c
        #@show ε̅t/ε̅c
    end
    #@show in_tension
    #@show (mat.ft/mat.E0)
    #@show abs(mat.εc0)
    #@show ε̅t/ε̅c
    #@show (mat.ft/mat.E0) / abs(mat.εc0)
    #@show σp
    #@show σ̅t/σ̅c
    #@show ipd.ε̅cmax

    ipd.ε̅tmax = max(ε̅t, ipd.ε̅tmax)
    ipd.ε̅cmax = max(ε̅c, ipd.ε̅cmax)

    # estimate tangent Young modulus
    if in_tension
        if ε̅t < ipd.ε̅tmax
            #@show "int linear tension"
            #E = σ̅t/ε̅t
            #E = ε̅t==0 ? mat.E0 : σ̅t/ε̅t
            E = ε̅t==0 ? mat.E0 : min(mat.E0, σ̅t/ε̅t)

            ipd.in_linear_range = true
        else
            #@show "int tension"
            E = Efun(ε̅t)
            #E = abs(σfun(ε̅t)-σfun(ett))/Δε̅tmax
            ipd.in_linear_range = false
        end
    else
        if ε̅c < ipd.ε̅cmax
            #@show "int linear compression"
            #@show σp
            #@show σ̅c
            #@show ε̅c
            #@show ipd.ε̅cmax
            #E = σ̅c/ε̅c
            #E = ε̅c==0 ? mat.E0 : σ̅c/ε̅c
            E = ε̅c==0 ? mat.E0 : min(mat.E0, σ̅c/ε̅c)
            ipd.in_linear_range = true
        else
            #@show "int compression"
            E = Efun(-ε̅c)
            #E = Efun(-ecc)
            #E = abs(σfun(-ε̅c)-σfun(-ecc))/Δε̅cmax
            ipd.in_linear_range = false
        end
    end

    # Fix negative E values to avoid stress signal change in compression
    if E<0
        #@show Δε̅cmax
        #@show Δε̅tmax
        #if σ̅c+E*Δε̅cmax < 0
        if σ̅c+E*Δε̅cmax < 0 && Δε̅cmax > Δε̅tmax
            E = -σ̅c/Δε̅cmax
        end
    end

    # update maximum strains
    #if in_tension
        #ipd.ε̅tmax = max(ε̅t, ipd.ε̅tmax)
    #else
        #ipd.ε̅cmax = max(ε̅c, ipd.ε̅cmax)
    #end

    Emin = mat.E0*1e-6
    abs(E)<Emin && (E=Emin)


    D  = calcDe(E, mat.ν, :general)
    Δσ = D*Δε
    ipd.σ .+= Δσ

    #@show σfun(ipd.ε̅tmax)/ipd.ε̅tmax
    #@show -σfun(-ipd.ε̅cmax)/ipd.ε̅cmax
    #@show σfun(ipd.ε̅tmax)/ipd.ε̅tmax/mat.E0
    #@show -σfun(-ipd.ε̅cmax)/ipd.ε̅cmax/mat.E0

    ipd.damt = 1.0 - σfun(ipd.ε̅tmax)/ipd.ε̅tmax/mat.E0
    ipd.damc = 1.0 + σfun(-ipd.ε̅cmax)/ipd.ε̅cmax/mat.E0

    #ipd.analysis_data.nstage==1 && ipd.analysis_data.ninc==2 && error()

    return Δσ
end

function ip_state_vals(mat::NIConcrete, ipd::NIConcreteIpState)
    dict = stress_strain_dict(ipd.σ, ipd.ε, ipd.analysis_data.ndim)
    dict[:damt] = ipd.damt
    dict[:damc] = ipd.damc
    return dict
end
