# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJoint

mutable struct MCJointIpState<:IpState
    shared_data::SharedAnalysisData
    σ  ::Array{Float64,1}
    w  ::Array{Float64,1}
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size
    function MCJointIpState(shared_data::SharedAnalysisData=SharedAnalysisData())
        this = new(shared_data)
        ndim = shared_data.ndim
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

mutable struct MCJoint<:Material
    E  ::Float64  # Young's modulus from bulk material
    ν  ::Float64  # Poisson ratio from bulk material
    σmax0::Float64  # tensile strength
    μ  ::Float64  # friction angle
    α  ::Float64  # elastic displacement scale factor
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    model::String # softening model ("bilinear" or "hordijk")

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alpha=NaN, wc=NaN, ws=NaN, model="bilinear")
        this = new(E, nu, ft, mu, alpha, wc, ws, model)
        return this
    end
end

# Returns the element type that works with this material model
@static if isdefined(:MechJoint)
    matching_elem_type(::MCJoint) = MechJoint
end

# Create a new instance of Ip data
new_ip_state(mat::MCJoint, shared_data::SharedAnalysisData) = MCJointIpState(shared_data)

function set_state(ipd::MCJointIpState, sig=zeros(0), eps=zeros(0))
    @assert(false)
end

function calc_σmax(mat::MCJoint, upa)
    if mat.model == "bilinear"
        σs = 0.25*mat.σmax0
        if upa<mat.ws
            a  = mat.σmax0 
            b  = (mat.σmax0 - σs)/mat.ws
        elseif upa<mat.wc
            a  = mat.wc*σs/(mat.wc-mat.ws)
            b  = σs/(mat.wc-mat.ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.model == "hordijk"
        if upa < mat.wc
            z = (1 + 27*(upa/mat.wc)^3)*e^(-6.93*upa/mat.wc) - 28*(upa/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.σmax0
    end
    return σmax
end

function σmax_derivs(mat::MCJoint, upa)
    # dσmax = ∂σmax/∂upa
    if mat.model == "bilinear"
        σs = 0.25*mat.σmax0
        if upa<mat.ws
            b  = (mat.σmax0 - σs)/mat.ws
        elseif upa<mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.model == "hordijk"
        if upa < mat.wc
            dz = ((81*upa^2*e^(-6.93*upa/mat.wc)/mat.wc^3) - (6.93*(1 + 27*upa^3/mat.wc^3)*e^(-6.93*upa/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.σmax0
    end
    return dσmax
end

function calc_kn_ks(ipd::MCJointIpState, mat::MCJoint)
    if mat.model == "bilinear"
        if ipd.upa<mat.ws
            α = mat.α - (mat.α - 2.0)*ipd.upa/mat.ws
        elseif ipd.upa<mat.wc
            α = 2.0 - 1.5*ipd.upa/mat.wc
        else
            α = 0.5
        end
    elseif mat.model == "hordijk"
        z = (1 + 27*(ipd.upa/mat.wc)^3)*e^(-6.93*ipd.upa/mat.wc) - 28*(ipd.upa/mat.wc)*e^(-6.93)
        if ipd.upa<mat.wc
            α = mat.α - (mat.α - 0.5)*(1-z)
        else
            α = 0.5
        end 
    end
    kn = mat.E*α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*α/ipd.h
    return kn, ks
end


function yield_func(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    ndim = ipd.shared_data.ndim
    σmax = calc_σmax(mat, upa)
    if ndim==3
        return abs(σ[2]) + abs(σ[3]) + (σ[1]-σmax)*mat.μ
        #return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1})
    ndim = ipd.shared_data.ndim
    if ndim==3
        return [ mat.μ, sign(σ[2]), sign(σ[3])]
        #return [ mat.μ, σ[2]/sqrt(σ[2]^2 + σ[3]^2), σ[3]/sqrt(σ[2]^2 + σ[3]^2)]
    else
        return [ mat.μ, sign(σ[2]) ]
    end
end


function potential_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    σmax = calc_σmax(mat, upa)
    ndim = ipd.shared_data.ndim
    if ndim==3
        if σmax>0.0
            if σ[1] >= 0.0 
                # G1:
                r = [ 2.0*σ[1]/σmax^2, 2.0*σ[2]/(σmax^2*mat.μ^2), 2.0*σ[3]/(σmax^2*mat.μ^2)]
            else
                # G2:
                r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
            end
        else # σmax==0.0    
                r = [ 0.0, 2*σ[2], 2*σ[3] ]
        end
    else
        if σmax>0.0
            if σ[1] >= 0.0 
                # G1:
                r = [ 2*σ[1]/σmax^2, 2*σ[2]/(σmax^2*mat.μ^2) ]
            else
                # G2:
                r = [ 0.0, sign(σ[2]) ]
            end
        else # σmax==0.0    
                r = [ 0.0, sign(σ[2]) ]
        end
    end
    return r/norm(r)
end

function find_intersection(mat::MCJoint, ipd::MCJointIpState, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)

    # Regula Falsi method
    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    maxits = 40
    for i=1:maxits
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, σ0 + α*Δσ, ipd.upa)
        abs(F)<1e-7 && break

        if F<0.0
            α1 = α
            F1 = F
        else
            α2 = α
            F2 = F
        end
    end

    return α
end

function mountD(mat::MCJoint, ipd::MCJointIpState)
    kn, ks = calc_kn_ks(ipd, mat)
    σmax = calc_σmax(mat, ipd.upa)
    ndim = ipd.shared_data.ndim

    if ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    if ipd.Δλ == 0.0 # Elastic 
        return De
    elseif σmax == 0.0 # No more tensile strength 
        return De*1e-10
    else # Elastic-plastic
        v = yield_derivs(mat, ipd, ipd.σ)
        r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        y = -mat.μ # ∂F/∂σmax
        m = σmax_derivs(mat, ipd.upa)  # ∂σmax/∂upa

        Dep  = De - De*r*v'*De/(v'*De*r - y*m)
        return Dep
    end
end

function stress_update(mat::MCJoint, ipd::MCJointIpState, Δw::Array{Float64,1})
   
    σini = copy(ipd.σ)

    μ = mat.μ
    kn, ks = calc_kn_ks(ipd, mat)
    σmax = calc_σmax(mat, ipd.upa)  
    ndim = ipd.shared_data.ndim

    if ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn  0.0
               0.0  ks ]
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa) 

    # Elastic and EP integration

    if σmax == 0.0 && ipd.w[1] >= 0.0 
        if ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = (abs(σtr[2]) + abs(σtr[3]) + σtr[1]*μ)/(kn*μ*r[1] + ks*abs(r[2]) + ks*abs(r[3]))
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks]
            r = r1/norm(r1)
            ipd.Δλ = (abs(σtr[2]) + σtr[1]*μ)/(kn*μ*r[1] + ks*abs(r[2]))     
        end

        ipd.upa += ipd.Δλ
        ipd.σ = σtr - ipd.Δλ*De*r

        # Plastic update of w and Δσ
        ipd.w += Δw
        Δσ = ipd.σ - σini

    elseif Ftr <= 0.0

        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr) 

        # Plastic update of w and Δσ
        ipd.w += Δw
        Δσ = ipd.σ - σini

    else
        # Pure elastic increment
        α  = 0.0

        # Find intersection with the yield surface
        Fini = yield_func(mat, ipd, ipd.σ, ipd.upa)

        if Ftr>1e-6 && Fini<0.0
            α = find_intersection(mat, ipd, Fini, Ftr, ipd.σ, De*Δw)

            # Elastic update of w and Δσ
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end

        σi = copy(ipd.σ)
        upai = copy(ipd.upa)

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw

        σtr = ipd.σ + De*Δwep
    
        v = yield_derivs(mat, ipd, ipd.σ)
        r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        m = calc_σmax(mat, ipd.upa)
    
            if ndim==3
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2]) + ks*abs(r[3]))
            else
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2]))
            end

            if ipd.Δλ<0.0
                warn("MCJoint: Negative plastic multiplier Δλ = $(ipd.Δλ).")
                @show ipd.σ
                @show σtr
            end
            
            if mat.model == "bilinear"
                upa = ipd.upa + ipd.Δλ
                
                if ipd.upa<mat.ws
                    a  = mat.σmax0 
                elseif ipd.upa<mat.wc
                    a  = mat.wc*0.25*mat.σmax0/(mat.wc-mat.ws)
                else
                    a = 0.0
                end

                if a  == mat.σmax0 && upa < mat.ws
                    ipd.upa += ipd.Δλ
                    ipd.σ = σtr - ipd.Δλ*De*r
                elseif a  == mat.σmax0 && upa > mat.ws   
                    σs = 0.25*mat.σmax0
                    Δλ1 = (mat.ws - ipd.upa)
                    Δw1 = (Δλ1/ipd.Δλ)*Δwep
                    Δw2 = Δwep-Δw1

                    ipd.upa += Δλ1
                    ipd.σ = σtr - Δλ1*De*r
                    σtr = ipd.σ + De*Δw2

                    v = yield_derivs(mat, ipd, ipd.σ)
                    r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                    a = mat.wc*σs/(mat.wc-mat.ws)
                    b = σs/(mat.wc-mat.ws)
                    m = -b 
                    
                    if ndim==3
                        ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2])+ ks*abs(r[3]))
                    else
                        ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2]))
                    end

                    ipd.upa += ipd.Δλ
                    ipd.σ = σtr - ipd.Δλ*De*r
                else
                    ipd.upa += ipd.Δλ
                    ipd.σ = σtr - ipd.Δλ*De*r
                end
            elseif mat.model == "hordijk"
                    ipd.upa += ipd.Δλ
                    ipd.σ = σtr - ipd.Δλ*De*r
            end

            # Return to surface:
            F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

            #if false 
            if F>1e-6 || σini[1] > ipd.σ[1] || σini[2] > ipd.σ[2] 
            	ipd.upa = upai
                ipd.σ = σi
                #Δwc = mat.wc/200

                #if ndim==3
                #    n = round(Int, max(Δwep[1]/Δwc + 1, Δwep[2]/Δwc + 1, Δwep[3]/Δwc + 1, 50))
                #    n = min(n, 100)
                #else
                #    n = round(Int, max(Δwep[1]/Δwc + 1, Δwep[2]/Δwc + 1, 50))
                #    n = min(n, 100)
                #end
            
                #Δwinc = Δwep/n
                n = 100
                Δwinc = Δwep/n

                for i= 1:n

            		σtr = ipd.σ + De*Δwinc

            		v   = yield_derivs(mat, ipd, ipd.σ)
        			r   = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        			m   = calc_σmax(mat, ipd.upa)

                	if ndim==3
                        ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2])+ ks*abs(r[3]))
                    else
                        ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] + m*μ + ks*abs(r[2]))
                    end

            		ipd.upa += ipd.Δλ
                	ipd.σ = σtr - ipd.Δλ*De*r
                end
            end
            
        #  Plastic update of w and Δσ
        ipd.w += Δwep
        Δσ = ipd.σ - σini
    end

    return Δσ
end

function ip_state_vals(mat::MCJoint, ipd::MCJointIpState)
    ndim = ipd.shared_data.ndim
    if ndim == 3
       return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] ,
          :upa   => ipd.upa
          )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa   => ipd.upa
          )
    end
end
