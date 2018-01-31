# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJoint

mutable struct MCJointIpState<:IpState
    shared_data::SharedAnalysisData
    σ  ::Array{Float64,1}
    w  ::Array{Float64,1}
    wpa::Array{Float64,1}
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size
    function MCJointIpState(shared_data::SharedAnalysisData=SharedAnalysisData())
        this = new(shared_data)
        ndim = shared_data.ndim
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.wpa = zeros(ndim)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

mutable struct MCJoint<:Material
    E  ::Float64  # Young's modulus
    ν  ::Float64  # Poisson ratio
    σmax0::Float64  # tensile strength
    μ  ::Float64  # friction angle
    α  ::Float64  # elastic scale factor
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alpha=NaN, wc=NaN, ws=NaN)
        this = new(E, nu, ft, mu, alpha, wc, ws)
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

function calc_a_b(mat::MCJoint, upa)
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
        #b = 1.0
        #b  = -σs/(mat.wc-mat.ws)/1.0e9 # important in pull tests
        b  = -σs/(mat.wc-mat.ws)/1.0e6 # important in pull tests
    end
    return a, b
end

function calc_kn_ks(ipd::MCJointIpState, mat::MCJoint)
    αmin = 1.0
    α = max(mat.α - (mat.α - αmin)/mat.ws*ipd.upa, αmin) 
    kn = mat.E*α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*α/ipd.h
    return kn, ks
end


function yield_func(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    ndim = ipd.shared_data.ndim
    a, b = calc_a_b(mat, upa)
    σmax = a - b*upa
    if ndim==3
        return √(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1})
    ndim = ipd.shared_data.ndim
    if ndim==3
        τ = √(σ[2]^2 + σ[3]^2)
        if τ==0.0
            return [ 1., 0., 0.]
        else
            return [ mat.μ, σ[2]/τ, σ[3]/τ ]
        end
    else
        τ = σ[2]
        if τ==0.0
            return [ 1., 0]
        else
            return [ mat.μ, sign(τ) ]
        end
    end
end


function potential_derivs_test(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    ndim = ipd.shared_data.ndim
    a, b = calc_a_b(mat, upa)
    σmax = a - upa*b
    if ndim==3
        if σ[1]>=0.0
            # G1:
            τ = √(σ[2]^2 + σ[3]^2)
            if τ==0.0
                return [ 1., 0., 0.]
            else
                r = [ mat.μ, σ[2]/τ, σ[3]/τ ]
                return r/norm(r)
            end
            #σm = √(σ[2]^2 + σ[3]^2 + σmax^2*mat.μ^2)
            #if σm == 0.0 return [ 1., 0., 0. ] end
            #r  = [ σ[1]*mat.μ^2/σm, σ[2]/σm, σ[3]/σm ]
            #return r/norm(r)
        else
            # G2:
            τ = √(σ[2]^2 + σ[3]^2)
            if τ != 0.0
                r = [ 0.0, σ[2]/τ, σ[3]/τ ]
                r = r/norm(r)
            else
                r = [ 1., 0., 0.]
            end
            return r
        end
    else
        if σ[1]>=0.0
            # G1:
            σm = √(σ[2]^2 + σmax^2*mat.μ^2)
            if σm == 0.0 return [ 1., 0. ] end
            r = [ σ[1]*mat.μ^2/σm, σ[2]/σm ]
            return r/norm(r)
        else
            # G2:
            τ = σ[2]
            if τ==0.0
                return [ 1., 0]
            else
                return r = [ mat.μ, sign(τ) ]
            end
            return r
        end

    end
end

function potential_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    ndim = ipd.shared_data.ndim
    a, b = calc_a_b(mat, upa)
    σmax = a - upa*b
    if ndim==3
        if σ[1]>=0.0
            # G1:
            σm = √(σ[2]^2 + σ[3]^2 + σmax^2*mat.μ^2)
            if σm == 0.0 return [ 1., 0., 0. ] end
            r  = [ σ[1]*mat.μ^2/σm, σ[2]/σm, σ[3]/σm ]
            return r/norm(r)
        else
            # G2:
            τ = √(σ[2]^2 + σ[3]^2)
            if τ != 0.0
                r = [ 0.0, σ[2]/τ, σ[3]/τ ]
                r = r/norm(r)
            else
                r = [ 1., 0., 0.]
            end
            return r
        end
    else

        if a>0
            if σ[1] >= 0.0
                # G1:
                a, b = calc_a_b(mat, upa)
                σmax = a - b*upa
                r = [ 2*σ[1]/σmax^2, 2*σ[2]/(σmax^2*mat.μ^2) ]
                return r/norm(r)
            else
                τ = σ[2]
                r = [ 0., sign(τ) ]
            end
        else  # a==0
            if ipd.w[1] > 0.0
                return [ 1.0, 0.0 ]
            else
                τ = σ[2]
                r = [ 0., sign(τ) ]
            end
        end


        #=
        if σ[1] >= 0.0  && a != 0.0
            # G1:
            a, b = calc_a_b(mat, upa)
            σmax = a - b*upa
            r = [ 2*σ[1]/σmax^2, 2*σ[2]/(σmax^2*mat.μ^2) ]
            return r/norm(r)
        else
            # G2:
            τ = σ[2]
            #if τ != 0.0
            if w[1]>0.0
                r = [ 0., sign(τ) ]
            else
                r = [ 1., 0. ]
            end
            return r
        end
        =#


    end
end


function mountD(mat::MCJoint, ipd::MCJointIpState)
    ndim = ipd.shared_data.ndim
    kn, ks = calc_kn_ks(ipd, mat)
    #kn = mat.E*mat.α/ipd.h
    #G  = mat.E/(2.0*(1.0+mat.ν))
    #ks = G*mat.α/ipd.h
    if ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    if ipd.Δλ==0.0 # Elastic 
        #@show De
        return De
    else
        v    = yield_derivs(mat, ipd, ipd.σ)
        r    = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        y    = -mat.μ # dF/dσmax
        a, b = calc_a_b(mat, ipd.upa)
        m    = -b   # dσmax/dupa

        Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))
        #@show De
        #@show Dep

        return Dep
    end
end


function find_intersection(mat::MCJoint, ipd::MCJointIpState, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)

    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    maxits = 40
    for i=1:maxits
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, σ0 + α*Δσ, ipd.upa)

        if abs(F)<1e-7 break end

        if F<0.0
            α1 = α
            F1 = F
        else
            α2 = α
            F2 = F
        end
        #@show F
    end

    return α
end


function SecantRoot(f::Function, x0::Float64, x1::Float64, eps::Float64=1e-6)
    # initial values
    x  = x0
    f0 = f(x0)
    maxits = 20
    i=0
    for i=1:maxits
        f1 = f(x1)
        d  = (f1-f0)/(x1-x0) 
        x  = x1 - f1/d
        err = abs(x-x1)

        if err < eps
            return x, true
        end

        # updating
        x0 = x1
        x1 = x
        f0 = f1
    end

    return x, false
end


function stress_update(mat::MCJoint, ipd::MCJointIpState, Δw::Array{Float64,1})
    ndim = ipd.shared_data.ndim
    σini = copy(ipd.σ)
    μ    = mat.μ

    # calculate De
    kn, ks = calc_kn_ks(ipd, mat)
    #kn = mat.E*mat.α/ipd.h
    #G  = mat.E/(2.0*(1.0+mat.ν))
    #ks = G*mat.α/ipd.h

    if ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn  0.0
               0.0  ks ]
    end

    # Trial
    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)


    # Elastic and EP integration
    if Ftr<=0.0
        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr)

        # update w
        ipd.w += Δw

        # calculate Δσ
        Δσ = ipd.σ - σini
    else

        # Pure elastic increment
        α  = 0.0

        #
        # Find intersection with the yield surface
        Fini = yield_func(mat, ipd, ipd.σ, ipd.upa)
        if Ftr>1e-6 && Fini<0.0
            α = find_intersection(mat, ipd, Fini, Ftr, ipd.σ, De*Δw)
            # Update w and Δσ up to the intersection
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end
        #

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw

        upa = ipd.upa
        nincs = 1
        Δwi = Δwep/nincs
        for i=1:nincs
            σtr = ipd.σ + De*Δwi
            r   = potential_derivs(mat, ipd, ipd.σ, ipd.upa)

            a, b = calc_a_b(mat, ipd.upa)

            # Direct return to origin
            if b<1 && ipd.w[1]>0
                r = Δwi/norm(Δwi)
            end

            if ndim==3
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ^2-b*μ)
                #ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1]-b*μ)
            else
                #ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ^2-b*μ)
                #ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] - b*μ + ks*r[2]*sign(σtr[2]))
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] - b*μ + ks*abs(r[2]))
            end


            if ipd.Δλ<0.0
                warn("MCJoint: Negative plastic multiplier Δλ = $(ipd.Δλ). Try with a bigger value for α.")
                @show ipd.σ
                @show σtr
            end


            # Treating transition at ipd.ws
            if ipd.upa<mat.ws && ipd.upa+ipd.Δλ>=mat.ws
                Δλ1 = mat.ws - ipd.upa
                Δw1 = Δλ1/ipd.Δλ*Δwi
                σtr = ipd.σ + De*Δw1
                ipd.σ = σtr - Δλ1*De*r
                ipd.upa += Δλ1

                r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                a, b = calc_a_b(mat, ipd.upa)

                Δw2 = Δwi-Δw1
                σtr = ipd.σ + De*Δw2
                Δλ2 = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ*r[1] - b*μ + ks*abs(r[2]))
                ipd.σ = σtr - Δλ2*De*r
                ipd.upa += Δλ2
            else
                ipd.σ = σtr - ipd.Δλ*De*r
                ipd.upa += ipd.Δλ
            end


            a, b = calc_a_b(mat, ipd.upa)
            σmax = a - b*ipd.upa
            if b<0 && σtr[1]>0
                ipd.σ .= 1e-11
            end

            # Return to surface:
            F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

            #if false 
            if F>1e-6
                # The return algorithm needs improving
                Δσ = ipd.σ
                #σ0 = [-0.000001, 0.0, 0.0]
                σ0 = [min(-0.00000001,ipd.σ[1]), 0.0, 0.0][1:ndim]
                F0 = yield_func(mat, ipd, σ0, ipd.upa)
                @assert(F0*F<0.0)
                α  = find_intersection(mat, ipd, F0, F, σ0, Δσ)

                ipd.σ = σ0 + α*Δσ
            end

        end

        # Elastic plastic update of w and Δσ
        ipd.w += Δw
        Δσ = ipd.σ - σini
        F  = yield_func(mat, ipd, ipd.σ, ipd.upa)
    end

    return Δσ
end



function ip_state_vals(mat::MCJoint, ipd::MCJointIpState)
    ndim = ipd.shared_data.ndim
    if ndim == 2
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa   => ipd.upa
          )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] ,
          :upa   => ipd.upa
          )
    end
end

