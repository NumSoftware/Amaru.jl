# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearSlipTip

mutable struct LinearSlipTipState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearSlipTipState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct LinearSlipTip<:Material
    k::Float64

    function LinearSlipTip(prms::Dict{Symbol,Float64})
        return  LinearSlipTip(;prms...)
    end

    function LinearSlipTip(;k=NaN)
        @check k>=0
        this = new(k)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearSlipTip}) = LinearSlipTipState

# Element types that work with this material
compat_elem_types(::Type{LinearSlipTip}) = (MechSlipTip,)


function calcD(mat::LinearSlipTip, state::LinearSlipTipState)
    return mat.k
end


function update_state!(mat::LinearSlipTip, state::LinearSlipTipState, Δw)
    Δf = mat.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function ip_state_vals(mat::LinearSlipTip, state::LinearSlipTipState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
