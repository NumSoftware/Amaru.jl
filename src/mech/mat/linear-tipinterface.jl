# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearTipContact

mutable struct LinearTipContactState<:IpState
    ctx::Context
    f ::Float64
    w ::Float64
    function LinearTipContactState(ctx::Context)
        this = new(ctx)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct LinearTipContact<:Material
    k::Float64

    function LinearTipContact(prms::Dict{Symbol,Float64})
        return  LinearTipContact(;prms...)
    end

    function LinearTipContact(;k=NaN)
        @check k>=0
        this = new(k)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearTipContact}) = LinearTipContactState

# Element types that work with this material
compat_elem_types(::Type{LinearTipContact}) = (MechTipJoint,)


function calcD(mat::LinearTipContact, state::LinearTipContactState)
    return mat.k
end


function update_state!(mat::LinearTipContact, state::LinearTipContactState, Δw)
    Δf = mat.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function ip_state_vals(mat::LinearTipContact, state::LinearTipContactState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
