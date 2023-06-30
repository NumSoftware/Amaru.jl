# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechProperties

struct MechProperties<:Properties
    ρ::Float16
    A::Float64
    th::Float64

    function MechProperties(;rho=0.0, A=0.0, th=0.0)
        return new(rho, A, th)
    end
end
