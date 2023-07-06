# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechProperties

struct MechProperties<:Properties
    ρ::Float16
    A::Float64
    th::Float64

    function MechProperties(;rho=NaN, A=NaN, th=NaN)
        !isnan(rho) && @check rho>0
        !isnan(A  ) && @check A>0
        !isnan(th ) && @check th>0

        return new(rho, A, th)
    end
end

# newproperties(::MechModel) = MechProperties()
