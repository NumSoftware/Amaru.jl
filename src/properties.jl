# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type ElemProperties end

@inline Base.:(<<)(a, b::ElemProperties) = return (a, b)
@inline Base.:(=>)(a, b::ElemProperties) = return (a, b)