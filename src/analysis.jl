# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type Analysis end

abstract type TransientAnalysis<:Analysis end
abstract type StaticAnalysis<:Analysis end