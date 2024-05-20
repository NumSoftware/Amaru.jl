# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

include("returnstatus.jl")
include("iteration.jl")

include("color.jl")
include("show.jl")

include("constants.jl")
include("math.jl")
include("linalg.jl")
include("path-function.jl")
include("tensors.jl")

include("signal.jl")
include("quaternion.jl")

include("expr.jl")
include("arguments.jl")
include("threads.jl")
include("table.jl")
include("book.jl")
include("utils.jl")
include("stopwatch.jl")
include("xml.jl")

include("array-pool.jl")

include("tex.jl")
include("error.jl")

Base.show(io::IO, obj::XmlElement) = _show(io, obj, 3, "")
Base.show(io::IO, obj::XmlDocument)  = _show(io, obj, 3, "")
