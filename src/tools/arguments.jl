# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Returns a list containing the arguments for each method of function
typeofargs(f) = [ ((t for t in fieldtypes(m.sig)[2:end])...,) for m in methods(f) ]


abstract type ArgObj
end

mutable struct FunInfo<:ArgObj
    key::Symbol
    desc::String
end


mutable struct ArgInfo<:ArgObj
    key::Symbol
    default::Any
    cond::Expr
    values::Any
    length::Int
    type::Union{DataType, Union, UnionAll}
    desc::String
    complement::String
end


function ArgInfo(key, desc, default=missing; cond=:(), values=(), length=0, type=Any, complement="")
    return ArgInfo(key, default, cond, values, length, type, desc, complement)
end


mutable struct KwArgInfo<:ArgObj
    key::Symbol
    aliases::Tuple
    default::Any
    cond::Expr
    values::Any
    length::Int
    type::Union{DataType, Union, UnionAll}
    desc::String
    complement::String
end

function KwArgInfo(key, desc, default=missing; cond=:(), values=(), length=0, type=Any, complement="")
    if key isa Tuple
        aliases = key[2:end]
        key = key[1]
    else
        aliases = ()
    end

    return KwArgInfo(key, aliases, default, cond, values, length, type, desc, complement)
end

mutable struct ArgCond<:ArgObj
    cond::Expr
    message::String

    function ArgCond(cond, message="")
        return new(cond, message)
    end
end

mutable struct ArgOpt<:ArgObj
    args1::Union{Symbol,Tuple}
    args2::Union{Symbol,Tuple}
end


function get_args_descriptions(fparams::AbstractArray)
    descs = String[]
    for item in fparams
        item isa KwArgInfo || continue
        desc = "$(item.key) "
        if length(item.aliases)>0
            desc = desc * "($(join(string.(item.aliases), ", "))) "
        end
        desc = desc*": $(item.desc)"
        push!(descs, desc)
    end
    return join(descs, "\n")
end


# Checks positional and keyword arguments for a function
function checkargs(args::AbstractArray, kwargs, fparams::AbstractArray; aliens=true, check_ps=true)
    PsInfos = ArgInfo[ arg for arg in fparams if arg isa ArgInfo ]
    KwInfos = KwArgInfo[ arg for arg in fparams if arg isa KwArgInfo ]
    Conds = ArgCond[ arg for arg in fparams if arg isa ArgCond ]
    Opts  = ArgOpt[ arg for arg in fparams if arg isa ArgOpt ]

    # Get function caller
    function funcname()
        st = stacktrace(backtrace())
        for frame in st
            str = string(frame.func)
            if !contains(str, r"(^_|funcname|checkargs|macro)")
                return match(r"#?(\w+)#?",  str)[1] # remove #s
            end
        end
        return "_"
    end

    # number of positional arguments
    n = length(args) # all args
    ninfos = length(PsInfos)
    check_ps = check_ps || ninfos==0
    
    if check_ps
        # check number of positional arguments
        if n != ninfos
            msg = "Wrong number of positional arguments. Expected $n, got $ninfos."
            throw(AmaruException("$(funcname()): $msg"))
        end

        # get a list of key=>val pairs of positional arguments
        ps_pairs = []
        for (i, item) in enumerate(PsInfos)
            key = item.key
            val = args[i]
            push!(ps_pairs, key=>val)
        end
    end

    # check for wrong arguments in kwargs
    if !aliens
        allkeys = Symbol[]
        for item in KwInfos
            push!(allkeys, item.key)
            append!(allkeys, item.aliases)
        end

        for key in keys(kwargs)
            if !(key in allkeys)
                msg = "Wrong argument $(key). The named arguments are:\n"*get_args_descriptions(fparams)
                throw(AmaruException("$(funcname()): $msg"))     
            end
        end
    end

    # replace aliases (update kwargs)
    pairs = []
    for (key, value) in kwargs
        newkey = key
        for item in KwInfos
            if key in item.aliases
                newkey = item.key
                break
            end
        end
        push!(pairs, newkey=>value)
    end
    kwargs = NamedTuple(pairs) # redefine kwargs replacing aliases
    argkeys = keys(kwargs)

    # list of optional keys to skip in kwargs
    skipkeys = Symbol[]
    for opt in Opts
        # check for more than one optionals per set
        args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1
        args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2

        if issubset(args1, argkeys) && issubset(args2, argkeys)
            s1 = join(args1, ", ")
            s2 = join(args2, ", ")
            throw(AmaruException("$(funcname()): Invalid argument combination. Provide either argument(s) ($s1) or ($s2)"))
        end

        if all(key in argkeys for key in args1)
            append!(skipkeys, args2)
        end
        if all(key in argkeys for key in args2)
            append!(skipkeys, args1)
        end
    end

    KwInfos = [ item for item in KwInfos if !(item.key in skipkeys) ] # skip non used optional keys

    # check for missing keys (skip optional args)
    mandatorykeys = Symbol[ item.key for item in KwInfos if ismissing(item.default) ]
    missingkeys = setdiff(mandatorykeys, keys(kwargs))

    if length(missingkeys)>0
        msg = "Missing arguments: $(join(missingkeys, ", ")). The named arguments are:\n"*get_args_descriptions(fparams)
        for opt in Opts
            args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1.args
            args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2.args
            msg = msg*"\nArgument(s) $(join(args1, ", ", " and ")) can be used instead of $(join(args2, ", ", " and "))."
        end
        throw(AmaruException("$(funcname()): $msg"))
    end

    # update kwargs with defaults args
    kw_pairs = []
    allkeys = Symbol[ item.key for item in KwInfos ]
    defaults  = NamedTuple([ item.key => item.default for item in KwInfos ])
    args_keys = keys(kwargs)
    for item in KwInfos
        key = item.key
        if key in args_keys
            value = kwargs[key] # use given value
        else
            default = item.default
            if default isa Symbol && default in allkeys # default value in terms of other values
                if default in args_keys
                    value = kwargs[default]
                else
                    value = defaults[default]
                end
            else
                value = default # use default value
            end
        end
        push!(kw_pairs, key=>value)
    end

    # get all arguments
    if check_ps
        args = NamedTuple([ps_pairs; kw_pairs])
        ArgInfos = [PsInfos; KwInfos]
    else
        args = NamedTuple(kw_pairs)
        ArgInfos = KwInfos
    end

    # check argument condition in both positional and keyword arguments
    for arginfo in ArgInfos
        key = arginfo.key
        val = args[key]

        # check condition
        if arginfo.cond != :() && !isnothing(val)
            func = Expr(:(->), key, arginfo.cond)
            if !invokelatest(eval(func), val)
                msg = "Invalid value for argument $key which must satisfy $(arginfo.cond).\nGot $key = $(repr(val))"
                throw(AmaruException("$(funcname()): $msg"))
            end
        end

        # check if in the given set
        if length(arginfo.values)>0 && !(val in arginfo.values)
            msg = "Invalid value for argument $key which must be one of $(arginfo.values).\nGot $key = $(repr(val))"
            throw(AmaruException("$(funcname()): $msg"))
        end

        # check length
        if arginfo.length>0
            haslength = hasmethod(length, (typeof(val),))
            rightlength = haslength && length(val)==arginfo.length
            if !haslength || !rightlength
                msg = "Invalid value for argument $key which must be of size $(arginfo.length).\nGot $key = $(repr(val))"
                throw(AmaruException("$(funcname()): $msg"))
            end
        end

        # check type
        if arginfo.type !== nothing && !(val isa arginfo.type)
            msg = "Invalid value for argument $key which must be of type $(arginfo.type).\nGot $key = $(repr(val))"
            throw(AmaruException("$(funcname()): $msg"))
        end
    end
    
    # check relational conditions between all arguments
    allconds = [ item for item in fparams if item isa ArgCond ]
    for argcond in allconds
        if !evaluate(argcond.cond; args...) # todo: update with eval
            msg = argcond.message!="" ? argcond.message : "Given arguments do not satisfy condition $(argcond.cond)"
            throw(AmaruException("$(funcname()): $msg"))
        end
    end

    return args
end


function checkargs(kwargs, fparams::AbstractArray; aliens=true)
    return checkargs([], kwargs, fparams::AbstractArray; aliens=aliens, check_ps=false)
end


function docstring(fparams)
    PsInfos = ArgInfo[ arg for arg in fparams if arg isa ArgInfo ]
    KwInfos = KwArgInfo[ arg for arg in fparams if arg isa KwArgInfo ]
    
    # find description
    idx = findfirst(x->isa(x,FunInfo), fparams)
    fdesc = fparams[idx]

    # add function description
    desc = String[]
    signature = "    $(fdesc.key)("
    if length(PsInfos)>0
        signature *= join([ "$(item.key)" for item in PsInfos ], ", ")
    end
    if length(KwInfos)>0
        signature *= "; kwargs..."
    end
    signature *= ")\n"
    push!(desc, signature)


    push!(desc, "$(fdesc.desc)")
    
    # add positional arguments descriptions
    if length(PsInfos)>0
        push!(desc, "# Positional arguments:\n")
        for item in PsInfos
            str = "- `$(item.key)` "
            str = str*": $(item.desc)."
            
            # condition
            if item.cond != :()
                str = str*" Requires $(item.cond)."
            end

            # suitable values
            if length(item.values)>0
                str = str*" Any of "*join(repr.(item.values), ", ", " and ")*"."
            end

            push!(desc, str)
        end
    end
    
    # add keyword arguments descriptions
    if length(KwInfos)>0
        push!(desc, "# Keyword arguments:\n")
        for item in KwInfos
            str = "- `$(item.key)` "
            
            # aliases
            if length(item.aliases)>0
                str_aliases = [ "`"*string(alias)*"`" for alias in item.aliases ]
                str = str*"(or "*join(str_aliases, ", ")*")."
            end
            str = str*": $(item.desc)."

            # condition
            if item.cond != :()
                str = str*" Requires $(item.cond)."
            end

            # suitable values
            if length(item.values)>0
                str = str*" Any of "*join(repr.(item.values), ", ", " and ")*"."
            end

            # default
            if !ismissing(item.default)
                str = str*" Default is $(repr(item.default))."
            end

            # complement
            if item.complement != ""
                str = str*" $(item.complement)."
            end

            push!(desc, str)
        end
    end

    return join(desc, "\n")
end

