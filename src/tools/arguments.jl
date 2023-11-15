# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Returns a list containing the arguments for each method of function
typeofargs(f) = [ ((t for t in fieldtypes(m.sig)[2:end])...,) for m in methods(f) ]


abstract type ArgObj
end

mutable struct FunInfo<:ArgObj
    key::Symbol
    desc::String
    signature::Tuple
end

mutable struct ArgInfo<:ArgObj
    key::Symbol
    aliases::Tuple
    default::Any
    cond::Expr
    values::Any
    length::Int
    type::Union{DataType, UnionAll}
    desc::String
end

function ArgInfo(key, desc, default=missing; condition=:(), values=(), length=0, type=Any)
    if key isa Tuple
        aliases = key[2:end]
        key = key[1]
    else
        aliases = ()
    end

    return ArgInfo(key, aliases, default, condition, values, length, type, desc)
end

mutable struct ArgCond<:ArgObj
    cond::Expr
end

mutable struct ArgOpt<:ArgObj
    args1::Union{Symbol,Tuple}
    args2::Union{Symbol,Tuple}
end

# todo: to be deprecated 
macro arginfo(args...)
    arg1 = args[1]
    if arg1 isa Symbol
        key = args[1]
        aliases = ()
        default = missing
    else # arg=value
        key = args[1].args[1]
        if key isa Tuple # has aliases
            key = key[1]
            aliases = key[2:end]
        else
            aliases = ()
        end

        default = args[1].args[2]
    end

    arg2 = args[2]
    if arg2 isa Tuple # has possible values instead of condition
        cond = :()
        values = arg2
    else
        cond = arg2
        values = ()
    end

    desc = args[end]

    return quote
        ArgInfo( $(Meta.quot(key)), $(Meta.quot(aliases)), $default, $(Meta.quot(cond)), $values, 0, Any, $desc ) 
    end
end

# todo: to be deprecated 
macro argcond(cond)
    return quote
        ArgCond( $(Meta.quot(cond)) ) 
    end
end

# todo: to be deprecated 
macro argopt(args1, args2)
    return quote
        ArgOpt( $(Meta.quot(args1)), $(Meta.quot(args2)) ) 
    end
end

function get_args_descriptions(args_params::AbstractArray)
    descs = String[]
    for item in args_params
        item isa ArgInfo || continue
        desc = "$(item.key) "
        if length(item.aliases)>0
            desc = desc * "($(join(string.(item.aliases), ", "))) "
        end
        desc = desc*": $(item.desc)"
        push!(descs, desc)
    end
    return join(descs, "\n")
end


function checkargs(args, args_params::AbstractArray; aliens=true)
    # Get function caller
    st = stacktrace(backtrace())
    fname = :_
    for frame in st
        if !startswith(string(frame.func), "_") && !contains(string(frame.func), "checkargs")
            fname = frame.func
            break
        end
    end

    # check for wrong arguments
    if !aliens
        allkeys = []
        for item in args_params 
            if item isa ArgInfo
                push!(allkeys, item.key)
                append!(allkeys, item.aliases)
            end
        end

        argkeys = keys(args)
        for key in argkeys
            if !(key in allkeys)
                msg = "Wrong argument $(key). The named arguments are:\n"*get_args_descriptions(args_params)
                throw(AmaruException(msg))     
            end
        end
    end

    # replace aliases (update args)
    pairs = []
    for (key, value) in args
        hasalias = false
        newkey = key
        for item in args_params
            item isa ArgInfo || continue
            if key in item.aliases
                hasalias = true
                newkey = item.key
                break
            end
        end
        push!(pairs, newkey=>value)
    end
    args = NamedTuple(pairs) # redefine args replacing aliases
    argkeys = keys(args)

    # list of optional keys to skip
    argopts  = [ item for item in args_params if item isa ArgOpt ]
    skipkeys = []
    for opt in argopts
        # check for more than one optionals per set
        args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1
        args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2

        if issubset(args1, argkeys) && issubset(args2, argkeys)
            s1 = join(args1, ", ")
            s2 = join(args2, ", ")
            throw(AmaruException("$fname: Invalid argument combination. Provide either argument(s) ($s1) or ($s2)"))
        end

        if all(key in argkeys for key in args1)
            append!(skipkeys, args2)
        end
        if all(key in argkeys for key in args2)
            append!(skipkeys, args1)
        end
    end

    arginfos = [ item for item in args_params if item isa ArgInfo && !(item.key in skipkeys) ] # skip non used optional keys

    # check for missing keys (skip optional args)
    mandatorykeys = [ item.key for item in arginfos if item.default===missing ]
    missingkeys = setdiff(mandatorykeys, keys(args))

    if length(missingkeys)>0
        msg = "Missing arguments: $(join(missingkeys, ", ")). The named arguments are:\n"*get_args_descriptions(args_params)
        for opt in argopts
            args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1.args
            args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2.args
            msg = msg*"\nArgument(s) $(join(args1, ", ", " and ")) can be used instead of $(join(args2, ", ", " and "))."
        end
        throw(AmaruException("$fname: $msg"))
    end

    # update args with defaults
    pairs = []
    allkeys = [ item.key for item in arginfos ]
    defaults  = NamedTuple([ item.key => item.default for item in arginfos ])
    args_keys = keys(args)
    for item in arginfos
        key = item.key
        if key in args_keys
            value = args[key] # use given value
        else
            default = item.default
            if default isa Symbol && default in allkeys # default value in terms of other values
                if default in args_keys
                    value = args[default]
                else
                    value = defaults[default]
                end
            else
                value = default # use default value
            end
        end
        push!(pairs, key=>value)
    end

    args = NamedTuple(pairs)

    # allkeys = [ item.key for item in arginfos ]
    # @show keys(args)
    # for item in arginfos
    #     @show item.key
    #     @show item.default
    #     if item.default isa Symbol && item.default in allkeys
    #         @show "xxxxxxxxxxxxxxxx"
    #         @show item.key
    #         item.default = getindex(args, item.default, nothing)
    #     end
    # end

    # alldefaults  = NamedTuple([ item.key => item.default for item in arginfos ])

    # @show alldefaults

    # args = merge(alldefaults, args)

    # check argument condition
    for arginfo in arginfos
        key = arginfo.key
        val = args[key]

        # check condition
        if arginfo.cond != :()
            func = Expr(:(->), key, arginfo.cond)
            if !invokelatest(eval(func), val)
                msg = "Invalid value for argument $key which must satisfy $(arginfo.cond).\nGot $key = $(repr(val))"
                throw(AmaruException("$fname: $msg"))
            end
        end

        # check if in the give set
        if length(arginfo.values)>0 && !(val in arginfo.values)
            msg = "Invalid value for argument $key which must be one of $(arginfo.values).\nGot $key = $(repr(val))"
            throw(AmaruException("$fname: $msg"))
        end

        # check length
        if arginfo.length>0
            haslength = hasmethod(length, (typeof(val),))
            rightlength = haslength && length(val)==arginfo.length
            if !haslength || !rightlength
                msg = "Invalid value for argument $key which must be of size $(arginfo.length).\nGot $key = $(repr(val))"
                throw(AmaruException("$fname: $msg"))
            end
        end

        # check type
        if arginfo.type!=Nothing && !(val isa arginfo.type)
            msg = "Invalid value for argument $key which must be of type $(arginfo.type).\nGot $key = $(repr(val))"
            throw(AmaruException("$fname: $msg"))
        end
    end
    
    # check relational conditions between arguments
    allconds = [ item for item in args_params if item isa ArgCond ]
    for argcond in allconds
        if !eval_arith_expr(argcond.cond, args...) # todo: update with eval
            msg = "Given arguments do not satisfy condition $(argcond.cond)"
            throw(AmaruException("$fname: $msg"))
        end
    end

    return args
end


function make_doc(typ)
    fparams = func_params(typ)
    
    # find description
    idx = findfirst(x->isa(x,FunInfo), fparams)
    fdesc = fparams[idx]

    # add function description
    desc = String[]
    push!(desc, "`$(fdesc.key)$(string(fdesc.signature))`\n")
    push!(desc, "$(fdesc.desc)")
    push!(desc, "# Optional arguments:\n")

    # add arguments description
    for item in fparams
        item isa ArgInfo || continue
        str = "- `$(item.key)` "
        
        # aliases
        if length(item.aliases)>0
            str_aliases = [ "`"*string(alias)*"`" for alias in item.aliases ]
            str = str*"("*join(str_aliases, ", ")*")"
        end
        str = str*": $(item.desc)"

        # suitable values
        if length(item.values)>0
            str = str*", any of: ("*join(repr.(item.values), ", ")*")"
        end

        push!(desc, str)
    end

    return join(desc, "\n")
end

