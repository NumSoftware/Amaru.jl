# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Returns a list containing the arguments for each method of function
typeofargs(f) = [ ((t for t in fieldtypes(m.sig)[2:end])...,) for m in methods(f) ]


abstract type ArgObj
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

function ArgInfo(key, desc; default=nothing, condition=:(), values=(), length=0, type=Any)
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
    args1::Union{Symbol,Expr}
    args2::Union{Symbol,Expr}
end

# todo: to be deprecated 
macro arginfo(args...)
    arg1 = args[1]
    if arg1 isa Symbol
        key = args[1]
        aliases = ()
        default = nothing
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

function checkargs(args, rules...)
    return checkargs(args, [rules...])
end

function checkargs(args, args_rules::AbstractArray)
     # Get function caller
     st = stacktrace(backtrace())
     fname = :_
     for frame in st
         if !startswith(string(frame.func), "_") && frame.func!=Symbol("checkargs")
             fname = frame.func
             break
         end
     end

    argkeys = keys(args)

    # replace aliases (update args)
    pairs = []
    for (key, value) in args
        hasalias = false
        newkey = key
        for item in args_rules
            item isa ArgInfo || continue
            if key in item.aliases
                hasalias = true
                newkey = item.key
                break
            end
        end
        push!(pairs, newkey=>value)
    end
    args = NamedTuple(pairs)

    # list of optional keys to skip
    argopts  = [ item for item in args_rules if item isa ArgOpt ]
    skipkeys = []
    for opt in argopts
        # check for more than one optionals per set
        args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1.args
        args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2.args
        if all(key in argkeys for key in args1)
            append!(skipkeys, args2)
        end
        if all(key in argkeys for key in args2)
            append!(skipkeys, args1)
        end
    end

    arginfos = [ item for item in args_rules if item isa ArgInfo && !(item.key in skipkeys) ] # skip non used optional keys
    alldescs = join( ["$(item.key): $(item.desc)" for item in arginfos], ", ", " and ")
    allconds = [ item for item in args_rules if item isa ArgCond ]

    # check mandatory keys (skip optional args)
    mandatorykeys = [ item.key for item in arginfos if item.default===nothing ]
    missingkeys = setdiff(mandatorykeys, argkeys)

    if length(missingkeys)>0
        msg = "Missing arguments: $(join(missingkeys, ", ")).\nPossible inputs are: $alldescs"
        for opt in argopts
            args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1.args
            args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2.args
            msg = msg*"\nArgument(s) $(join(args1, ", ", " and ")) can be used instead of $(join(args2, ", ", " and "))."
        end
        throw(AmaruException("$fname: $msg"))
    end

    # update args with defaults
    alldefaults  = NamedTuple([ item.key => item.default for item in arginfos ])
    args = merge(alldefaults, args)

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
    for argcond in allconds
        if !eval_arith_expr(argcond.cond, args...)
            msg = "Given arguments do not satisfy condition $(argcond.cond)"
            throw(AmaruException("$fname: $msg"))
        end
    end

    return args
end