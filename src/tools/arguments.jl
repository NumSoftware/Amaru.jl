# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct ArgInfo
    key::Symbol
    default::Any
    cond::Expr
    desc::String
end

mutable struct ArgCond
    cond::Expr
end

mutable struct ArgOpt
    args1::Union{Symbol,Expr}
    args2::Union{Symbol,Expr}
end

macro arginfo(args...)
    if typeof(args[1])==Symbol
        key = args[1]
        default = nothing
    else # arg=value
        key = args[1].args[1]
        default = args[1].args[2]
    end
    cond = args[2]
    desc = args[end]

    return quote
        ArgInfo( $(Meta.quot(key)), $default, $(Meta.quot(cond)), $desc ) 
    end
end

macro argcond(cond)
    return quote
        ArgCond( $(Meta.quot(cond)) ) 
    end
end

macro argopt(args1, args2)
    return quote
        ArgOpt( $(Meta.quot(args1)), $(Meta.quot(args2)) ) 
    end
end

function checkargs(args, args_rules::Array)
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

    # list of optional keys to skip
    argopts  = [ item for item in args_rules if item isa ArgOpt ]
    skipkeys = []
    for opt in argopts
        args1 = opt.args1 isa Symbol ? (opt.args1,) : opt.args1.args
        args2 = opt.args2 isa Symbol ? (opt.args2,) : opt.args2.args
        if all(key in argkeys for key in args1)
            append!(skipkeys, args2)
        end
        if all(key in argkeys for key in args2)
            append!(skipkeys, args1)
        end
    end

    arginfos = [ item for item in args_rules if item isa ArgInfo && !(item.key in skipkeys) ]
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
        arginfo.cond == () && continue
        key = arginfo.key
        val = args[key]
        if !eval_arith_expr(arginfo.cond; (; arginfo.key => val)...)
            msg = "Invalid value for argument $key which must satisfy $(arginfo.cond). Got $key = $val"
            throw(AmaruException("$fname: $msg"))
        end
    end
    
    # check relational conditions
    for argcond in allconds
        if !eval_arith_expr(argcond.cond, args...)
            msg = "Given arguments do not satisfy condition $(argcond.cond)"
            throw(AmaruException("$fname: $msg"))
        end
    end

    return args
end