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

    # check mandatory keys
    arginfos = [ item for item in args_rules if item isa ArgInfo ]

    alldescs = join( ["$(item.key): $(item.desc)" for item in arginfos], ", ")

    mandatorykeys = [ item.key for item in arginfos if item.default===nothing ]
    missingkeys = setdiff(mandatorykeys, argkeys)

    # @show mandatorykeys
    # @show missingkeys

    if length(missingkeys)>0
        msg = "Missing arguments found: $(join(missingkeys, ", ")). Possible inputs are: $alldescs"
        throw(AmaruException("$fname: $msg"))
    end

    # check consistency of optional arguments
    # todo

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
    
    # check independent conditions
    allconds = [ item for item in args_rules if item isa ArgCond ]
    
    for argcond in allconds
        if !eval_arith_expr(argcond.cond, args...)
            msg = "Given arguments do not satisfy condition $(argcond.cond)"
            throw(AmaruException("$fname: $msg"))
        end
    end

    return args
end