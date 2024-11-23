# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type Analysis 
    # stages::Array{Stage,1}
    # loggers::Array{AbstractLogger,1}
    # monitors::Array{AbstractMonitor,1}
end

abstract type TransientAnalysis<:Analysis end
abstract type StaticAnalysis<:Analysis end


export addstage!

addstage!_params = [
    FunInfo(:addstage!, "Add a stage to a FEModel instance"),
    ArgInfo(:analysis, "A FE analysis instance"),
    ArgInfo(:bcs, "Array of boundary conditions"),
    KwArgInfo(:nincs, "Number of increments", 1),
    KwArgInfo(:nouts, "Number of output files", 0),
    KwArgInfo(:tspan, "Time span", 0.0),
    KwArgInfo(:toactivate, "Array of elements to activate", Element[]),
    KwArgInfo(:todeactivate, "Array of elements to deactivate", Element[]),
]
@doc docstring(addstage!_params) addstage!

function addstage!(analysis::Analysis, bcs::AbstractArray; kwargs...)
    args = checkargs([analysis, bcs], kwargs, addstage!_params)
    stage = Stage(args.bcs; args.nincs, args.nouts, args.tspan, args.toactivate, args.todeactivate)
    stage.id = length(analysis.stages) + 1
    push!(analysis.stages, stage)
end


function addlogger!(analysis::Analysis, logpair::Tuple)
    push!(analysis.loggers, logpair[2])
    setup_logger!(analysis, logpair[1], logpair[2])
end

"""
    $(TYPEDSIGNATURES)

Register the loggers from the array `loggers` into `analysis`.

"""
function addloggers!(analysis::Analysis, loggers::Array{<:Tuple,1})
    analysis.loggers = []
    for (filter,logger) in loggers
        push!(analysis.loggers, logger)
        setup_logger!(analysis, filter, logger)
    end
end
setloggers! = addloggers!

function addmonitor!(analysis::Analysis, monpair::Pair)
    push!(analysis.monitors, monpair.second)
    setup_monitor!(analysis, monpair.first, monpair.second)
end

"""
    addmonitors!(analysis, loggers)

Register monitors from the array `loggers` into `analysis`.

"""
function addmonitors!(analysis::Analysis, monitors::Vector{<:Pair})
    analysis.monitors = []
    for (filter,monitor) in monitors
        push!(analysis.monitors, monitor)
        setup_monitor!(analysis, filter, monitor)
    end
end
setmonitors! = addmonitors!


function update_records!(ana::Analysis; checkpoint=true, force=false)
    sctx = ana.sctx
    outdir = sctx.outdir

    flushinterval = 5.0
    flush = time()-sctx.flushtime>flushinterval || checkpoint || force || sctx.T >= 1-1e-8
    flush && (sctx.flushtime = time())

    if checkpoint
        rm.(glob("*conflicted*.log"), force=true)
        rm.(glob("*conflicted*.*", "$outdir/"), force=true)
                
        update_output_data!(ana.model) # need to be before group loggers
        save(ana.model, "$outdir/$(sctx.outkey)-$(sctx.out).vtu", quiet=true)

        # update multiloggers
        for logger in ana.loggers
            if isa(logger, MultiLogger)
                update_logger!(logger, ana)
                logger.filename!="" && save(logger.book, logger.filename, quiet=true)
            end
        end

    end

    flush && Base.flush(sctx.log)

    # update single loggers
    sctx.Tupdate > sctx.T && (sctx.Tupdate=0.0) # for subsequent stages
    update_single_loggers = sctx.T-sctx.Tupdate >= 0.00025 || sctx.T==0
    update_single_loggers && (sctx.Tupdate = sctx.T)
    
    for logger in ana.loggers
        if isa(logger, SingleLogger)
            update_single_loggers && update_logger!(logger, ana)
            flush && logger.filename!="" && save(logger.table, logger.filename, quiet=true)
        end
    end

    # update monitors
    for monitor in ana.monitors
        rstatus = update_monitor!(monitor, ana)
        failed(rstatus) && return rstatus
        flush && monitor.filename!="" && save(monitor.table, monitor.filename, quiet=true)
    end

    return success()
end