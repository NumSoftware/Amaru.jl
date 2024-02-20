# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Solves a system with unknowns in U and F vectors
function solve_system!(
                       K ::SparseMatrixCSC{Float64, Int},
                       U ::Vect,
                       F ::Vect,
                       nu::Int,
                      )
    #  ┌  K11   K12 ┐  ┌ U1? ┐    ┌ F1  ┐
    #  │            │  │     │ =  │     │
    #  └  K21   K22 ┘  └ U2  ┘    └ F2? ┘

    msg = ""

    # Decomposing the coefficients matrix
    if nu>0
        nu1 = nu+1
        K11 = K[1:nu, 1:nu]
        K12 = K[1:nu, nu1:end]
        K21 = K[nu1:end, 1:nu]
    end
    K22 = K[nu+1:end, nu+1:end]

    F1  = F[1:nu]
    U2  = U[nu+1:end]

    # @showm K11
    
    # Solve linear system
    F2 = K22*U2
    U1 = zeros(nu)
    # @showm F1
    # @showm F2
    # @showm U2
    if nu>0
        RHS = F1 - K12*U2

        try
            # try
                LUfact = lu(K11)
                U1 = LUfact\RHS
            # catch err
            #     err isa InterruptException && rethrow(err)
            #     if typeof(err)==SingularException
            #         # Regularization attempt
            #         msg = "$msg\nsolve_system!: Syngular matrix - regularization attempt"
            #         S = spdiagm([ 1/maximum(abs, K11[i,:]) for i in 1:nu ])
            #         LUfact = lu(S*K11)
            #         U1  = (LUfact\(S*RHS))
            #     else
            #         return failure("$msg\nsolve_system!: $err")
            #     end
            # end

            F2 += K21*U1
        catch err
            err isa InterruptException && rethrow(err)
            if any(isnan.(K11)) 
                msg = "$msg\nsolve_system!: NaN values in coefficients matrix"
            end
            # U1 .= NaN
            return failure("$msg\nsolve_system!: $err")
        end
    end

    maxU = 1e8 # maximum essential value
    if maximum(abs, U1)>maxU 
        return failure("$msg\nsolve_system!: Possible syngular matrix")
    end

    # Completing vectors
    U[1:nu]     .= U1
    F[nu+1:end] .= F2

    yield()
    return success(msg)
end


function stage_iterator!(name::String, stage_solver!::Function, model::Model; args...)
    autoinc = get(args, :autoinc, false)
    outdir  = get(args, :outdir, ".")
    outkey  = get(args, :outkey, "out")
    quiet   = get(args, :quiet, false)
    env     = model.env
    
    cstage = findfirst(st->st.status!=:done, model.stages)
    cstage === nothing && throw(AmaruException("stage_iterator!: No stages have been set for $name"))

    solstatus = success()

    if !quiet && cstage==1 
        printstyled(name, "\n", bold=true, color=:cyan)
        println("  active threads: ", Threads.nthreads())
        # println("  model type: ", env.ana.stressmodel)
    end

    outdir = rstrip(outdir, ['/', '\\'])
    env.outdir = outdir
    env.outkey = outkey

    if !isdir(outdir)
        info("solve!: creating output directory ./$outdir")
        mkpath(outdir)
    end


    if cstage==1
        env.log = open("$outdir/solve.log", "w")
    else
        env.log = open("$outdir/solve.log", "a")
    end

    
    for stage in model.stages[cstage:end]
        stage.status = :solving

        nincs  = stage.nincs
        nouts  = stage.nouts
        
        env.stage = stage.id
        env.inc   = 0
        env.T = 0.0

        if !quiet
            # if stage.id>1 || length(model.stages)>1
                printstyled("Stage $(stage.id)\n", bold=true, color=:cyan)
            # end
        end

        save_outs = stage.nouts > 0
        if save_outs
            if !autoinc
                if nouts > nincs
                    nincs = nouts
                    quiet || info("nincs changed to $(nincs) to match nouts")
                end
                if nincs%nouts != 0
                    stage.nincs = nincs - (nincs%nouts) + nouts
                    quiet || info("nincs changed to $nincs to be multiple of nouts")
                end
            end
            stage.nincs = nincs
            stage.nouts = nouts
        end

        sw = StopWatch() # timing
        # sline = StatusLine(model, stage, sw)
        if !quiet
            status_cycler_task = Threads.@spawn :interactive status_cycler(model, sw)
        end

        local runerror
        local error_st
        try
            solstatus = stage_solver!(model, stage; args...)
            if succeeded(solstatus)
                stage.status = :done
            else
                stage.status = :failed
            end
        catch err            
            runerror = err
            flush(env.log)
            if err isa InterruptException
                stage.status = :interrupted
            else
                stage.status = :error
                error_st = stacktrace(catch_backtrace())
            end
        end
        close(env.log)

        if !quiet
            wait(status_cycler_task)
            solstatus.message != "" && println(solstatus.message)
        end

        if stage.status == :interrupted 
            throw(AmaruException("The analysis was interrupted"))
        elseif stage.status == :error
            # trim not important frames; search for the frame that contains current function
            idx = findfirst(contains("iterator"), string(frame) for frame in error_st)
            if idx!==nothing
                error_st = error_st[1:idx-1]
            end

            alert("Amaru internal error", level=1)
            showerror(stdout, runerror, error_st)
            println()
            # error()
            # @show runerror.backtrace()
            throw(runerror) 
        end

        getlapse(sw)>60 && sound_alert()

    end
    return solstatus

end

# Main function to call specific solvers
function solve!(model::FEModel; args...)
    solve!(model, model.env.ana; args...)
end


function progress_bar(T::Float64)
    dwidth  = displaysize(stdout)[2]-2
    width   = max(2, min(25, dwidth-35))
    ch_done = T*width
    frac    = ch_done - floor(ch_done)

    barl = repeat(['━'], floor(Int, ch_done))
    barr = Char[]

    if frac<=0.25
        # @show 10
        push!(barr, '╶')
        push!(barr, '─')
    elseif 0.25<frac<0.75
        push!(barl, '╸')
        push!(barr, '─')
    else
        push!(barl, '━')
        push!(barr, '╶')
    end

    if length(barl)>=width
        barl = barl[1:width]
        barr = Char[]
    end
    
    if length(barl)+length(barr)==width+1
        barr = Char[barr[1]]
    end

    append!(barr, repeat(['─'], width -length(barl) -length(barr) ))
    barls = reduce(*, barl) 
    barrs = reduce(*, barr) 

    iscolor = get(stdout, :color, false)
    if iscolor
        color        = :blue
        enable_color = get(Base.text_colors, color, Base.text_colors[:default])
        enable_bold  = get(Base.text_colors, :bold, Base.text_colors[:default])
        normal_color  = get(Base.disable_text_style, :normal, Base.text_colors[:default])
        disable_bold = get(Base.disable_text_style, :bold, Base.text_colors[:default])
        barls        = string(enable_color, enable_bold, barls, disable_bold, normal_color)
        enable_color = get(Base.text_colors, :light_black, Base.text_colors[:default])
        normal_color = get(Base.disable_text_style, :bold, Base.text_colors[:default])
        barrs        = string(enable_color, barrs, normal_color)
    end

    return barls*barrs

end


function status_cycler(model::FEModel, sw::StopWatch)
    print("\e[?25l") # disable cursor

    stage     = model.stages[model.env.stage]
    last_loop = false
    alerts    = String[]
    while true
        nlines = 0

        nlines += print_info(model)
        nlines += print_alerts(model, alerts)
        nlines += print_summary(model, sw)
        
        last_loop && break

        print("\e[$(nlines)A")
        stage.status != :solving && (last_loop=true)
        sleep(0.05)
        # yield()
    end

    print("\e[?25h") # enable cursor
end


function print_info(model::FEModel)
    str = strip(String(take!(model.env.info)))
    str!="" && println("  ", str, "\e[K")

    return 0
end


function print_alerts(model::FEModel, alerts::Array{String,1})
    str = strip(String(take!(model.env.alerts)))
    
    if str!=""
        list = split(str, "\n")
        list = String[ string("  ", Time(now()), "  ", m) for m in list ]
        append!(alerts, list)
    end
    
    n = length(alerts)
    if n>5
        splice!(alerts, 1:n-5)
        alerts[1] = "  ⋮"
    end

    for m in alerts
        printstyled(m, "\e[K\n", color=Base.warn_color())
    end

    return length(alerts)
end


function print_summary(model::FEModel, sw::StopWatch)
    env = model.env
    nlines = 2

    # line 1:
    T  = env.T
    ΔT = env.ΔT
    printstyled("  inc $(env.inc) output $(env.out)", bold=true, color=:light_blue)
    if env.transient
        t = round(env.t, sigdigits=3)
        printstyled(" t=$t", bold=true, color=:light_blue)
    end
    dT  = round(ΔT,sigdigits=4)
    res = round(env.residue,sigdigits=4)

    printstyled(" dT=$dT res=$res\e[K\n", bold=true, color=:light_blue)
    
    # line 2:
    bar = progress_bar(T)
    progress = @sprintf("%5.3f", T*100)
    printstyled("  $(see(sw)) ", bold=true, color=:light_blue)
    print(bar)
    printstyled(" $(progress)% \e[K\n", bold=true, color=:light_blue)

    # print monitors
    for mon in model.monitors
        str     = output(mon)
        nlines += count("\n", str)
        printstyled(str, color=:light_blue)
    end

    return nlines
end
