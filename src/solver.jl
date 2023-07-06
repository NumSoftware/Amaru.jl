# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct StatusLine
    model::Model
    stage::Stage
    sw::StopWatch
    msg_queue::Array{String,1}
    function StatusLine(model, stage, sw)
        return new(model, stage, sw, [])
    end
end


function message(sline::StatusLine, msg::String, color::Symbol=:default)
    iscolor = get(stdout, :color, false)
    if iscolor
        enable_ansi  = get(Base.text_colors, color, Base.text_colors[:default])
        disable_ansi = get(Base.disable_text_style, color, Base.text_colors[:default])
        msg = string(enable_ansi, msg, disable_ansi)
    end
    push!(sline.msg_queue, msg)
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


function run_status_line(sline::StatusLine)
    env    = sline.model.env
    stage  = sline.stage

    nlines = 0
    last_loop = false
    while true
        # sleep(0.1) # !fixme

        print("\r")
        nlines>0 && print("\e[$(nlines)A")

        # Print messages
        if length(sline.msg_queue)>0
            println.(sline.msg_queue, "\e[K")
            sline.msg_queue = []
        end
        
        # Print status
        T = env.T
        ΔT = env.ΔT
        progress = @sprintf("%5.3f", T*100)
        bar = progress_bar(T)
        
        # line 1:
        printstyled("  inc $(env.inc) output $(env.out)", bold=true, color=:light_blue)
        if env.transient
            t = round(env.t, sigdigits=3)
            printstyled(" t=$t", bold=true, color=:light_blue)
        end
        dT = round(ΔT,sigdigits=4)
        res = round(env.residue,sigdigits=4)

        printstyled(" dT=$dT res=$res\e[K\n", bold=true, color=:light_blue)
        
        # line 2:
        printstyled("  $(see(sline.sw)) ", bold=true, color=:light_blue)
        print(bar)
        printstyled(" $(progress)% \e[K\n", bold=true, color=:light_blue)

        # Print monitors
        nlines = 2
        for mon in sline.model.monitors
            str     = output(mon)
            nlines += count("\n", str)
            printstyled(str, color=:light_blue)
        end
        last_loop && break
        stage.status != :solving && (last_loop=true)
        yield()
    end
end


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

    # Solve linear system
    F2 = K22*U2
    U1 = zeros(nu)
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
            U1 .= NaN
            return failure("$msg\nsolve_system!: $err")
        end
    end

    maxU = 1e5 # maximum essential value
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
    quiet   = get(args, :quiet, false)
    env     = model.env
    
    cstage = findfirst(st->st.status!=:done, model.stages)
    cstage === nothing && throw(AmaruException("stage_iterator!: No stages have been set for $name"))

    solstatus = success()

    if !quiet && cstage==1 
        printstyled(name, "\n", bold=true, color=:cyan)
        println("  active threads: ", Threads.nthreads())
        # println("  model type: ", env.anaprops.stressmodel)
    end

    outdir = rstrip(outdir, ['/', '\\'])
    env.outdir = outdir

    if !isdir(outdir)
        info("solve!: creating output directory ./$outdir")
        mkpath(outdir)
    end


    if cstage==1
        logfile = open("solve.log", "w")
    else
        logfile = open("solve.log", "a")
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
            if nouts > nincs
                nincs = nouts
                quiet || info("nincs changed to $(nincs) to match nouts")
            end
            if nincs%nouts != 0 && !autoinc
                stage.nincs = nincs - (nincs%nouts) + nouts
                quiet || info("nincs changed to $nincs to be multiple of nouts")
            end
            stage.nincs = nincs
            stage.nouts = nouts
        end

        sw = StopWatch() # timing
        sline = StatusLine(model, stage, sw)
        if !quiet
            print("\e[?25l") # disable cursor
            sline_task = Threads.@spawn run_status_line(sline)
        end

        local runerror
        local error_st
        try
            solstatus = stage_solver!(model, stage, logfile, sline; args...)
            if succeeded(solstatus)
                stage.status = :done
            else
                stage.status = :failed
            end
        catch err            
            runerror = err
            flush(logfile)
            if err isa InterruptException
                stage.status = :interrupted
            else
                stage.status = :error
                error_st = stacktrace(catch_backtrace())
            end
        end
        close(logfile)

        if !quiet
            wait(sline_task)
            print("\e[?25h") # enable cursor
            solstatus.message != "" && println(solstatus.message)
        end

        if stage.status == :interrupted 
            throw(AmaruException("The analysis was interrupted"))
        elseif stage.status == :error
            # trim not important frames; search for the frame that contains current function
            idx = findfirst(contains("iterator"), string(frame) for frame in error_st)
            error_st = error_st[1:idx-1]

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

function solve!(model::FEModel; args...)
    solve!(model, model.env.anaprops; args...)
end