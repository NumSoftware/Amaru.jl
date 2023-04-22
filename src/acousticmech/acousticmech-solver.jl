# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


function am_mount_M(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_mass, (ty,)) || continue

            Me, rmap, cmap = elem_acoustic_mass(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end

function am_mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_stiffness, (ty,)) || continue

            Ke, rmap, cmap = elem_acoustic_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


# function complete_ut_T(model::Model)
#     haskey(model.node_data, "ut") || return
#     Ut = model.node_data["ut"]
#     T0 = get(model.env.params, :T0, 0.0)

#     for elem in model.elems
#         elem.shape.family==BULKCELL || continue
#         elem.shape==elem.shape.basic_shape && continue
#         npoints  = elem.shape.npoints
#         nbpoints = elem.shape.basic_shape.npoints
#         map = [ elem.nodes[i].id for i in 1:nbpoints ]
#         Ute = Ut[map]
#         C = elem.shape.nat_coords
#         for i in nbpoints+1:npoints
#             id = elem.nodes[i].id
#             R = C[i,:]
#             N = elem.shape.basic_shape.func(R)
#             Ut[id] = dot(N,Ute)
#         end
#     end

#     model.node_data["T"] = Ut .+ T0
# end


"""
    am_solve!(D, bcs, options...) -> Bool

Performs one stage finite element thermo-mechanical analysis of a `domain`
subjected to a list of boundary conditions `bcs`.

# Arguments

`model` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`time_span = NaN` : Time lapse for the transient analysis in the current stage

`end_time = NaN` : Final time for the transient analysis in the current stage

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxincs = 1000000` : Maximum number of increments

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`Ttol     = 1e-9` : Pseudo-time tolerance

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available scheme is `:FE`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps

`silent = false` : If true, no information is printed
"""
function am_solve!(model::Model; args...)
    name = "Acoustic solver"
    st = stage_iterator!(name, am_stage_solver!, model; args...)
    return st
end


function am_stage_solver!(model::Model, stage::Stage, logfile::IOStream, sline::StatusLine; 
    tol     :: Number  = 1e-2,
    Ttol    :: Number  = 1e-9,
    rspan   :: Number  = 1e-2,
    maxits  :: Int     = 5,
    autoinc :: Bool    = false,
    maxincs :: Int     = 1000000,
    outdir  :: String  = ".",
    outkey  :: String  = "out",
    quiet  :: Bool    = false
                  )

    println(logfile, "Acoustic FE analysis: Stage $(stage.id)")
    stage.status = :solving

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    env       = model.env
    save_outs = stage.nouts > 0

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs    = length(dofs)
    umap     = 1:nu         # map for unknown bcs
    pmap     = nu+1:ndofs   # map for prescribed bcs
    model.ndofs = length(dofs)
    println(logfile, "unknown dofs: $nu")
    message(sline, "  unknown dofs: $nu")

    quiet || nu==ndofs && message(sline, "solve_system!: No essential boundary conditions", Base.warn_color)

    if stage.id == 1
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end

        # Save initial file and loggers
        update_output_data!(model)
        update_single_loggers!(model)
        update_multiloggers!(model)
        update_monitors!(model)
        # complete_ut_T(model)
        save_outs && save(model, "$outdir/$outkey-0.vtu", quiet=true)
    end


    # Get the domain current state and backup
    State = [ ip.state for elem in model.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t    = env.t     # current time

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT,0.01))

    ΔTbk = 0.0
    ΔTcheck = save_outs ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    inc  = 0             # increment counter
    iout = env.stagebits.out # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    sysstatus = ReturnStatus()

    # Get boundary conditions
    Uex, Fex = get_bc_vals(model, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # Get unbalanced forces
    Fin = zeros(ndofs)
    # for elem in model.elems
    #     elem_internal_forces(elem, Fin)
    # end
    Fex .-= Fin # add negative forces to external forces vector

    # Get global vectors from values at dofs
    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while T < 1.0-Ttol
        env.stagebits.ΔT = ΔT

        # Update counters
        inc += 1
        env.stagebits.inc = inc

        println(logfile, "  inc $inc")

        if inc > maxincs
            quiet || message(sline, "solver maxincs = $maxincs reached (try maxincs=0)", Base.default_color_error)
            return failure("$maxincs reached")
        end

        # Get forces and displacements from boundary conditions
        Δt = tspan*ΔT
        env.t = t + Δt
        UexN, FexN = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt

        ΔUex = UexN - U
        ΔFex = FexN - F

        # ΔTcr = min(rspan, 1-T)    # time span to apply cumulated residues
        # αcr  = min(ΔT/ΔTcr, 1.0)  # fraction of cumulated residues to apply
        # T<1-rspan && (ΔFex .+= αcr.*Rc) # addition of residuals

        ΔUex[umap] .= 0.0
        ΔFex[pmap] .= 0.0

        # R   .= ΔFex
        Fex_Fin = Fex-Fina    # residual
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        nfails    = 0  # counter for iteration fails
        nits      = 0
        residue1  = 0.0
        converged = false
        errored   = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            K = am_mount_K(model.elems, ndofs)
            M = am_mount_M(model.elems, ndofs)

            KK = K + (4/(Δt^2))*M # pseudo-stiffness matrix
            ΔQQ = Fex_Fin + M*(A + 4*V/Δt - 4*ΔUa/Δt^2) 

            # Solve
            solve_system!(KK, ΔUi, ΔQQ, nu)

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUt   = ΔUa + ΔUi
            ΔFin, sysstatus = update_state!(model.elems, ΔUt, Δt)

            Fina = Fin + ΔFin            

            TFin = Fina + C*Va + M*Aa  # Internal force including dynamic effects

            residue = maximum(abs, (Fex-TFin)[umap] )
            
            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            Fex_Fin .= Fex .- Fina  # Check this variable, it is not the residue actually
            Fex_Fin[pmap] .= 0.0  # Zero at prescribed positions

            @printf(logfile, "    it %d  residue: %-10.4e\n", it, residue)

            it==1 && (residue1=residue)
            residue > tol && (Fina -= ΔFin)
            residue < tol  && (converged=true; break)
            isnan(residue) && break
            it>maxits      && break
            it>1 && residue>lastres && break
            residue>0.9*lastres && (nfails+=1)
            nfails==maxfails    && break
        end

        q = 0.0 # increment size factor for autoinc

        if errored
            println(logfile, sysstatus.message)
            converged = false
        end
        quiet || sysstatus.message!="" && message(sline, sysstatus.message, Base.default_color_warn)

        if converged
            # Update nodal natural and essential values for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Uex .= UexN
            Fex .= FexN

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
                if dof.name==:ut
                    dof.vals[:T] = U[i] + T0
                end
            end

            # Update time
            t += Δt
            T += ΔT
            env.t = t
            env.stagebits.T = T
            env.stagebits.residue = residue

            # Check for saving output file
            if T>Tcheck-Ttol && save_outs
                env.stagebits.out += 1
                iout = env.stagebits.out
                rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                
                update_output_data!(model)
                # update_embedded_disps!(active_elems, model.node_data["U"])

                update_multiloggers!(model)
                save(model, "$outdir/$outkey-$iout.vtu", quiet=true)

                Tcheck += ΔTcheck # find the next output time
            end

            update_single_loggers!(model)
            update_monitors!(model)
            flush(logfile)

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        q = (1+tanh(log10(tol/residue1)))^1
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
                    if T+ΔTtr>Tcheck-Ttol
                        ΔTbk = ΔT
                        ΔT = Tcheck-T
                        @assert ΔT>=0.0
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            # Restore counters
            inc -= 1
            env.stagebits.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(logfile, "      increment failed")
                q = (1+tanh(log10(tol/residue1)))
                q = clamp(q, 0.2, 0.9)
                errored && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Ttol
                    solstatus = failure("solver did not converge")
                    break
                end
            else
                solstatus = failure("solver did not converge")
                break
            end
        end
    end

    if !save_outs
        update_output_data!(model)
        complete_ut_T(model)
        update_multiloggers!(model)
    end

    return solstatus

end
