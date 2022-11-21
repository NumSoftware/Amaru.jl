# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Assemble the global stiffness matrix
function mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            Ke, rmap, cmap = elem_stiffness(elem)

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

    local K
    try
        K = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show err
    end

    yield()

    return K
end


function mount_RHS(model::Model, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in active_elems
        F, map = elem_RHS(elem) # RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = F
    end
    return RHS
end


# Get internal forces and update data at integration points
function update_state!(active_elems::Array{<:Element,1}, ΔUt::Vect, t::Float64)

    ndofs = length(ΔUt)
    @withthreads begin 
        ΔFin     = zeros(ndofs)
        statuses = ReturnStatus[]
        for elem in active_elems
            ΔF, map, status = elem_update!(elem, ΔUt, t)
            if failed(status) 
                push!(statuses, status)
                break
            end
            ΔFin[map] .+= ΔF
        end
    end

    yield()
    
    length(statuses)>0 && return ΔFin, statuses[1]
    any(isnan.(ΔFin)) && return failure("solve_system!: NaN values in internal forces vector")
    return ΔFin, success()
end


function update_embedded_disps!(active_elems::Array{<:Element,1}, U::Matx)
    for elem in active_elems.embeddeds
        Ue, nodemap, dimmap = elem_displacements(elem)
        U[nodemap, dimmap] .= Ue
    end
end


"""
    solve!(model,options...) :: Bool

Performs one stage static finite element analysis of a domain `model`
subjected to a set of boundary conditions `bcs`.

# Arguments

`model` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`Ttol     = 1e-9` : Pseudo-time tolerance

`scheme  = "FE"` : Predictor-corrector scheme at each increment. Available schemes are "FE", "ME", "BE", "Ralston"

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`outkey = ""` : File key for output files

`quiet = false` : verbosity level from 0 (silent) to 2 (verbose)

"""
function mech_solve!(model::Model; args...)
    name = "Mechanical solver"
    st = stage_iterator!(name, mech_stage_solver!, model; args...)
    return st
end
solve! = mech_solve!


function mech_stage_solver!(model::Model, stage::Stage, logfile::IOStream, sline::StatusLine; 
    tol     :: Number  = 1e-2,
    Ttol    :: Number  = 1e-8,
    rspan   :: Number  = 1e-2,
    scheme  :: String  = "FE",
    maxits  :: Int     = 5,
    autoinc :: Bool    = false,
    outdir  :: String  = ".",
    outkey  :: String  = "out",
    quiet   :: Bool    = false
    )

    println(logfile, "Mechanical FE analysis: Stage $(stage.id)")
    stage.status = :solving

    solstatus = success()
    scheme in ("FE", "ME", "BE", "Ralston") || error("solve! : invalid scheme \"$(scheme)\"")

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    env       = model.env
    save_outs = stage.nouts > 0

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
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
        update_composed_loggers!(model)
        update_monitors!(model)
        save_outs && save(model, "$outdir/$outkey-0.vtu", quiet=true)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
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

    # Get forces and displacements from boundary conditions
    Uex, Fex = get_bc_vals(model, bcs)

    # Get unbalanced forces from activated elements
    if length(stage.toactivate)>0
        Fin = zeros(ndofs)
        for elem in stage.toactivate
            elem_internal_forces(elem, Fin)
        end
        Fex .-= Fin # add negative forces to external forces vector
    end

    local K::SparseMatrixCSC{Float64,Int64}

    while T < 1.0-Ttol
        # sleep(0.0)
        env.stagebits.ΔT = ΔT

        # Update counters
        inc += 1
        env.stagebits.inc = inc

        println(logfile, "  inc $inc")

        ΔUex, ΔFex = ΔT*Uex, ΔT*Fex     # increment of external vectors

        ΔTcr = min(rspan, 1-T)    # time span to apply cumulated residues
        αcr  = min(ΔT/ΔTcr, 1.0)  # fraction of cumulated residues to apply
        T<1-rspan && (ΔFex .+= αcr.*Rc) # addition of residuals

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex  # essential values at iteration i

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

            # Predictor step for FE, ME and BE
            if scheme in ("FE", "ME", "BE")
                K = mount_K(active_elems, ndofs)
                sysstatus = solve_system!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (errored=true; break)
                
                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            # Corrector step for ME and BE
            if residue > tol && scheme in ("ME", "BE")
                K2 = mount_K(active_elems, ndofs)
                if scheme=="ME"
                    K = 0.5*(K + K2)
                elseif scheme=="BE"
                    K = K2
                end
                sysstatus = solve_system!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            if scheme=="Ralston"
                # Predictor step
                K = mount_K(active_elems, ndofs)
                ΔUit = 2/3*ΔUi
                sysstatus = solve_system!(K, ΔUit, 2/3*R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt = ΔUa + ΔUit
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (errored=true; break)

                # Corrector step
                K2 = mount_K(active_elems, ndofs)
                K = 0.25*K + 0.75*K2
                sysstatus = solve_system!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            @printf(logfile, "    it %d  residue: %-10.4e\n", it, residue)

            it==1 && (residue1=residue)
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
            # Update forces and displacement for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Rc .= (1.0-αcr).*Rc .+ R  # update cumulated residue

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end
            
            # Update time
            T += ΔT
            env.stagebits.T = T
            env.stagebits.residue = residue

            # Check for saving output file
            if T>Tcheck-Ttol && save_outs
                env.stagebits.out += 1
                iout = env.stagebits.out
                rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                
                update_output_data!(model)
                update_embedded_disps!(active_elems, model.node_data["U"])

                update_composed_loggers!(model)
                save(model, "$outdir/$outkey-$iout.vtu", quiet=true) #!

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
        update_composed_loggers!(model)
    end

    return solstatus

end