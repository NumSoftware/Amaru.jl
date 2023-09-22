# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ThermoAnalysis, ThermomechAnalysis

mutable struct ThermomechAnalysisProps<:Analysis
    stressmodel::String # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    T0::Float64 # reference temperature
    
    function ThermomechAnalysisProps(;stressmodel="3d", thickness=1.0, g=0.0, T0=0)
        @check stressmodel in ("plane-stress", "plane-strain", "axisymmetric", "3d")
        @check thickness>0
        @check g>=0
        @check T0>=-273.15
        return new(stressmodel, thickness, g, T0)
    end
end

ThermoAnalysis = ThermomechAnalysis = ThermomechAnalysisProps

# Assemble the global stiffness matrix
function tm_mount_global_matrices(model::Model,
                                  ndofs::Int,
                                  Δt::Float64,
                                 )
    # Assembling matrix G

    α = 1.0 # time integration factor
    T0k = model.env.ana.T0 + 273.15

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]
        RHS = zeros(ndofs)

        for elem in model.elems
            ty                      = typeof(elem)
            has_stiffness_matrix    = hasmethod(elem_stiffness, (ty,))
            has_coupling_matrix     = hasmethod(elem_coupling_matrix, (ty,))
            has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
            has_mass_matrix         = hasmethod(elem_mass_matrix, (ty,))

            # Assemble the stiffness matrix
            if has_stiffness_matrix
                K, rmap, cmap = elem_stiffness(elem)
                nr, nc = size(K)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, K[i,j])
                    end
                end
            end

            # Assemble the coupling matrices
            if has_coupling_matrix
                Cut, rmap, cmap = elem_coupling_matrix(elem)
                nr, nc = size(Cut)
                for i in 1:nr
                    for j in 1:nc
                        # matrix Cut
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, Cut[i,j])

                        # matrix Cut'
                        push!(R, cmap[j])
                        push!(C, rmap[i])
                        push!(V, T0k*Cut[i,j]) # transposed and multiplied by T0k
                    end
                end
            end

            # Assemble the conductivity matrix
            if has_conductivity_matrix
                H, rmap, cmap =  elem_conductivity_matrix(elem)
                nr, nc = size(H)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, α*Δt*H[i,j])
                    end
                end

                # Assembling RHS components
                Ut = [ dof.vals[:ut] for node in elem.nodes for dof in node.dofs if dof.name==:ut ]
                RHS[rmap] -= Δt*(H*Ut)
            end

            # Assemble the mass matrix
            if has_mass_matrix
                M, rmap, cmap =  elem_mass_matrix(elem)
                nr, nc = size(M)
                for i in 1:nr
                    for j in 1:nc
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, M[i,j])
                    end
                end
            end
        end
    end

    # generating sparse matrix G
    local G
    try
        G = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show ndofs
        @show err
    end

    yield()

    return G, RHS
end


function complete_ut_T(model::Model)
    haskey(model.node_data, "ut") || return
    Ut = model.node_data["ut"]
    T0 = model.env.ana.T0

    for elem in model.elems
        elem.shape.family==BULKCELL || continue
        elem.shape==elem.shape.basic_shape && continue
        npoints  = elem.shape.npoints
        nbpoints = elem.shape.basic_shape.npoints
        map = [ elem.nodes[i].id for i in 1:nbpoints ]
        Ute = Ut[map]
        C = elem.shape.nat_coords
        for i in nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i,:]
            N = elem.shape.basic_shape.func(R)
            Ut[id] = dot(N,Ute)
        end
    end

    model.node_data["T"] = Ut .+ T0
end


"""
    tm_solve!(D, bcs, options...) -> Bool

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

`Tmin     = 1e-9` : Pseudo-time tolerance

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available scheme is `:FE`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps

`silent = false` : If true, no information is printed
"""
function solve!(model::Model, ana::ThermomechAnalysis; args...)
    name = "Solver for thermal and thermomechanical analyses"
    status = stage_iterator!(name, tm_stage_solver!, model; args...)
    return status
end


function tm_stage_solver!(model::Model, stage::Stage; 
    tol     :: Number  = 1e-2,
    rtol    :: Number  = 1e-2,
    Tmin    :: Number  = 1e-9,
    rspan   :: Number  = 1e-2,
    scheme  :: String  = "FE",
    maxits  :: Int     = 5,
    autoinc :: Bool    = false,
    maxincs :: Int     = 1000000,
    outdir  :: String  = ".",
    outkey  :: String  = "out",
    quiet   :: Bool    = false
    )

    env = model.env
    println(env.log, "Hydromech FE analysis: Stage $(stage.id)")

    solstatus = success()
    scheme in ("FE", "ME", "BE") || error("solve! : invalid scheme \"$(scheme)\"")

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    env       = model.env
    save_outs = stage.nouts > 0
    T0        = env.ana.T0
    ftol      = tol


    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown bcs
    pmap  = nu+1:ndofs   # map for prescribed bcs
    model.ndofs = length(dofs)

    println(env.info,"unknown dofs: $nu")
    println(env.log, "unknown dofs: $nu")

    quiet || nu==ndofs && println(env.alerts, "No essential boundary conditions")

    if stage.id == 1
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name==:ut
                dof.vals[:T] = T0 # real temperature
            end
        end

        update_records!(model)

        # Save initial file and loggers
        # update_output_data!(model)
        # update_single_loggers!(model)
        # update_multiloggers!(model)
        # update_monitors!(model)
        complete_ut_T(model)
        # save_outs && save(model, "$outdir/$outkey-0.vtu", quiet=true)
    end
    lastflush = time()
    flushinterval = 5

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
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
    iout = env.out       # file output counter
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
    for elem in model.elems
        elem_internal_forces(elem, Fin)
    end
    Fex .-= Fin # add negative forces to external forces vector

    # Get global vectors from values at dofs
    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while T < 1.0-Tmin
        env.ΔT = ΔT

        # Update counters
        inc += 1
        env.inc = inc

        println(env.log, "  inc $inc")

        # Get forces and displacements from boundary conditions
        Δt = tspan*ΔT
        env.t = t + Δt
        UexN, FexN = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt

        ΔUex = UexN - U
        ΔFex = FexN - F

        ΔTcr = min(rspan, 1-T)    # time span to apply cumulated residues
        αcr  = min(ΔT/ΔTcr, 1.0)  # fraction of cumulated residues to apply
        T<1-rspan && (ΔFex .+= αcr.*Rc) # addition of residuals

        ΔUex[umap] .= 0.0
        ΔFex[pmap] .= 0.0

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        nits      = 0
        err       = 0.0
        res       = 0.0
        converged = false
        syserror  = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Predictor step for FE, ME and BE
            if scheme in ("FE", "ME", "BE")
                G, RHS = tm_mount_global_matrices(model, ndofs, Δt)
                R .+= RHS

                sysstatus = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (syserror=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(model.elems, ΔUt, Δt)
                failed(sysstatus) && (syserror=true; break)

                res = maximum(abs, (ΔFex-ΔFin)[umap] )
            end

            # Corrector step for ME and BE
            if res > ftol && scheme in ("ME", "BE")
                G2, RHS = tm_mount_global_matrices(model, ndofs, Δt)
                if scheme=="ME"
                    G = 0.5*(G + G2)
                elseif scheme=="BE"
                    G = G2
                end
                sysstatus = solve_system!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (syserror=true; break)

                copyto!.(State, StateBk)
                ΔUt    = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(model.elems, ΔUt, Δt)
                failed(sysstatus) && (syserror=true; break)

                res = maximum(abs, (ΔFex-ΔFin)[umap])
            end
            
            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            err = norm(ΔUi, Inf)/norm(ΔUa, Inf)

            @printf(env.log, "    it %d  residue: %-10.4e\n", it, res)

            res<ftol  && (converged=true; break)
            err<rtol  && (converged=true; break)

            isnan(res) && break
            it>maxits  && break

            it>1 && err>lastres && break
        end

        q = 0.0 # increment size factor for autoinc

        if syserror
            println(env.log, sysstatus.message)
            converged = false
        end
        quiet || sysstatus.message!="" && println(env.alerts, sysstatus.message, Base.default_color_warn)

        if converged
            # Update nodal natural and essential values for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Uex .= UexN
            Fex .= FexN
            Rc .= (1.0-αcr).*Rc .+ R  # update cumulated residue

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
            env.T = T
            env.residue = res

            # Check for saving output file
            # if T>Tcheck-Tmin && save_outs
            checkpoint = T>Tcheck-Tmin 
            if checkpoint


                env.out += 1
                # iout = env.out
                # rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                
                # update_output_data!(model)
                # update_multiloggers!(model)
                # save(model, "$outdir/$outkey-$iout.vtu", quiet=true)
                
                # update_embedded_disps!(active_elems, model.node_data["U"])
                Tcheck += ΔTcheck # find the next output time
            end

            flush = time()-lastflush>flushinterval || T >= 1.0-Tmin

            rstatus = update_records!(model, checkpoint=checkpoint, flush=flush)
            if failed(rstatus)
                println(env.alerts, rstatus.message)
                return rstatus
            end


            # update_single_loggers!(model, flush=flush)
            # update_monitors!(model, flush=flush)
            if flush
                # flush(env.log)
                lastflush = time()
                GC.gc()
            end

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        # q = 1+tanh(log10(ftol/res))
                        q = 1.333
                    else
                        q = 1+tanh(log10(rtol/err))
                    end
                    q = max(q, 1.1)

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
                    if T+ΔTtr>Tcheck-Tmin
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
            env.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(env.log, "      increment failed")
                if nits==1
                    # q = 1+tanh(log10(ftol/res))
                    q = 0.666
                else
                    q = 1+tanh(log10(rtol/err))
                end
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Tmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Solver did not converge. Try `autoinc=true`. ")
                break
            end
        end
    end

    failed(solstatus) && update_records!(model, solverfailed=true)
    # if !save_outs
    #     rstatus = update_records!(model, checkpoint=checkpoint, flush=false)

    #     update_output_data!(model)
    #     update_multiloggers!(model)
    # end

    return solstatus

end
