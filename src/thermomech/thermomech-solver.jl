# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Assemble the global stiffness matrix
function tm_mount_global_matrices(dom::Domain,
                                  ndofs::Int,
                                  Δt::Float64,
                                  printlog=false
                                 )
    printlog && print("    assembling... \e[K \r")

    # Assembling matrix G

    R, C, V = Int64[], Int64[], Float64[]
    RHS = zeros(ndofs)

    α = 1.0 # time integration factor

    for elem in dom.elems

        ty = typeof(elem)
        has_stiffness_matrix    = hasmethod(elem_stiffness, (ty,))
        has_coupling_matrix     = hasmethod(elem_coupling_matrix, (ty,))
        has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
        has_mass_matrix         = hasmethod(elem_mass_matrix, (ty,)) # M

        # Assemble the stiffness matrix
        if has_stiffness_matrix
            K, rmap, cmap = elem_stiffness(elem)
            nr, nc = size(K)
            for i=1:nr
                for j=1:nc
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
            T0     = dom.env.T0 + 273.15
            for i=1:nr
                for j=1:nc
                    # matrix Cut
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, Cut[i,j])

                    # matrix Cut'
                    push!(R, cmap[j])
                    push!(C, rmap[i])
                    push!(V, T0*Cut[i,j]) # transposed and multiplied by T0
                end
            end
        end


        # Assemble the conductivity matrix
        if has_conductivity_matrix
            H, rmap, cmap =  elem_conductivity_matrix(elem)
            nr, nc = size(H)
            for i=1:nr
                for j=1:nc
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
            for i=1:nr
                for j=1:nc
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, M[i,j])
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

    return G, RHS
end


# Solves for a load/displacement increment
function tm_solve_system!(
                          G  :: SparseMatrixCSC{Float64, Int},
                          DU :: Vect,
                          DF :: Vect,
                          nu :: Int,
                          printlog=false
                         )
    printlog && print("    solving... \e[K \r")

    #  [  G11   G12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  G21   G22 ]  [ U2  ]    [ F2? ]

    ndofs = length(DU)
    umap  = 1:nu
    pmap  = nu+1:ndofs
    if nu == ndofs
        warn("tm_solve_system!: No essential boundary conditions")
    end

    # Global stifness matrix
    if nu>0
        nu1 = nu+1
        G11 = G[1:nu, 1:nu]
        G12 = G[1:nu, nu1:end]
        G21 = G[nu1:end, 1:nu]
    end
    G22 = G[nu+1:end, nu+1:end]

    F1  = DF[1:nu]
    U2  = DU[nu+1:end]

    # Solve linear system
    F2 = G22*U2
    U1 = zeros(nu)
    if nu>0
        RHS = F1 - G12*U2
        try
            LUfact = lu(G11)
            U1  = LUfact\RHS
            F2 += G21*U1
        catch err
            U1 .= NaN
            return failure("tm_solve!: $err")
        end
    end

    # Completing vectors
    DU[1:nu]     .= U1
    DF[nu+1:end] .= F2
    return success()
end


function complete_ut_T(dom::Domain)
    haskey(dom.node_data, "ut") || return
    Ut = dom.node_data["ut"]

    for elem in dom.elems
        elem.shape.family==SOLID_CELL || continue
        elem.shape==elem.shape.basic_shape && continue
        npoints  = elem.shape.npoints
        nbpoints = elem.shape.basic_shape.npoints
        map = [ elem.nodes[i].id for i=1:nbpoints ]
        Ute = Ut[map]
        C = elem.shape.nat_coords
        for i=nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i,:]
            N = elem.shape.basic_shape.func(R)
            Ut[id] = dot(N,Ute)
        end
    end
    dom.node_data["T"] = Ut .+ dom.env.T0
end


function tm_update_state!(dom::Domain, ΔUt::Vect, ΔFin::Vect, Δt::Float64, printlog=false)
    # Update
    printlog && print("    updating... \r")

    # Get internal forces and update data at integration points (update ΔFin)
    ΔFin .= 0.0
    for elem in dom.elems
        status = elem_update!(elem, ΔUt, ΔFin, Δt)
        failed(status) && return status
    end
    return success()
end


"""
    tm_solve!(D, bcs, options...) -> Bool

Performs one stage finite element thermo-mechanical analysis of a `domain`
subjected to a list of boundary conditions `bcs`.

# Arguments

`dom` : A finite element domain

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
function tm_solve!(
                   dom       :: Domain,
                   bcs       :: Array;
                   time_span :: Real    = NaN,
                   end_time  :: Real    = NaN,
                   nincs     :: Int     = 1,
                   maxits    :: Int     = 5,
                   autoinc   :: Bool    = false,
                   maxincs   :: Int     = 1000000,
                   tol       :: Number  = 1e-2,
                   Ttol      :: Number  = 1e-9,
                   rspan     :: Number  = 1e-3,
                   scheme    :: Union{String,Symbol} = "FE",
                   nouts     :: Int     = 0,
                   outdir    :: String  = ".",
                   filekey   :: String  = "out",
                   printlog = false,
                   verbose = false
                  )

    # Arguments checking
    verbosity = 0
    printlog && (verbosity=1)
    printlog && verbose && (verbosity=2)

    scheme = string(scheme)
    scheme in ("FE",) || error("tm_solve! : invalid scheme \"$scheme\"")

    tol>0 || error("tm_solve! : tolerance `tol `should be greater than zero")
    Ttol>0 || error("tm_solve! : tolerance `Ttol `should be greater than zero")

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    env.transient = true
    sw = StopWatch() # timing

    if !isnan(end_time)
        end_time > env.t || error("tm_solve! : end_time ($end_time) is greater that current time ($(env.t))")
        time_span = end_time - env.t
    end
    isnan(time_span) && error("tm_solve!: neither time_span nor end_time were set.")

    if verbosity>0
        printstyled("Thermomechanical FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
        println("  from t=$(round(dom.env.t,digits=4)) to t=$(round(dom.env.t+time_span,digits=3))")
    end

    verbosity>1 && println("  model type: ", env.modeltype)

    save_outs = nouts>0
    if save_outs && !autoinc
        if nouts>nincs
            nincs = nouts
            info("nincs changed to $nincs to match nouts")
        end
        if nincs%nouts != 0
            nincs = nincs - (nincs%nouts) + nouts
            info("nincs changed to $nincs to be a multiple of nouts")
        end
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements and pw
    pmap  = nu+1:ndofs   # map for prescribed displacements and pw
    dom.ndofs = length(dofs)
    verbosity>0 && message("unknown dofs: $nu")

    # Setup quantities at dofs
    if env.cstage==1
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name==:ut
                dof.vals[:T] = env.T0 # real temperature
            end
        end
    end

    outdir = rstrip(outdir, ['/', '\\'])
    env.outdir = outdir
    if !isdir(outdir)
        info("tm_solve!: creating output directory ./$outdir")
        mkpath(outdir)
    end

    # Save initial file and loggers
    if env.cstage==1
        update_output_data!(dom)
        complete_ut_T(dom)
        update_single_loggers!(dom)
        update_composed_loggers!(dom)
        save_outs && save(dom, "$outdir/$filekey-0.vtu", printlog=printlog)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t    = env.t     # current time
    tend = t + time_span # end time

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT,0.01))
    ΔTbk = 0.0

    ΔTcheck = save_outs ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    inc  = 0             # increment counter
    iout = env.cout      # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    status = ReturnStatus()

    Fex  = zeros(ndofs)  # vector of external loads
    Uex  = zeros(ndofs)  # vector of external essential values

    Uex, Fex = get_bc_vals(dom, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # Get unbalanced forces
    Fin = zeros(ndofs)
    for elem in dom.elems
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

    while T < 1.0-Ttol
        # Update counters
        inc += 1
        env.cinc += 1
        env.T = T

        if inc > maxincs
            alert("solver maxincs = $maxincs reached (try maxincs=0)\n")
            return false
        end

        progress = @sprintf("%4.2f", T*100)
        verbosity>0 && printstyled("  stage $(env.cstage) $(see(sw)) progress $(progress)% increment $inc dT=$(round(ΔT,sigdigits=4))\e[K\r", bold=true, color=:blue) # color 111
        verbosity>1 && println()

        # Get forces and displacements from boundary conditions
        Δt = time_span*ΔT
        env.t = t + Δt
        UexN, FexN = get_bc_vals(dom, bcs, t+Δt) # get values at time t+Δt

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
        maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        nfails    = 0  # counter for iteration fails
        nits      = 0
        residue1  = 0.0
        converged = false
        errored   = false
        for it=1:maxits
            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            G, RHS = tm_mount_global_matrices(dom, ndofs, Δt, printlog)

            R .+= RHS

            # Solve
            status = tm_solve_system!(G, ΔUi, R, nu, printlog)   # Changes unknown positions in ΔUi and R
            failed(status) && (errored=true; break)

            copyto!.(State, StateBk)
            ΔUt = ΔUa + ΔUi
            status = tm_update_state!(dom, ΔUt, ΔFin, Δt, printlog)
            failed(status) && (errored=true; break)

            residue = maximum(abs, (ΔFex-ΔFin)[umap] )

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            if verbosity>1
                printstyled("    it $it  ", bold=true)
                @printf(" residue: %-10.4e\n", residue)
            end

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
            verbosity>1 && notify(status.message, level=3)
            converged = false
        end

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
                    dof.vals[:T] = U[i] + env.T0
                end
            end

            update_single_loggers!(dom)

            # Update time
            t += Δt
            T += ΔT

            # Check for saving output file
            if T>Tcheck-Ttol && save_outs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                complete_ut_T(dom)
                update_composed_loggers!(dom)
                save(dom, "$outdir/$filekey-$iout.vtu", printlog=printlog)
                Tcheck += ΔTcheck # find the next output time
            end

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
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            # Restore counters
            inc -= 1
            env.cinc -= 1

            if autoinc
                verbosity>1 && notify("increment failed", level=3)
                q = (1+tanh(log10(tol/residue1)))
                q = clamp(q, 0.2, 0.9)
                errored && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Ttol
                    alert("tm_solve!: solver did not converge")
                    return false
                end
            else
                alert("tm_solve!: solver did not converge")
                return false
            end
        end

    end

    if !save_outs
        update_output_data!(dom)
        complete_ut_T(dom)
        update_composed_loggers!(dom)
    end

    # time spent
    progress = @sprintf("%4.2f", T*100)
    verbosity>1 && printstyled("  stage $(env.cstage) $(see(sw)) progress $(progress)%\e[K\n", bold=true, color=:blue) # color 111
    if verbosity>0 
        message("valid increments: ", inc)
        message("time spent: ", see(sw, format=:hms))
    end
    getlapse(sw)>60 && sound_alert()

    return success()
end
