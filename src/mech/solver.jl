# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Assemble the global stiffness matrix
function mount_K(dom::Domain, 
                 ndofs::Int,
                 verbosity::Int
                )
    verbosity>1 && print("    assembling... \e[K \r")

    Threads.nthreads()>1 && return mount_K_threads(dom, ndofs, verbosity)

    R, C, V = Int64[], Int64[], Float64[]

    for elem in dom.elems
        Ke, rmap, cmap = elem_stiffness(elem)

        nr, nc = size(Ke)
        for i=1:nr
            for j=1:nc
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, Ke[i,j])
            end
        end
    end

    local K
    try
        K = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show ndofs
        @show err
    end

    return K
end


function mount_K_threads(
    dom::Domain, 
    ndofs::Int,
    verbosity::Int
)

    verbosity>1 && print("    assembling... \e[K \r")

    nelems = length(dom.elems)
    Rs = Array{Int64,1}[ [] for i=1:nelems  ]
    Cs = Array{Int64,1}[ [] for i=1:nelems  ]
    Vs = Array{Float64,1}[ [] for i=1:nelems  ]

    let Rs=Rs, Cs=Cs, Vs=Vs, dom=dom

        Threads.@threads for elem in dom.elems
            #id = elem.id
            Ke, rmap, cmap = elem_stiffness(elem)
            #IDs[elem.id] = Threads.threadid()

            nr, nc = size(Ke)
            for i=1:nr
                for j=1:nc
                    push!(Rs[elem.id], rmap[i])
                    push!(Cs[elem.id], cmap[j])
                    push!(Vs[elem.id], Ke[i,j])
                end
            end
        end

    end

    R = reduce(vcat, Rs)
    C = reduce(vcat, Cs)
    V = reduce(vcat, Vs)

    K = sparse(R, C, V, ndofs, ndofs)

    return K
end


function mount_RHS(dom::Domain, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in dom.elems
        F, map = elem_RHS(elem) # RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = F
    end
    return RHS
end


# Solves for a load/displacement increment
function solve_system!(
                       K  :: SparseMatrixCSC{Float64, Int},
                       DU :: Vect,
                       DF :: Vect,
                       nu :: Int,
                       verbosity :: Int
                      )
    verbosity>1 && print("    solving... \e[K \r")

    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    ndofs = length(DU)
    umap  = 1:nu
    pmap  = nu+1:ndofs
    if nu == ndofs
        warn("solve_system!: No essential boundary conditions")
    end 

    # Global stifness matrix
    if nu>0
        nu1 = nu+1
        K11 = K[1:nu, 1:nu]
        K12 = K[1:nu, nu1:end]
        K21 = K[nu1:end, 1:nu]
    end
    K22 = K[nu+1:end, nu+1:end]

    F1  = DF[1:nu]
    U2  = DU[nu+1:end]

    # Solve linear system
    F2 = K22*U2
    U1 = zeros(nu)
    if nu>0
        RHS = F1 - K12*U2
        
        try
            # Regularization attempt
            # S = spdiagm([ 1/maximum(abs, K11[:,i]) for i in 1:nu ])
            # LUfact = lu(K11*S)
            # U1  = S*(LUfact\RHS)

            LUfact = lu(K11)
            U1  = LUfact\RHS

            F2 += K21*U1
        catch err
            any(isnan.(K11)) && warn("solve_system!: NaN values in coefficients matrix")
            # isnan(det(K11)) && warn("solve_system!: Determinant of coefficients matrix is NaN")
            U1 .= NaN
            warn("solve_system!: $err")
            return failure("solve!: $err")
        end
    end

    maxU = 1e5 # maximum essential value
    # maximum(abs, U1)>maxU && warn("solve_system!: Possible syngular matrix")
    if maximum(abs, U1)>maxU 
        #  warn("solve_system!: Possible syngular matrix")
        return failure("solve!: Possible syngular matrix")
    end

    # Completing vectors
    DU[1:nu]     .= U1
    DF[nu+1:end] .= F2

    return success()
end


function update_state!(dom::Domain, ΔUt::Vect, ΔFin::Vect, t::Float64, verbosity::Int)
    # Update
    verbosity>1 && print("    updating... \r")

    # Get internal forces and update data at integration points (update ΔFin)
    ΔFin .= 0.0
    for elem in dom.elems
        status = elem_update!(elem, ΔUt, ΔFin, 0.0)
        failed(status) && return status
    end
    return success()
end


function update_embedded_disps!(dom::Domain)
    for elem in dom.elems.embedded
        Ue, nodemap, dimmap = elem_displacements(elem)
        U = dom.node_data["U"]
        U[nodemap, dimmap] .= Ue
    end
end


function update_status_line(dom::Domain, sw::StopWatch, inc::Int, T::Float64, ΔT::Float64, verbosity::Int, resetcursor=true)
    verbosity==0 && return
    env = dom.env

    # Print status
    progress = @sprintf("%4.2f", T*100)
    printstyled("  stage $(env.cstage) $(see(sw)) out $(env.cout) progress $(progress)% increment $inc dT=$(round(ΔT,sigdigits=4))\e[K\n", bold=true, color=:light_blue) # color 111

    # Print monitors
    nlines = 1
    for mon in dom.monitors
        str = output(mon)
        printstyled(str, color=:light_blue)
        verbosity==1 && (nlines+=count("\n", str))
    end
    resetcursor && verbosity==1 && print("\e[$(nlines)A")
end


function clean_status_line(dom::Domain, verbosity::Int)
    verbosity!=1 && return
    println("\e[K")
    nlines = 1
    for mon in dom.monitors
        str = output(mon)
        monlines = count("\n", str)
        for i in 1:monlines
            println("\e[K")
        end
        nlines+=count("\n", str)
    end
    print("\e[$(nlines)A")
end


"""
    solve!(dom, bcs, options...) :: Bool

Performs one stage static finite element analysis of a domain `dom`
subjected to a set of boundary conditions `bcs`.

# Arguments

`dom` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxincs = 1000000` : Maximum number of increments

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`Ttol     = 1e-9` : Pseudo-time tolerance

`scheme  = "FE"` : Predictor-corrector scheme at each increment. Available schemes are "FE", "ME", "BE", "Ralston"

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`printlog = false` : verbosity level from 0 (silent) to 2 (verbose)

"""
function solve!(
                dom     :: Domain,
                bcs     :: Array;
                nincs   :: Int     = 1,
                maxits  :: Int     = 5,
                autoinc :: Bool    = false,
                maxincs :: Int     = 1000000,
                tol     :: Number  = 1e-2,
                Ttol    :: Number  = 1e-9,
                rspan   :: Number  = 1e-2,
                scheme  :: Union{String,Symbol} = "FE",
                nouts   :: Int     = 0,
                outdir  :: String  = ".",
                filekey :: String  = "out",
                printlog = false,
                verbose = false,
               )
    # Arguments checking
    verbosity = 0
    printlog && (verbosity=1)
    printlog && verbose && (verbosity=2)

    scheme = string(scheme)
    scheme in ("FE", "ME", "BE", "Ralston") || error("solve! : invalid scheme \"$scheme\"")

    tol>0 || error("solve! : tolerance should be greater than zero")
    Ttol>0 || error("solve! : tolerance `Ttol `should be greater than zero")

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    sw = StopWatch() # timing

    if verbosity>0
        printstyled("Mechanical FE analysis: Stage $(env.cstage)\e[K\n", bold=true, color=:cyan)
    end

    verbosity>1 && println("  model type: ", env.modeltype)

    save_outs = nouts>0
    if save_outs #&& !autoinc
        if nouts>nincs
            nincs = nouts
            info("nincs changed to $nincs to match nouts")
        end
        if nincs%nouts != 0 && !autoinc
            nincs = nincs - (nincs%nouts) + nouts
            info("nincs changed to $nincs to be a multiple of nouts")
        end
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbosity>0 && message("unknown dofs: $nu")

    # Setup quantities at dofs
    if env.cstage==1
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end
    end

    outdir = rstrip(outdir, ['/', '\\'])
    env.outdir = outdir
    if !isdir(outdir)
        info("solve!: creating output directory ./$outdir")
        mkpath(outdir)
    end

    # Save initial file and loggers
    if env.cstage==1
        update_output_data!(dom)
        update_single_loggers!(dom)
        update_composed_loggers!(dom)
        update_monitors!(dom)
        save_outs && save(dom, "$outdir/$filekey-0.vtu", printlog=false)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT,0.01))

    ΔTbk = 0.0

    ΔTcheck = save_outs ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    inc       = 0             # increment counter
    iout      = env.cout      # file output counter
    F         = zeros(ndofs)  # total internal force for current stage
    U         = zeros(ndofs)  # total displacements for current stage
    R         = zeros(ndofs)  # vector for residuals of natural values
    ΔFin      = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa       = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi       = zeros(ndofs)  # vector of essential values for current iteration
    Rc        = zeros(ndofs)  # vector of cumulated residues
    sysstatus = ReturnStatus()
    solstatus = success()

    # Get forces and displacements from boundary conditions
    Uex, Fex = get_bc_vals(dom, bcs)

    # Get unbalanced forces
    if env.cstage==1
        Fin = zeros(ndofs)
        for elem in dom.elems
            elem_internal_forces(elem, Fin)
        end
        Fex .-= Fin # add negative forces to external forces vector
    end

    local K::SparseMatrixCSC{Float64,Int64}

    while T < 1.0-Ttol
        # Update counters
        inc += 1
        env.cinc += 1

        if inc > maxincs
            alert("solver maxincs = $maxincs reached (try maxincs=0)\n")
            return failure("$maxincs reached")
        end

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
            update_status_line(dom, sw, inc, T, ΔT, verbosity)
            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Predictor step for FE, ME and BE
            if scheme in ("FE", "ME", "BE")
                K = mount_K(dom, ndofs, verbosity)
                sysstatus = solve_system!(K, ΔUi, R, nu, verbosity)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                sysstatus = update_state!(dom, ΔUt, ΔFin, 0.0, verbosity)
                failed(sysstatus) && (errored=true; break)
                
                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            # Corrector step for ME and BE
            if residue > tol && scheme in ("ME", "BE")
                K2 = mount_K(dom, ndofs, verbosity)
                if scheme=="ME"
                    K = 0.5*(K + K2)
                elseif scheme=="BE"
                    K = K2
                end
                sysstatus = solve_system!(K, ΔUi, R, nu, verbosity)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                sysstatus = update_state!(dom, ΔUt, ΔFin, 0.0, verbosity)
                failed(sysstatus) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

            if scheme=="Ralston"
                # Predictor step
                K = mount_K(dom, ndofs, verbosity)
                ΔUit = 2/3*ΔUi
                sysstatus = solve_system!(K, ΔUit, 2/3*R, nu, verbosity)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt = ΔUa + ΔUit
                sysstatus = update_state!(dom, ΔUt, ΔFin, 0.0, verbosity)
                failed(sysstatus) && (errored=true; break)

                # Corrector step
                K2 = mount_K(dom, ndofs, verbosity)
                K = 0.25*K + 0.75*K2
                sysstatus = solve_system!(K, ΔUi, R, nu, verbosity)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (errored=true; break)

                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                sysstatus = update_state!(dom, ΔUt, ΔFin, 0.0, verbosity)
                failed(sysstatus) && (errored=true; break)

                residue = maximum(abs, (ΔFex-ΔFin)[umap])
            end

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
            verbosity>1 && notify(sysstatus.message, level=3)
            converged = false
        end

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
            env.T = T
            
            update_single_loggers!(dom)
            update_monitors!(dom)

            # Check for saving output file
            if T>Tcheck-Ttol && save_outs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                update_composed_loggers!(dom)
                update_embedded_disps!(dom)

                rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                save(dom, "$outdir/$filekey-$iout.vtu", printlog=false)
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
            env.cinc -= 1

            if autoinc
                verbosity>1 && notify("increment failed", level=3)
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
        update_output_data!(dom)
        update_composed_loggers!(dom)
    end

    update_status_line(dom, sw, inc, T, ΔT, verbosity, false)
    getlapse(sw)>60 && sound_alert()

    return solstatus
end
