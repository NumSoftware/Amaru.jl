# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Assemble the global stiffness matrix
function mount_K(dom::Domain, ndofs::Int)
    Threads.nthreads()>1 && return mount_K_threads(dom, ndofs)

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

function mount_K_threads(dom::Domain, ndofs::Int)

    nelems = length(dom.elems)
    Rs = Array{Int64,1}[ [] for i=1:nelems  ]
    Cs = Array{Int64,1}[ [] for i=1:nelems  ]
    Vs = Array{Float64,1}[ [] for i=1:nelems  ]
    #IDs = zeros(Int64, nelems)

    #GC.gc()

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
function solve_step!(
                    K  :: SparseMatrixCSC{Float64, Int},
                    DU :: Vect,
                    DF :: Vect,
                    nu :: Int
                   )
    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    ndofs = length(DU)
    umap  = 1:nu
    pmap  = nu+1:ndofs
    if nu == ndofs
        @warn "solve!: No essential boundary conditions."
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
            LUfact = lu(K11)
            U1  = LUfact\RHS
            F2 += K21*U1
        catch err
            @warn "solve!: $err"
            U1 .= NaN
        end
    end

    # Completing vectors
    DU[1:nu]     .= U1
    DF[nu+1:end] .= F2
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

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available schemes are `:FE` and `:ME`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps

`silent = false` : If true, no information is printed
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
                scheme  :: Symbol  = :FE,
                nouts   :: Int     = 0,
                outdir  :: String  = "",
                filekey :: String  = "out",
                verbose :: Bool    = false,
                silent  :: Bool    = false,
               )

    # Arguments checking
    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    tol>0 || error("solve! : tolerance should be greater than zero")
    Ttol>0 || error("solve! : tolerance `Ttol `should be greater than zero")

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    sw = StopWatch() # timing

    if verbosity>0
        printstyled("Mechanical FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
    end

    verbosity>1 && println("  model type: ", env.modeltype)

    save_outs = nouts>0
    if save_outs
        if nouts>nincs
            nincs = nouts
            @info "  nincs changed to $nincs to match nouts"
        end
        if nincs%nouts != 0
            nincs = nincs - (nincs%nouts) + nouts
            @info "  nincs changed to $nincs to be a multiple of nouts"
        end
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbosity>0 && println("  unknown dofs: $nu")

    # Get forces and displacements from boundary conditions
    Uex, Fex = get_bc_vals(dom, bcs)

    # Setup quantities at dofs
    if env.cstage==1
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end
    end

    # Save initial file and loggers
    if env.cstage==1
        outdir = strip(outdir, ['/', '\\'])
        strip(outdir) == "" && (outdir = ".")
        isdir(outdir) || error("solve!: output directory <$outdir> not found")

        update_output_data!(dom)
        update_single_loggers!(dom)
        update_composed_loggers!(dom)
        save_outs && save(dom, "$outdir/$filekey-0.vtu", silent=silent)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    ΔT_bk = 0.0

    ΔTout = 1.0/nouts    # output time increment for saving output file
    Tout  = ΔTout        # output time for saving the next output file

    inc  = 0             # increment counter
    iout = env.cout      # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration

    # Get unbalanced forces
    Fin = zeros(ndofs)
    for elem in dom.elems
        elem_internal_forces(elem, Fin)
    end
    Fex .-= Fin # add negative forces to external forces vector

    local K::SparseMatrixCSC{Float64,Int64}

    while T < 1.0 - Ttol
        # Update counters
        inc += 1
        env.cinc += 1

        if inc > maxincs
            printstyled("  solver maxincs = $maxincs reached (try maxincs=0)\n", color=:red)
            return false
        end

        progress = @sprintf("%4.2f", T*100)
        verbosity>0 && printstyled("  stage $(env.cstage) $(see(sw)) progress $(progress)% increment $inc dT=$(round(ΔT,sigdigits=4))\033[K\r", bold=true, color=:blue) # color 111
        verbosity>1 && println()

        ΔUex, ΔFex = ΔT*Uex, ΔT*Fex     # increment of external vectors
        R   .= ΔFex    # residual
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        converged = false
        maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        nfails    = 0  # counter for iteration fails
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            verbosity>1 && print("    assembling... \r")
            K = mount_K(dom, ndofs)

            # Solve
            verbosity>1 && print("    solving...   \r")
            solve_step!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

            # Update
            verbosity>1 && print("    updating... \r")

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUt   = ΔUa + ΔUi
            for elem in dom.elems
                elem_update!(elem, ΔUt, ΔFin, 0.0)
            end

            residue = maximum(abs, (ΔFex-ΔFin)[umap] )

            # use ME scheme
            if residue > tol && scheme == :ME
                verbose && print("    assembling... \r")
                K2 = mount_K(dom, ndofs)
                K  = 0.5*(K + K2)
                verbose && print("    solving...   \r")
                solve_step!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                copyto!.(State, StateBk)

                ΔFin .= 0.0
                ΔUt   = ΔUa + ΔUi
                for elem in dom.elems
                    elem_update!(elem, ΔUt, ΔFin, 0.0)
                end

                residue = maximum(abs, (ΔFex-ΔFin)[umap] )
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R = ΔFex - ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions

            if verbosity>1
                printstyled("    it $it  ", bold=true)
                @printf(" residue: %-10.4e\n", residue)
            end

            if residue < tol;        converged = true ; break end
            if isnan(residue);       converged = false; break end
            if it > maxits;          converged = false; break end
            if residue > 0.9*lastres;  nfails += 1 end
            if nfails == maxfails;     converged = false; break end
        end

        if converged
            # Update forces and displacement for the current stage
            U .+= ΔUa
            F .+= ΔFin

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            update_single_loggers!(dom)

            # Update time
            T += ΔT

            # Check for saving output file
            if abs(T - Tout) < Ttol && save_outs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                update_composed_loggers!(dom)
                save(dom, "$outdir/$filekey-$iout.vtu", silent=silent)
                Tout += ΔTout # find the next output time
            end

            if autoinc
                if ΔT_bk>0.0
                    ΔT = ΔT_bk
                else
                    ΔT = min(1.5*ΔT, 1.0/nincs)
                end
            end
            ΔT_bk = 0.0

            # Fix ΔT in case T+ΔT>Tout
            if T+ΔT>Tout
                ΔT_bk = ΔT
                ΔT = Tout-T
            end
        else
            # Restore counters
            inc -= 1
            env.cinc -= 1
            ΔT_bk = ΔT

            # Restore the state to last converged increment
            if autoinc
                verbosity>1 && println("    increment failed.")
                ΔT *= 0.5
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Ttol
                    printstyled("solve!: solver did not converge \033[K \n", color=:red)
                    return false
                end
            else
                printstyled("solve!: solver did not converge \033[K \n", color=:red)
                return false
            end
        end
    end

    if !save_outs
        update_output_data!(dom)
        update_composed_loggers!(dom)
    end

    # time spent
    verbosity>1 && printstyled("  stage $(env.cstage) $(see(sw)) progress 100%\033[K\n", bold=true, color=:blue) # color 111
    verbosity==1 && println("  time spent: ", see(sw, format=:hms), "\033[K")
    getlapse(sw)>60 && sound_alert()

    return true
end
