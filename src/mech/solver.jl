# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export solve!

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
            #Ke = rand(24,24)
            #rmap=rand(1:ndofs, 24)
            #cmap=rand(1:ndofs, 24)


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
    #@show IDs

    R = reduce(vcat, Rs)
    C = reduce(vcat, Cs)
    V = reduce(vcat, Vs)

    #local K
    #try
        K = sparse(R, C, V, ndofs, ndofs)
    #catch err
        #@show Rs
        #@show R
        #error("Error using sparse function")
        #@show ndofs
        #@show err
    #end

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
function solve_step(
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

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available schemes are `:FE` and `:ME`

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps
"""
function solve!(
                dom     :: Domain,
                bcs     :: Array;
                nincs   :: Int     = 1,
                maxits  :: Int     = 5,
                autoinc :: Bool    = false,
                maxincs :: Int     = 1000000,
                tol     :: Number  = 1e-2,
                scheme  :: Symbol  = :FE,
                nouts   :: Int     = 0,
                outdir  :: String  = "",
                filekey :: String  = "out",
                verbose :: Bool    = true,
               )

    tol>0 || error("solve! : tolerance should be greater than zero")
    env = dom.env
    env.cstage += 1
    env.cinc    = 0

    if verbose
        printstyled("FEM analysis:\n", bold=true, color=:cyan)
        println("  model type: ", env.modeltype)
    end
    
    save_incs = nouts>0
    if save_incs
        if nouts>nincs
            nincs = nouts
            @info "  nincs changed to $nincs to match nouts"
        end
        if nincs%nouts != 0
            nincs = nincs - (nincs%nouts) + nouts
            @info "  nincs changed to $nincs to be a multiple of nouts"
        end

        strip(outdir) == "" && (outdir = ".")
        isdir(outdir) || error("solve!: output directory <$outdir> not fount")
        outdir[end] in ('/', '\\')  && (outdir = outdir[1:end-1])
    end


    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs)
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbose && println("  unknown dofs: $nu")
    
    # Get forces and displacements from boundary conditions
    Uex, Fex = get_bc_vals(dom, bcs)

    # Setup quantities at dofs
    if env.cstage==1
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end
    end

    # Timing
    sw = StopWatch()

    # Initial logging
    update_loggers!(dom)  # Tracking nodes, ips, elements, etc.

    # Save initial file
    if env.cstage==1 && save_incs
        update_output_data!(dom)
        save(dom, "$outdir/$filekey-0.vtk", verbose=false)
        verbose && printstyled("  $outdir/$filekey-0.vtk file written (Domain)\n", color=:green)
    end

    # Get the domain current state and backup
    State = [ ip.data for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    t  = 0.0
    dt = 1.0/nincs # initial dt value
    dt_bk = 0.0

    dT = 1.0/nouts  # output time increment for saving vtk file
    T  = dT        # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 0       # increment counter
    iout = dom.nouts     # file output counter
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
    #@show Fin
    #@show Fex

    local K::SparseMatrixCSC{Float64,Int64}

    #remountK = true

    while t < 1.0 - ttol
        # Update counters
        inc += 1
        env.cinc += 1

        #verbose && printstyled("  stage $(env.cstage) increment $inc from t=$(round(t-dt,digits=10)) to t=$(round(t,digits=10)) (dt=$(round(dt,digits=10))):", bold=true, color=:blue) # color 111
        verbose && printstyled("  stage $(env.cstage) increment $inc from t=$(round(t,digits=10)) to t=$(round(t+dt,digits=10)) (dt=$(round(dt,digits=10))):", bold=true, color=:blue) # color 111
        verbose && println()

        ΔUex, ΔFex = dt*Uex, dt*Fex     # increment of external vectors
        R   .= ΔFex    # residual
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        converged = false
        maxfails  = 3    # maximum number of it. fails with residual change less than 90%
        nfails    = 0    # counter for iteration fails
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            #if it>1; remountK=true end 
            lastres = residue # residue from last iteration

            # Try FE step
            verbose && print("    assembling... \r")
            #remountK && (K = mount_K(dom, ndofs))
            K = mount_K(dom, ndofs)
            #K = mount_K_threads(dom, ndofs)

            # Solve
            verbose && print("    solving...   \r")
            solve_step(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

            # Update
            verbose && print("    updating... \r")

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUt   = ΔUa + ΔUi
            for elem in dom.elems  
                elem_update!(elem, ΔUt, ΔFin, dt)
            end

            residue = maximum(abs, (ΔFex-ΔFin)[umap] )

            # use ME scheme
            if residue > tol && scheme == :ME
                verbose && print("    assembling... \r")
                K2 = mount_K(dom, ndofs)
                K  = 0.5*(K + K2)
                verbose && print("    solving...   \r")
                solve_step(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                copyto!.(State, StateBk)

                ΔFin .= 0.0
                ΔUt   = ΔUa + ΔUi
                for elem in dom.elems
                    elem_update!(elem, ΔUt, ΔFin, dt)
                end

                residue = maximum(abs, (ΔFex-ΔFin)[umap] )
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R = ΔFex - ΔFin  
            R[pmap] .= 0.0  # zero at prescribed positions

            if verbose
                printstyled("    it $it  ", bold=true)
                @printf(" residue: %-10.4e\n", residue)
            end

            if residue < tol;        converged = true ; remountK=false; break end
            if isnan(residue);       converged = false; break end
            if it > maxits;          converged = false; break end
            if residue > 0.9*lastres;  nfails += 1 end
            if nfails == maxfails;     converged = false; break end
        end

        if converged
            # Update forces and displacement for the current stage
            F .+= ΔFin
            U .+= ΔUa

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Update time t and dt
            t += dt

            if inc > maxincs
                printstyled("  solver maxincs = $maxincs reached (try maxincs=0)\n", color=:red)
                return false
            end

            # Check for saving output file
            if abs(t - T) < ttol
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                save(dom, "$outdir/$filekey-$iout.vtk", verbose=false)
                T += dT # find the next output time
                verbose && printstyled("  $outdir/$filekey-$iout.vtk file written (Domain)\n", color=:green)
            end

            if autoinc
                if dt_bk>0.0
                    dt = dt_bk
                else
                    dt = min(1.5*dt, 1.0/nincs)
                    dt = round(dt, digits=-ceil(Int, log10(dt))+3)  # round to 3 significant digits
                end
            end
            dt_bk = 0.0

            # Fix dt in case d+dt>T
            if t+dt>T
                dt_bk = dt
                dt = T-t
            end
        else
            # Restore counters
            env.cinc -= 1
            inc -= 1

            # Restore the state to last converged increment
            if autoinc
                verbose && println("    increment failed.")
                dt *= 0.5
                dt = round(dt, digits =-ceil(Int, log10(dt))+3)  # round to 3 significant digits
                if dt < ttol
                    printstyled("solve!: solver did not converge\n", color=:red)
                    return false
                end
            else
                printstyled("solve!: solver did not converge\n", color=:red)
                return false
            end
        end
    end

    # time spent
    verbose && println("  time spent: ", see(sw, format=:hms))

    update_output_data!(dom)

    # Update number of used increments at domain
    #dom.nincs += inc
    #dom.nouts = iout

    return true

end
