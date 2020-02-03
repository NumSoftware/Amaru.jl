# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

#export solve!
# Assemble the global stiffness matrix
function tm_mount_G_RHS(dom::Domain, ndofs::Int, Δt::Float64)

    # Assembling matrix G

    R, C, V = Int64[], Int64[], Float64[]
    RHS = zeros(ndofs)

    α = 1.0 # time integration factor

    for elem in dom.elems









        # Assemble the stiffness matrix
        K, rmap, cmap = elem_stiffness_matrix(elem)
        nr, nc = size(K)
        for i=1:nr
            for j=1:nc
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, K[i,j])
            end
        end

        # Assemble the coupling matrices
        Cup, rmap, cmap = elem_coupling_matrix(elem)
        nr, nc = size(Cup)
        for i=1:nr
            for j=1:nc
                # matrix Cup
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, Cup[i,j])

                # matrix Cup'
                push!(R, cmap[j])
                push!(C, rmap[i])
                push!(V, Cup[i,j])
            end
        end


        # Assemble the conductivity matrix
        H, rmap, cmap = elem_conductivity_matrix(elem)
        nr, nc = size(H)
        for i=1:nr
            for j=1:nc
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, α*Δt*H[i,j])
            end
        end

        # Assemble the mass matrix
        M, rmap, cmap = elem_mass_matrix(elem)
        nr, nc = size(M)
        for i=1:nr
            for j=1:nc
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, M[i,j])
            end
        end

        # Assembling RHS components
        Ut = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes ]
        RHS[rmap] -= Δt*(H*Ut)

        # Assemble ramaining RHS vectors
        #Q, map = elem_RHS_vector(elem)
        #RHS[map] += Δt*Q
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
function tm_solve_step!(G::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, nu::Int)
    #  [  G11   G12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  G21   G22 ]  [ U2  ]    [ F2? ]

    ndofs = length(DU)
    umap  = 1:nu
    pmap  = nu+1:ndofs
    if nu == ndofs
        @warn "solve!: No essential boundary conditions."
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
            @warn "solve!: $err"
            U1 .= NaN
        end
    end

    # Completing vectors
    DU[1:nu]     .= U1
    DF[nu+1:end] .= F2
end

"""
    tm_solve!(D, bcs, options...) -> Bool

Performs one stage finite element analysis of a domain `D`
subjected to an array of boundary conditions `bcs`.

# Arguments

`dom` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)


# Keyword arguments

`time_span` : The simulated time

`end_time` : The end time of the simulation

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxincs = 1000000` : Maximum number of increments

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available schemes are `:FE` and `:ME`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps
"""
function tm_solve!(
                   dom       :: Domain,
                   bcs       :: Array;
                   time_span :: Float64 = NaN,
                   end_time  :: Float64 = NaN,
                   nincs     :: Int     = 1,
                   maxits    :: Int     = 5,
                   autoinc   :: Bool    = false,
                   maxincs   :: Int     = 1000000,
                   tol       :: Number  = 1e-2,
                   scheme    :: Symbol  = :FE,
                   nouts     :: Int     = 0,
                   outdir    :: String  = "",
                   filekey   :: String  = "out",
                   verbose   :: Bool    = false,
                   silent    :: Bool    = false
                  )
    env = dom.env
    env.cstage += 1
    env.cinc    = 0

    if !silent
        printstyled("Thermomechanical FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
        sw = StopWatch() # timing
    end


    # Arguments checking
    silent && (verbose=false)

    tol>0 || error("solve! : tolerance should be greater than zero")

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


    if !isnan(end_time)
        time_span = end_time - dom.env.t
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements and pw
    pmap  = nu+1:ndofs   # map for prescribed displacements and pw
    dom.ndofs = length(dofs)
    silent || println("  unknown dofs: $nu")


    # Get array with all integration points
    ips = [ ip for elem in dom.elems for ip in elem.ips ]
    # Get the domain current state and backup
    State = [ ip.data for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Save initial file and loggers
    if env.cstage==1
        # Setup initial quantities at dofs
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end

        update_loggers!(dom)  # Tracking nodes, ips, elements, etc.
        update_output_data!(dom) # Updates data arrays in domain

        if save_incs
            save(dom, "$outdir/$filekey-0.vtk", verbose=false)
            verbose && printstyled("  $outdir/$filekey-0.vtk file written (Domain)\n", color=:green)
        end
    end


    # Incremental analysis
    t    = dom.env.t # current time
    tini = t # initial time
    tend = t + time_span  # end time
    Δt = time_span/nincs # initial Δt value

    dT = time_span/nouts  # output time increment for saving vtk file
    T  = t + dT        # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 1       # increment counter
    iout = env.cout     # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration

    Fex  = zeros(ndofs)  # vector of external loads
    Uex  = zeros(ndofs)  # vector of external essential values

    Uex, Fex = get_bc_vals(dom, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    while t < tend - ttol

        #verbose && printstyled("  increment $inc from t=$(round(t,sigdigits=9)) to t=$(round(t+Δt,sigdigits=9)) (dt=$(round(Δt,sigdigits=9))):\n", bold=true, color=:blue) # color 111
        progress = round((t-tini)/time_span*100, digits=3)
        silent || printstyled("  stage $(env.cstage) progress $progress% increment $inc dt=$(round(Δt,sigdigits=2))\033[K\r", bold=true, color=:blue) # color 111

        # Get forces and displacements from boundary conditions
        dom.env.t = t + Δt
        UexN, FexN = get_bc_vals(dom, bcs, t+Δt) # get values at time t+Δt

        ΔUex = UexN - U
        ΔFex = FexN - F

        ΔUex[umap] .= 0.0
        ΔFex[pmap] .= 0.0

        R   .= ΔFex    # residual
        ΔUa .= 0.0
        ΔUi .= ΔUex    # essential values at iteration i

        # Newton Rapshon iterations
        residue   = 0.0
        converged = false
        maxfails  = 3    # maximum number of it. fails with residual change less than 90%
        nfails    = 0    # counter for iteration fails
        local G::SparseMatrixCSC{Float64,Int64}
        local RHS::Array{Float64,1}
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            verbose && print("    assembling... \r")
            G, RHS = tm_mount_G_RHS(dom, ndofs, it==1 ? Δt : 0.0 ) # TODO: check for Δt after iter 1

            R .+= RHS

            # Solve
            verbose && print("    solving...   \r")
            tm_solve_step!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

            # Update
            verbose && print("    updating... \r")

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUt   = ΔUa + ΔUi
            for elem in dom.elems
                elem_update!(elem, ΔUt, ΔFin, Δt)
            end

            residue = maximum(abs, (ΔFex-ΔFin)[umap] )

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R = ΔFex - ΔFin  #
            R[pmap] .= 0.0  # Zero at prescribed positions

            if verbose
                printstyled("    it $it  ", bold=true)
                @printf(" residue: %-10.4e\n", residue)
            else
                if !silent
                    printstyled("  increment $inc: ", bold=true, color=:blue)
                    printstyled("  it $it  ", bold=true)
                    @printf("residue: %-10.4e  \r", residue)
                end
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
            Uex .= UexN
            Fex .= FexN

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Check for saving output file
            Tn = t + Δt
            if Tn+ttol>=T && saveincs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                save(dom, "$outdir/$filekey-$iout.vtk", verbose=false)
                T = Tn - mod(Tn, dT) + dT
                silent || printstyled("  $outdir/$filekey-$iout.vtk file written (Domain) \033[K \n",color=:green)
            end

            # Update time t and Δt
            inc += 1
            t   += Δt

            # Get new Δt
            if autoinc
                Δt = min(1.5*Δt, 1.0/nincs)
                Δt = round(Δt, sigdigits=3)
                Δt = min(Δt, 1.0-t)
            end
        else
            # Restore counters
            env.cinc -= 1
            inc -= 1

            if autoinc
                verbose && println("    increment failed.")
                Δt = round(0.5*Δt, sigdigits=3)
                if Δt < ttol
                    printstyled("solve!: solver did not converge \033[K \n", color=:red)
                    return false
                end
            else
                printstyled("solve!: solver did not converge \033[K \n", color=:red)
                return false
            end
        end
    end

    # time spent
    silent || println("  time spent: ", see(sw, format=:hms), "\033[K")

    update_output_data!(dom)

    return true

end
