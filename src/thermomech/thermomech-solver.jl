# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Assemble the global stiffness matrix
function tm_mount_global_matrices(dom::Domain, ndofs::Int, Δt::Float64)

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


function complete_ut_T(dom::Domain)
    haskey(dom.point_data, "ut") || return
    Ut = dom.point_data["ut"]

    for elem in dom.elems
        elem.shape.family==SOLID_SHAPE || continue
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

    dom.point_data["T"] = Ut .+ dom.env.T0
end


"""
    tm_solve!(D, bcs, options...) -> Bool

Performs one stage finite element analysis of a domain `D`
subjected to an array of boundary conditions `bcs`.

Available options are:

`verbose=true` : If true, provides information of the analysis steps

`tol=1e-2` : Tolerance for the absolute error in forces

`nincs=1` : Number of increments

`autoinc=false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxits=5` : The maximum number of Newton-Rapson iterations per increment

`nouts=0` : Number of output files per analysis

`scheme= :FE` : Predictor-corrector scheme at iterations. Available schemes are `:FE` and `:ME`

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
                   scheme    :: Symbol  = :FE,
                   nouts     :: Int     = 0,
                   outdir    :: String  = "",
                   filekey   :: String  = "out",
                   verbose   :: Bool    = false,
                   silent    :: Bool    = false,
                  )

    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    env.transient = true

    # Arguments checking
    if !isnan(end_time)
        end_time > env.t || error("tm_solve! : end_time ($end_time) is greater that current time ($(env.t))")
        time_span = end_time - env.t
    end

    if !silent
        printstyled("Thermomechanical FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
        println("  from t=$(round(env.t,digits=4)) to t=$(round(env.t+time_span,digits=3))")
        sw = StopWatch() # timing
    end

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

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements and pw
    pmap  = nu+1:ndofs   # map for prescribed displacements and pw
    dom.ndofs = length(dofs)
    verbosity>0 && println("  unknown dofs: $nu")

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
            if dof.name==:ut
                dof.vals[:T] = env.T0 # real temperature
            end
        end

        update_loggers!(dom)  # Tracking nodes, ips, elements, etc.
        update_output_data!(dom) # Updates data arrays in domain
        complete_ut_T(dom)

        if save_incs
            save(dom, "$outdir/$filekey-0.vtu", verbose=false)
        end
    end

    # Incremental analysis
    t    = env.t     # current time
    tend = t + time_span # end time
    Δt = time_span/nincs # initial Δt value

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

    Fex  = zeros(ndofs)  # vector of external loads
    Uex  = zeros(ndofs)  # vector of external essential values

    Uex, Fex = get_bc_vals(dom, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # Get unbalanced forces
    Fin = zeros(ndofs)
    for elem in dom.elems
        elem_internal_forces(elem, Fin)
    end
    Fex .-= Fin # add negative forces to external forces vector

    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

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

        # Get forces and displacements from boundary conditions
        env.t = t + Δt
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
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            verbosity>1 && print("    assembling... \r")
            G, RHS = tm_mount_global_matrices(dom, ndofs, Δt)

            R .+= RHS

            # Solve
            verbosity>1 && print("    solving...   \r")
            tm_solve_step!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

            # Update
            verbosity>1 && print("    updating... \r")

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
            R .= ΔFex - ΔFin
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
                    dof.vals[:T] = U[i] + env.T0
                end
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Update time
            t += Δt
            T += ΔT

            # Check for saving output file
            if abs(T - Tout) < Ttol && save_incs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                complete_ut_T(dom)
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

        # Fix Δt according to ΔT
        Δt = ΔT*time_span
    end

    # time spent
    verbosity>1 && printstyled("  stage $(env.cstage) $(see(sw)) progress 100%\033[K\n", bold=true, color=:blue) # color 111
    verbosity==1 && println("  time spent: ", see(sw, format=:hms), "\033[K")

    update_output_data!(dom)
    save_incs || complete_ut_T(dom)

    return true

end
