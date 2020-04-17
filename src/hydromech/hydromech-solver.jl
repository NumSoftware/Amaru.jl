# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# Assemble the global stiffness matrix
function hm_mount_global_matrices(dom::Domain, ndofs::Int, Δt::Float64)

    # Assembling matrix G

    R, C, V = Int64[], Int64[], Float64[]
    RHS = zeros(ndofs)

    α = 1.0 # time integration factor

    for elem in dom.elems

        ty = typeof(elem)
        has_stiffness_matrix    = hasmethod(elem_stiffness, (ty,))
        has_coupling_matrix     = hasmethod(elem_coupling_matrix, (ty,))
        has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
        has_compressibility_matrix = hasmethod(elem_compressibility_matrix, (ty,))
        has_RHS_vector          = hasmethod(elem_RHS_vector, (ty,))

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
            Uw = [ dof.vals[:uw] for node in elem.nodes for dof in node.dofs if dof.name==:uw ]
            #@show Uw
            #@show H
            #@show elem.shape.name
            #@show elem.nodes[end].dofs
            RHS[rmap] -= Δt*(H*Uw)
        end

        # Assemble the compressibility matrix
        if has_compressibility_matrix
            Cpp, rmap, cmap =  elem_compressibility_matrix(elem)
            nr, nc = size(Cpp)
            for i=1:nr
                for j=1:nc
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, Cpp[i,j])
                end
            end
        end

        # Assemble ramaining RHS vectors
        if has_RHS_vector
            Q, map = elem_RHS_vector(elem)
            RHS[map] += Δt*Q
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
function hm_solve_step!(G::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, nu::Int)
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


function complete_uw_h(dom::Domain)
    haskey(dom.point_data, "uw") || return
    Uw = dom.point_data["uw"]
    H  = dom.point_data["h"]

    for elem in dom.elems
        elem.shape.family==SOLID_SHAPE || continue
        elem.shape==elem.shape.basic_shape && continue
        npoints  = elem.shape.npoints
        nbpoints = elem.shape.basic_shape.npoints
        map = [ elem.nodes[i].id for i=1:nbpoints ]
        Ue = Uw[map]
        He = H[map]
        C = elem.shape.nat_coords
        for i=nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i,:]
            N = elem.shape.basic_shape.func(R)
            Uw[id] = dot(N,Ue)
            H[id] = dot(N,He)
        end
    end
end


"""
    hm_solve!(domain, bcs, options...) -> Bool

Performs one stage finite element hydro-mechanical analysis of a `domain`
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
function hm_solve!(
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

    # Arguments checking
    verbosity = 1
    verbose && (verbosity=2)
    silent && (verbosity=0)

    tol>0 || error("hm_solve! : tolerance `tol `should be greater than zero")
    Ttol>0 || error("hm_solve! : tolerance `Ttol `should be greater than zero")

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    env.transient = true
    sw = StopWatch() # timing

    if !isnan(end_time)
        end_time > env.t || error("hm_solve! : end_time ($end_time) is greater that current time ($(env.t))")
        time_span = end_time - env.t
    end
    isnan(time_span) && error("hm_solve!: neither time_span nor end_time were set.")

    if verbosity==0
        printstyled("Hydromechanical FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
        println("  from t=$(round(dom.env.t,digits=4)) to t=$(round(dom.env.t+time_span,digits=3))")
    end

    verbosity>0 && println("  model type: ", env.modeltype)

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

        strip(outdir) == "" && (outdir = ".")
        isdir(outdir) || error("hm_solve!: output directory <$outdir> not fount")
        outdir[end] in ('/', '\\')  && (outdir = outdir[1:end-1])
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements and pw
    pmap  = nu+1:ndofs   # map for prescribed displacements and pw
    dom.ndofs = length(dofs)
    verbosity>0 && println("  unknown dofs: $nu")

    # get elevation Z for all Dofs
    Z = zeros(ndofs)
    for node in dom.nodes
        for dof in node.dofs
            Z[dof.eq_id] = node.X[env.ndim]
        end
    end

    # Get global parameters
    gammaw = get(dom.env.params, :gammaw, NaN)
    isnan(gammaw) && error("hm_solve!: gammaw parameter was not set in Domain")
    gammaw > 0 || error("hm_solve: invalid value for gammaw: $gammaw")

    # Save initial file and loggers
    if env.cstage==1
        # Setup initial quantities at dofs
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name==:uw
                dof.vals[:h] = 0.0 # water head
            end
        end

        update_output_data!(dom) # Updates data arrays in domain
        update_single_loggers!(dom)  # Tracking nodes, ips, elements, etc.
        update_composed_loggers!(dom)
        complete_uw_h(dom)

        if save_outs
            save(dom, "$outdir/$filekey-0.vtu", verbose=false)
        end
    end

    # Get the domain current state and backup
    State = [ ip.data for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

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
        maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        nfails    = 0  # counter for iteration fails
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            verbosity>1 && print("    assembling... \r")
            G, RHS = hm_mount_global_matrices(dom, ndofs, Δt)

            R .+= RHS

            # Solve
            verbosity>1 && print("    solving...   \r")
            hm_solve_step!(G, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

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
                if dof.name==:uw
                    dof.vals[:h] = Z[i] + U[i]/gammaw
                end
            end

            update_single_loggers!(dom) # Tracking nodes, ips, etc.

            # Update time
            t += Δt
            T += ΔT

            # Check for saving output file
            if abs(T - Tout) < Ttol && save_outs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                update_composed_loggers!(dom)
                complete_uw_h(dom)
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

        # Set Δt according to ΔT
        Δt = ΔT*time_span
    end

    if !save_outs
	update_output_data!(dom)
	complete_uw_h(dom)
	update_composed_loggers!(dom)
    end

    # time spent
    verbosity>1 && printstyled("  stage $(env.cstage) $(see(sw)) progress 100%\033[K\n", bold=true, color=:blue) # color 111
    verbosity==1 && println("  time spent: ", see(sw, format=:hms), "\033[K")

    return true
end
