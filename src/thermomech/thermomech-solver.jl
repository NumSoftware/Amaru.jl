# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

#export solve!
# Assemble the global stiffness matrix
function tm_mount_G_RHS(dom::Domain, ndofs::Int, Δt::Float64)

    # Assembling matrix G

    R, C, V = Int64[], Int64[], Float64[]
    RHS = zeros(ndofs)

    α = 1.0 # time integration factor

    for elem in dom.elems

        ty = typeof(elem)
        has_stiffness_matrix    = hasmethod(elem_stiffness, (ty,))
        has_coupling_matrix     = hasmethod(elem_coupling_matrix, (ty,))
        has_mass_matrix = hasmethod(elem_mass_matrix, (ty,))
        has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
        has_RHS_vector          = hasmethod(elem_RHS_vector, (ty,))

        # Assemble the stiffness matrix
        if has_stiffness_matrix
            K, rmap, cmap = elem_stiffness_matrix(elem)
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
            H, rmap, cmap = elem_conductivity_matrix(elem)
            nr, nc = size(H)
            for i=1:nr
                for j=1:nc
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, α*Δt*H[i,j])
                end
            end

            # Assembling RHS components
            Ut = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes ]
            RHS[rmap] -= Δt*(H*Ut)
        end

        # Assemble the conductivity matrix
        if has_mass_matrix
            M, rmap, cmap = elem_mass_matrix(elem)
            nr, nc = size(M)
            for i=1:nr
                for j=1:nc
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, M[i,j])
                end
            end
        end

        # Assembling RHS components
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

        function complete_ut_h(dom::Domain)
            haskey(dom.point_data, "ut") || return
            Ut = dom.point_data["ut"]
            H  = dom.point_data["h"]
            for ele in dom.elems
                ele.shape.family==SOLID_SHAPE || continue
                ele.shape==ele.shape.basic_shape && continue
                npoints = ele.shape.npoints
                nbpoints = ele.shape.basic_shape.npoints
                map = [ ele.nodes[i].id for i=1:nbpoints ]
                Ue = Ut[map]
                He = H[map]
                C = ele.shape.nat_coords
                for i=nbpoints+1:npoints
                    id = ele.nodes[i].id
                    R = C[i,:]
                    N = ele.shape.basic_shape.func(R)
                    Ut[id] = dot(N,Ue)
                    H[id] = dot(N,He)
                end
            end
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

    # get elevation Z for all Dofs
    #Z = zeros(ndofs)
    #for node in dom.nodes
    #    for dof in node.dofs
    #        Z[dof.eq_id] = node.X[env.ndim]
    #    end
#    end

# Get global parameters
#gammaw = get(dom.env.params, :gammaw, NaN)
#isnan(gammaw) && error("hm_solve!: gammaw parameter was not set in Domain")
#gammaw > 0 || error("hm_solve: invalid value for gammaw: $gammaw")

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
            #if dof.name==:uw
        #    dof.vals[:h] = 0.0 # water head
            #end
        end

        update_loggers!(dom)  # Tracking nodes, ips, elements, etc.
        update_output_data!(dom) # Updates data arrays in domain
        #complete_uw_h(dom)

        if save_incs
            save(dom, "$outdir/$filekey-0.vtk", verbose=false)
            verbose && printstyled("  $outdir/$filekey-0.vtk file written (Domain)\n", color=:green)
        end
    end

    # Incremental analysis
    t    = dom.env.t # current time
    tend = t + time_span  # end time
    Δt = time_span/nincs # initial Δt value

    dT = time_span/nouts  # output time increment for saving vtk file
    T  = t + dT        # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 0       # increment counter
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

    # Get unbalanced forces
    Fin = zeros(ndofs)
    for elem in dom.elems
        elem_internal_forces(elem, Fin)
    end
    Fex .-= Fin # add negative forces to external forces vecto

    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while t < tend - ttol
        # Update counters
        inc += 1
        env.cinc += 1

        if inc > maxincs
            printstyled("  solver maxincs = $maxincs reached (try maxincs=0)\n", color=:red)
            return false
        end

        silent || printstyled("  increment $inc from t=$(round(t,sigdigits=9)) to t=$(round(t+Δt,sigdigits=9)) (dt=$(round(Δt,sigdigits=9))):"," "^10,"\r", bold=true, color=:blue) # color 111
        #silent || printstyled("  increment $inc  (progress=$(round), dt=$(round(Δt,sigdigits=9))):"," "^10,"\r", bold=true, color=:blue) # color 111
        verbose && println()

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
        #        if dof.name==:uw
            #        dof.vals[:h] = Z[i] + U[i]/gammaw
    #            end
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Check for saving output file
            Tn = t + Δt

            if Tn+ttol>=T && saveincs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
        #       complete_uw_h(dom)
                save(dom, "$outdir/$filekey-$iout.vtk", verbose=false)
                T = Tn - mod(Tn, dT) + dT
                silent || verbose || print(" "^70, "\r")
                silent || printstyled("  $outdir/$filekey-$iout.vtk file written (Domain) \033[K \n",color=:green)
            end

            # Update time t
            t   += Δt

            # Get new Δt
            if autoinc
                Δt = min(1.5*Δt, 1.0/nincs)
                Δt = round(Δt, sigdigits=3)
                Δt = min(Δt, tend-t)
            end
        else
            # Restore counters
            inc -= 1
            env.cinc -= 1

            # Restore the state to last converged increment
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
    #complete_uw_h(dom)

    return true

end
