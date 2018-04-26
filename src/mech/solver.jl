# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export solve!

# Assemble the global stiffness matrix
function mount_K(dom::Domain, ndofs::Int)

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


function mount_RHS(dom::Domain, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in dom.elems
        F, map = elem_RHS(elem) # RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = F
    end
    return RHS
end


# Solves for a load/displacement increment
function solve_step(K::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, nu::Int)
    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    ndofs = length(DU)
    umap  = 1:nu
    pmap  = nu+1:ndofs
    if nu == ndofs 
        warn("solve!: No essential boundary conditions.")
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
            LUfact = lufact(K11)
            U1  = LUfact\RHS
            F2 += K21*U1
        catch err
            warn("solve!: $err")
            U1 .= NaN
        end
    end

    # Completing vectors
    DU[1:nu]     .= U1
    DF[nu+1:end] .= F2
end

"""
    solve!(D, bcs, options...) -> Bool

Performs one stage finite element analysis of a mechanical domain `D`
subjected to an array of boundary conditions `bcs`.

Available options are:

`verbose=true` : If true, provides information of the analysis steps

`tol=1e-2` : Tolerance for the absolute error in forces

`nincs=1` : Number of increments

`autoinc=false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxits=5` : The maximum number of Newton-Rapson iterations per increment

`save_incs=false` : If true, saves output files according to `nouts` option

`nouts=0` : Number of output files per analysis

`scheme= :FE` : Predictor-corrector scheme at iterations. Available schemes are `:FE` and `:ME`

`saveips=false` : If true, saves corresponding output files with ip information

"""
function solve!(dom::Domain, bcs::Array; nincs=1::Int, maxits::Int=5, autoinc::Bool=false, maxincs::Int=0,
    tol::Number=1e-2, verbose::Bool=true, nouts::Int=0,
    scheme::Symbol = :FE, save_ips::Bool=false)::Bool

    if verbose
        print_with_color(:cyan,"FEM analysis:\n", bold=true) 
        tic()
    end
    
    save_incs = nouts>0
    if save_incs
        if nouts>nincs
            nincs = nouts
            info("  nincs changed to $nincs to match nouts")
        end
        if nincs%nouts != 0
            nincs = nincs - (nincs%nouts) + nouts
            info("  nincs changed to $nincs to be a multiple of nouts")
        end
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs)
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbose && println("  unknown dofs: $nu")
    
    # Get array with all integration points
    ips = [ ip for elem in dom.elems for ip in elem.ips ]

    # Get forces and displacements from boundary conditions
    Uex, Fex = get_bc_vals(dom, bcs)

    # Global RHS vector 
    RHS   = mount_RHS(dom, ndofs, 0.0)
    Fex .+= RHS

    # Setup quantities at dofs
    if dom.nincs == 0
        for (i,dof) in enumerate(dofs)
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end
    end

    update_loggers!(dom)  # Tracking nodes, ips, elements, etc.

    # Save initial file
    if dom.nincs == 0 && save_incs 
        save(dom, "$(dom.filekey)-0.vtk", verbose=false, save_ips=save_ips)
        verbose && print_with_color(:green, "  $(dom.filekey)-0.vtk file written (Domain)\n")
    end

    # Backup the last converged state at ips. TODO: make backup to a vector of states
    for ip in ips
        ip.data0 = deepcopy(ip.data)
    end

    # Incremental analysis
    t  = 0.0
    dt = 1.0/nincs # initial dt value
    dt_bk = 0.0

    dT = 1.0/nouts  # output time increment for saving vtk file
    T  = dT        # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 1       # increment counter
    iout = dom.nouts     # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration

    remountK = true

    while t < 1.0 - ttol
        verbose && print_with_color(:blue, "  increment $inc from t=$(round(t,10)) to t=$(round(t+dt,10)) (dt=$(round(dt,10))):", bold=true) # color 111
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
        local K::SparseMatrixCSC{Float64,Int64}
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            if it>1; remountK=true end 
            lastres = residue # residue from last iteration

            # Try FE step
            verbose && print("    assembling... \r")
            remountK && (K = mount_K(dom, ndofs))

            # Solve
            verbose && print("    solving...   \r")
            solve_step(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R

            # Update
            verbose && print("    updating... \r")

            # Restore the state to last converged increment
            for ip in ips; ip.data = deepcopy(ip.data0) end

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
                for ip in ips; ip.data = deepcopy(ip.data0) end

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
            R[pmap] .= 0.0  # Zero at prescribed positions

            if verbose
                print_with_color(:bold, "    it $it  ")
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
            for ip in ips; ip.data0 = deepcopy(ip.data) end

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Update time t and dt
            inc += 1
            t   += dt

            if maxincs != 0 && inc > maxincs
                print_with_color(:red, "  solver maxincs = $maxincs reached (try maxincs=0)\n")
                break
            end

            # Check for saving output file
            if abs(t - T) < ttol
                iout += 1
                save(dom, "$(dom.filekey)-$iout.vtk", verbose=false, save_ips=save_ips)
                T += dT # find the next output time
                verbose && print_with_color(:green, "  $(dom.filekey)-$iout.vtk file written (Domain)\n")
            end

            if autoinc
                if dt_bk>0.0
                    dt = dt_bk
                else
                    dt = min(1.5*dt, 1.0/nincs)
                    #dt = 1.5*dt
                    dt = round(dt, -ceil(Int, log10(dt))+3)  # round to 3 significant digits
                end
            end
            dt_bk = 0.0

            # Fix dt in case d+dt>T
            if t+dt>T
                dt_bk = dt
                dt = T-t
            end
        else
            if autoinc
                verbose && println("    increment failed.")
                dt *= 0.5
                dt = round(dt, -ceil(Int, log10(dt))+3)  # round to 3 significant digits
                if dt < ttol
                    print_with_color(:red, "solve!: solver did not converge\n",)
                    return false
                end
            else
                print_with_color(:red, "solve!: solver did not converge\n",)
                return false
            end
        end
    end

    # time spent
    if verbose
        h, r = divrem(toq(), 3600)
        m, r = divrem(r, 60)
        println("  time spent: $(h)h $(m)m $(round(r,3))s")
    end

    # Update number of used increments at domain
    dom.nincs += inc
    dom.nouts = iout

    return true

end