# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export dynsolve!

# Assemble the global stiffness matrix
#=
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
=#


# Assemble the global mass matrix
function mount_M(dom::Domain, ndofs::Int)

    R, C, V = Int64[], Int64[], Float64[]

    for elem in dom.elems
        Me, rmap, cmap = elem_mass(elem)
        nr, nc = size(Me)
        for i=1:nr
            for j=1:nc
                push!(R, rmap[i])
                push!(C, cmap[j])
                push!(V, Me[i,j])
            end
        end
    end

    local M
    try
        M = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show ndofs
        @show err
    end

    return M
end





# Solves for a load/displacement increment
function solve_system!(K::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, nu::Int)
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



function frequencies(dom::Domain, bcs::Array; nmods::Int=5, savemods=true, verbose=true)
    if verbose
        print_with_color(:cyan,"FEM modal analysis:\n", bold=true) 
        tic()
    end

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs)

    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbose && println("  unknown dofs: $nu")

    K = mount_K(dom, ndofs)
    M = mount_M(dom, ndofs)

    if nu>0
        K11 = K[1:nu, 1:nu]
        M11 = M[1:nu, 1:nu]
    end

    # Inverse of the lumped matrix in vector form
    m11 = inv.(sum(M11, 2)) 
    P = m11.*K11

    # Eigenvalues and eigenvectors
    neign = min(20, length(m11)-2)
    L, V =eigs(P, nev=neign, which=:SM)


    p = sortperm(real(L))
    q = [ i for (i,v) in enumerate(L) if imag(v)!=0.0 ] 
    r = [ i for (i,v) in enumerate(L) if real(v)<=0.0 ] 
    pp = setdiff(p, union(q,r))

    L = L[pp][1:nmods].^0.5
    V = V[pp][1:nmods]

    idefmod = 1    
    for j=1:nmod

        Umod = V[:,j]

        for (k,dof) in enumerate(dofs)
            dof.vals[dof.name] = Umod[k] 
        end

        save(dom, dom.filekey * "-defmod$idefmod.vtk", verbose=false) # saves current domain state for modal deformated
        verbose && print_with_color(:green, "  ", dom.filekey * "-defmod$idefmod.vtk file written (Domain)\n")
        idefmod += 1

    end

    # reset displacement values
    for (k,dof) in enumerate(dofs)
        dof.vals[dof.name] = 0.0
    end

    return L, V

end

function rayleigh_coeffs(w1::Float64, w2::Float64, xi1::Float64, xi2::Float64=xi1)
    alpha = 2 * ((w1 ^ 2) * w2 * xi2 - w1 * (w2 ^ 2) * xi1) / ((w1 ^ 2)-(w2 ^ 2))
    beta  = 2 * (w1 * xi1 - w2 * xi2)/((w1 ^ 2) - (w2 ^ 2))
    return alpha, beta
end




"""
    solve!(D, bcs, options...) -> Bool

Performs one stage finite element analysis of a mechanical domain `D`
subjected to an array of boundary conditions `bcs`.

Available options are:

`verbose=true` : If true, provides information of the analysis steps

`tol=1e-2` : Tolerance for the absolute error in forces

`nincs=1` : Number of increments

`auto_inc=false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxits=5` : The maximum number of Newton-Rapson iterations per increment

`save_incs=false` : If true, saves output files according to `nouts` option

`nouts=0` : Number of output files per analysis

`scheme= :FE` : Predictor-corrector scheme at iterations. Available schemes are `:FE` and `:ME`

`saveips=false` : If true, saves corresponding output files with ip information

"""
function dynsolve!(dom::Domain, bcs::Array; time_span::Real=0.0, nincs::Int=1, maxits::Int=5, auto_inc::Bool=false, 
                   nouts::Int=0, nmods::Int=10, alpha::Real=0.0, beta::Real=0.0,
                   tol::Number=1e-2, verbose::Bool=true, save_incs::Bool=false, 
                   scheme::Symbol = :FE, save_ips::Bool=false, filekey="out")::Bool

    if verbose
        print_with_color(:cyan,"FEM dynamic analysis:\n", bold=true) 
        tic()
    end
    
    (save_incs && nouts==0) && (nouts=min(nincs,10))  # default value for nouts
    save_incs = nouts>0
    if save_incs
        if nouts>nincs
            nincs = nouts
        end
        if nincs%nouts != 0
            nincs = nincs - (nincs%nouts) + nouts
        end
        #info("  updating nincs to $nincs")
    end


    # Dictionary of data keys related with a dof
    components_dict = Dict(:ux => (:ux, :fx, :vx, :ax), 
                           :uy => (:uy, :fy, :vy, :ay),
                           :uz => (:uz, :fz, :vz, :az))


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
    #Uex, Fex = get_bc_vals(dom, bcs)

    # Global RHS vector 
    #RHS   = mount_RHS(dom, ndofs, 0.0)
    #Fex .+= RHS

    # Setup quantities at dofs
    if dom.nincs == 0
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end
    end



    # Backup the last converged state at ips. TODO: make backup to a vector of states
    for ip in ips
        ip.data0 = deepcopy(ip.data)
    end

    # Initial accelerations
    K = mount_K(dom, ndofs)
    M = mount_M(dom, ndofs)
    A = zeros(ndofs) 
    V = zeros(ndofs) 
    Uex, Fex = get_bc_vals(dom, bcs) # get values at time t
    solve_system!(M, A, Fex, nu)

    for (i,dof) in enumerate(dofs)
        us, fs, vs, as = components_dict[dof.name]
        dof.vals[vs] = V[i]
        dof.vals[as] = A[i]
        dof.vals[fs] = Fex[i]
    end

    update_loggers!(dom)  # Tracking nodes, ips, elements, etc.

    # Save initial file
    if save_incs
        save(dom, "$(dom.filekey)-0.vtk", verbose=false, save_ips=save_ips)
        verbose && print_with_color(:green, "  $(dom.filekey)-0.vtk file written (Domain)\n")
    end

    # Incremental analysis
    Dt = time_span
    t  = dom.shared_data.t
    tend = t + Dt # end time
    dt = Dt/nincs # initial dt value

    dT = Dt/nouts  # output time increment for saving vtk file
    T  = t + dT         # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 1       # increment counter
    iout = dom.nouts     # file output counter
    #F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    Fin  = zeros(ndofs)  # total internal force
    ΔFin = zeros(ndofs)  # internal forces for current increment
    ΔUa  = zeros(ndofs)  # essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # essential values for current iteration obtained from NR algorithm
    Fina = zeros(ndofs)  # current internal forces
    TFin = zeros(ndofs)  
    Aa   = zeros(ndofs)
    Va   = zeros(ndofs)
    #lastFin = zeros(ndofs)

    remountK = true



    while t < tend - ttol
        Uex, Fex = get_bc_vals(dom, bcs, t+dt) # get values at time t+dt

        verbose && print_with_color(:blue, "  increment $inc from t=$(round(t,10)) to t=$(round(t+dt,10)) (dt=$(round(dt,10))):", bold=true) # color 111
        verbose && println()
        #R   .= FexN - F    # residual
        Fex_Fin = Fex-Fin    # residual
        ΔUa .= 0.0
        ΔUi .= Uex    # essential values at iteration i # TODO: check for prescribed disps

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
            verbose && print("    assembling K... \r")
            remountK && (K = mount_K(dom, ndofs))

            C   = alpha*M + beta*K # Damping matrix
            Kp  = K + (4/(dt^2))*M + (2/dt)*C # pseudo-stiffness matrix
            ΔFp = Fex_Fin + M*(A + 4*V/dt + 4*(-ΔUa)/(dt^2)) + C*(V + (2*(-ΔUa)/dt))

            # Solve
            verbose && print("    solving...   \r")
            solve_system!(Kp, ΔUi, ΔFp, nu)

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

            Fina = Fin + ΔFin

            # Update V and A
            # Optional (V and A may be used from the interval beginning)
            Va = -V + 2*(ΔUa + ΔUi)/dt;
            Aa = -A + 4*((ΔUa + ΔUi) - V*dt)/(dt^2);               


            TFin = Fina + C*Va + M*Aa  # Internal force including dynamic effects 

            residue = maximum(abs, (Fex-TFin)[umap] ) 

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            Fex_Fin .= Fex .- Fina  # Check this variable, it is not the residue actually
            Fex_Fin[pmap] .= 0.0  # Zero at prescribed positions

            if verbose
                print_with_color(:bold, "    it $it  ")
                @printf(" residue: %-10.4e\n", residue)
            end

            if residue > tol; Fina -= ΔFin end
            if residue < tol;        converged = true ; remountK=false; break end
            if isnan(residue);       converged = false; break end
            if it > maxits;          converged = false; break end
            if residue > 0.9*lastres;  nfails += 1 end
            if nfails == maxfails;     converged = false; break end
        end

        if converged
            # Update forces and displacement for the current stage
            Fin = Fina
            U .+= ΔUa

            # Backup converged state at ips
            for ip in ips; ip.data0 = deepcopy(ip.data) end

            # Update vectors for velocity and acceleration
            A = Aa
            V = Va

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                us, fs, vs, as = components_dict[dof.name]
                dof.vals[us] = U[i]
                dof.vals[fs] = TFin[i]
                dof.vals[vs] = V[i]
                dof.vals[as] = A[i]
            end

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Update time t and dt
            inc += 1
            t   += dt

            # Check for saving output file
            if abs(t - T) < ttol
                iout += 1
                save(dom, "$(dom.filekey)-$iout.vtk", verbose=false, save_ips=save_ips)
                T += dT # find the next output time
                verbose && print_with_color(:green, "  $(dom.filekey)-$iout.vtk file written (Domain)\n")
            end


            if auto_inc
                dt = min(1.5*dt, Dt/nincs)
                dt = round(dt, -ceil(Int, log10(dt))+3)  # round to 3 significant digits
            end

            # Fix dt in case d+dt>T
            if t+dt>T
                dt = T-t
            end
        else
            if auto_inc
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
