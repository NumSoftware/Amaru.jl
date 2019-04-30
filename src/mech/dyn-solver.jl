# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export dynsolve!


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


function sismic_force(dom::Domain, bcs, M::SparseMatrixCSC{Float64, Int}, F::Vect, AS::Array{Float64,2}, keysis::Symbol, tfia::Float64, tds::Float64)

    dofs, nu = configure_dofs!(dom, bcs)
    ndofs = length(dofs)

    ndat = length(AS)    # quantity of aceleration data
    AS2  = zeros(ndat+1) # sismic aceleration data how vector
    c = 0
    for i=1:size(AS,1)
        for j=1:size(AS,2)
            c += 1
            AS2[c+1] = AS[i,j]
        end
    end

    vts = zeros(ndat+1) # time vetor correspond to acelerations
    
    for i=1:ndat+1
        vts[i] = (i-1)*tds/ndat
    end
    
    FAS = hcat(vts,AS2) # Function of aceleration

    # Interpolation of the aceleration value

    inic = 0
    fin = 0

    for i=1:ndat
        if FAS[i,1]<=tfia<=FAS[i+1,1]
            inic = i
            fin = i+1
        end
        if inic!=0 & fin!=0; break end
    end

    m = (FAS[fin,2]-FAS[inic,2])/(FAS[fin,1]-FAS[inic,1])
    acel = FAS[inic,2] + m*(tfia - FAS[inic,1])

    #Dof sismic aceleration    

    VAS  = zeros(ndofs) #Sismic aceleration vector accord dof
    
    for node in dom.nodes
        dof = node.dofdict[keysis]
        VAS[dof.eq_id] += acel
    end    
    
    #Dof sismic force
    FS = M*VAS
    F += FS

    return F

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
function dynsolve!(dom::Domain, bcs::Array; time_span::Real=0.0, sism=false, tds::Float64=0.0, tss::Float64=0.0, 
                   nincs::Int=1, maxits::Int=5, autoinc::Bool=false, 
                   nouts::Int=0, nmods::Int=10, alpha::Real=0.0, beta::Real=0.0,
                   tol::Number=1e-2, verbose::Bool=true, save_incs::Bool=false, 
                   scheme::Symbol = :FE, filekey="out")::Bool

    if verbose
        printstyled("FEM dynamic analysis:\n", bold=true, color=:cyan)
        tic = time()
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

    # Set model environment as transient
    dom.env.transient = true

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs)
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbose && println("  unknown dofs: $nu")
    
    # Get array with all integration points
    ips = [ ip for elem in dom.elems for ip in elem.ips ]

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

    # Get the domain current state and backup
    State = [ ip.data for elem in dom.elems for ip in elem.ips ]
    StateBk = copy.(State)

    #If the problem is sismic, read the sismic acelerations asking to user the file's name AS:SeismicAcelerations
    if sism==true
         print("What is the .dat file name of the sismic acelerations?")
         AS = readdlm("$(chomp(readline())).dat")
         AS= 9.81*AS
         print("What is the key correspond to sismic direction (fx, fy, fz)?")
         keysis = Symbol(readline())
    end 

    # Timing
    sw = StopWatch()

    # Initial accelerations
    K = mount_K(dom, ndofs)
    M = mount_M(dom, ndofs)
    A = zeros(ndofs) 
    V = zeros(ndofs) 
    Uex, Fex = get_bc_vals(dom, bcs) # get values at time t
                
    #If the problem has a sism, the force sismic is add
    if sism && tss<=0 #tss:time when seismic activity starts tds: time of seismic duration 0:current time =0s
        #M = mount_M(dom,ndofs)
        Fex = sismic_force(dom, bcs, M, Fex, AS, keysis, 0.0, tds)
    end                
    solve_system!(M, A, Fex, nu)

    # Initial values at nodes
    for (i,dof) in enumerate(dofs)
        us, fs, vs, as = components_dict[dof.name]
        dof.vals[vs] = V[i]
        dof.vals[as] = A[i]
        dof.vals[fs] = Fex[i]
    end
    update_loggers!(dom)  # Tracking nodes, ips, elements, etc.

    # Save initial file
    if save_incs
        save(dom, "$(dom.filekey)-0.vtk", verbose=false)
        verbose && printstyled("  $(dom.filekey)-0.vtk file written (Domain)\n", color=:green)
    end

    # Incremental analysis
    Dt = time_span
    t  = dom.env.t
    tend = t + Dt # end time
    dt = Dt/nincs # initial dt value

    dT = Dt/nouts  # output time increment for saving vtk file
    T  = t + dT         # output time for saving the next vtk file

    ttol = 1e-9    # time tolerance
    inc  = 1       # increment counter
    iout = dom.nouts     # file output counter
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

    #remountK = true

    while t < tend - ttol
        Uex, Fex = get_bc_vals(dom, bcs, t+dt) # get values at time t+dt
                    
        # If the problem has a sism, the force sismic is added
        if sism && tss<=t+dt && tds+tss>=t+dt
            M = mount_M(dom,ndofs)
            Fex = sismic_force(dom, bcs, M,Fex,AS,keysis,t+dt,tds)
        end
                    

        verbose && printstyled("  increment $inc from t=$(round(t,digits=10)) to t=$(round(t+dt,digits=10)) (dt=$(round(dt,digits=10))):", bold=true, color=:blue) # color 111
        verbose && println()
        #R   .= FexN - F    # residual
        Fex_Fin = Fex-Fina    # residual
        ΔUa .= 0.0
        ΔUi .= Uex    # essential values at iteration i # TODO: check for prescribed disps

        # Newton Rapshon iterations
        residue   = 0.0
        converged = false
        maxfails  = 3    # maximum number of it. fails with residual change less than 90%
        nfails    = 0    # counter for iteration fails
        #local K::SparseMatrixCSC{Float64,Int64}

        for it=1:maxits
            if it>1; ΔUi .= 0.0 end # essential values are applied only at first iteration
            #if it>1; remountK=true end 
            lastres = residue # residue from last iteration

            # Try FE step
            verbose && print("    assembling K... \r")
            #remountK && (K = mount_K(dom, ndofs))
            K = mount_K(dom, ndofs)

            C   = alpha*M + beta*K # Damping matrix
            Kp  = K + (4/(dt^2))*M + (2/dt)*C # pseudo-stiffness matrix
            ΔFp = Fex_Fin + M*(A + 4*V/dt + 4*(-ΔUa)/(dt^2)) + C*(V + (2*(-ΔUa)/dt))

            # Solve
            verbose && print("    solving...   \r")
            solve_system!(Kp, ΔUi, ΔFp, nu)

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
                printstyled("    it $it  ", bold=true)
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
            copyto!.(StateBk, State)

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


            # Update time t and dt
            inc += 1
            t   += dt
            dom.env.t = t

            update_loggers!(dom) # Tracking nodes, ips, elements, etc.

            # Check for saving output file
            if abs(t - T) < ttol
                iout += 1
                save(dom, "$(dom.filekey)-$iout.vtk", verbose=false)
                T += dT # find the next output time
                verbose && printstyled("  $(dom.filekey)-$iout.vtk file written (Domain)\n", color=:green)
            end


            if autoinc
                dt = min(1.5*dt, Dt/nincs)
                dt = round(dt, digits=-ceil(Int, log10(dt))+3)  # round to 3 significant digits
            end

            # Fix dt in case d+dt>T
            if t+dt>T
                dt = T-t
            end
        else
            if autoinc
                verbose && println("    increment failed.")
                dt *= 0.5
                dt = round(dt, digits=-ceil(Int, log10(dt))+3)  # round to 3 significant digits
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

    # Update number of used increments at domain
    dom.nincs += inc
    dom.nouts = iout

    return true

end
