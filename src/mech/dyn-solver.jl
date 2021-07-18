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
function dym_solve_system!(K::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, nu::Int)
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
    dynsolve!(dom, bcs, options...) :: Bool

Performs one stage finite element dynamic analysis of a mechanical domain `dom`
subjected to a set of boundary conditions `bcs` and a time span.


# Arguments

`dom` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`time_span` = 0.0 : Simulated time span

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`maxincs = 1000000` : Maximum number of increments

`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`filekey = ""` : File key for output files

`verbose = true` : If true, provides information of the analysis steps

"""
function dynsolve!(
                   dom       :: Domain,
                   bcs       :: Array;
                   time_span :: Real    = NaN,
                   end_time  :: Real    = NaN,
                   nincs     :: Int     = 1,
                   maxits    :: Int     = 5,
                   autoinc   :: Bool    = false,
                   tol       :: Number  = 1e-2,
                   Ttol      :: Number  = 1e-9,
                   nouts     :: Int     = 0,
                   sism      :: Bool    = false,
                   tds       :: Float64 = 0.0,
                   tss       :: Float64 = 0.0,
                   nmods     :: Int     = 10,
                   alpha     :: Real    = 0.0,
                   beta      :: Real    = 0.0,
                   outdir    :: String  = "",
                   filekey   :: String  = "out",
                   printlog = false,
                   verbose = false
                  )

    # Arguments checking
    verbosity = 0
    printlog && (printlog=false)
    printlog && verbose && (verbosity=2)

    tol>0 || error("solve! : tolerance should be greater than zero")
    Ttol>0 || error("solve! : tolerance `Ttol `should be greater than zero")

    env = dom.env
    env.cstage += 1
    env.cinc    = 0
    sw = StopWatch() # timing

    if !isnan(end_time)
        end_time > env.t || error("dynsolve!! : end_time ($end_time) is greater that current time ($(env.t))")
        time_span = end_time - env.t
    end
    isnan(time_span) && error("dynsolve!: neither time_span nor end_time were set.")

    verbosity>0 && println("  model type: ", env.modeltype)

    if verbosity==0
        printstyled("Dynamic FE analysis: Stage $(env.cstage)\n", bold=true, color=:cyan)
    end

    tol>0 || error("solve! : tolerance should be greater than zero.")
    verbosity>0 || println("  model type: ", env.modeltype)

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
        isdir(outdir) || error("solve!: output directory <$outdir> not fount")
        outdir[end] in ('/', '\\')  && (outdir = outdir[1:end-1])
    end

    # Dictionary of data keys related with a dof
    components_dict = Dict(:ux => (:ux, :fx, :vx, :ax),
                           :uy => (:uy, :fy, :vy, :ay),
                           :uz => (:uz, :fz, :vz, :az),
                           :rx => (:rx, :mx, :vrx, :arx),
                           :ry => (:ry, :my, :vry, :ary),
                           :rz => (:rz, :mz, :vrz, :arz))

    # Set model environment as transient
    dom.env.transient = true

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(dom, bcs)
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    dom.ndofs = length(dofs)
    verbosity>0 && println("  unknown dofs: $nu")

    # Setup quantities at dofs
    if env.cstage==1
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in dom.elems for ip in elem.ips ]
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
    K = mount_K(dom, ndofs, verbosity)
    M = mount_M(dom, ndofs)
    A = zeros(ndofs)
    V = zeros(ndofs)
    Uex, Fex = get_bc_vals(dom, bcs) # get values at time t

    #If the problem has a sism, the force sismic is add
    if sism && tss<=0 #tss:time when seismic activity starts tds: time of seismic duration 0:current time =0s
        #M = mount_M(dom,ndofs)
        Fex = sismic_force(dom, bcs, M, Fex, AS, keysis, 0.0, tds)
    end
    dym_solve_system!(M, A, Fex, nu)

    # Initial values at nodes
    for (i,dof) in enumerate(dofs)
        us, fs, vs, as = components_dict[dof.name]
        dof.vals[vs] = V[i]
        dof.vals[as] = A[i]
        dof.vals[fs] = Fex[i]
    end

    # Save initial file
    if env.cstage==1 && save_outs
        update_output_data!(dom)
        update_single_loggers!(dom)
        update_composed_loggers!(dom)
        save(dom, "$outdir/$filekey-0.vtu", printlog=true)
    end

    # Incremental analysis
    Dt = time_span
    t  = dom.env.t
    tend = t + Dt # end time
    dt = Dt/nincs # initial dt value

    dT = Dt/nouts  # output time increment for saving output file
    T  = t + dT         # output time for saving the next output file

    ttol = 1e-9    # time tolerance
    inc  = 0       # increment counter
    iout = env.cout      # file output counter
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
        # Update counters
        inc += 1
        env.cinc += 1

        Uex, Fex = get_bc_vals(dom, bcs, t+dt) # get values at time t+dt

        # If the problem has a sism, the force sismic is added
        if sism && tss<=t+dt && tds+tss>=t+dt
            M = mount_M(dom,ndofs)
            Fex = sismic_force(dom, bcs, M,Fex,AS,keysis,t+dt,tds)
        end

        verbosity>0 || printstyled("  stage $(env.cstage)  increment $inc from t=$(round(t,digits=10)) to t=$(round(t+dt,digits=10)) (dt=$(round(dt,digits=10)))\e[K\r", bold=true, color=:blue) # color 111
        verbosity>0 && println()
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
            #verbosity>0 && print("    assembling K... \r")
            #remountK && (K = mount_K(dom, ndofs))
            K = mount_K(dom, ndofs, verbosity)

            C   = alpha*M + beta*K # Damping matrix
            Kp  = K + (4/(dt^2))*M + (2/dt)*C # pseudo-stiffness matrix
            ΔFp = Fex_Fin + M*(A + 4*V/dt + 4*(-ΔUa)/(dt^2)) + C*(V + (2*(-ΔUa)/dt))

            # Solve
            verbosity>0 && print("    solving...   \r")
            dym_solve_system!(Kp, ΔUi, ΔFp, nu)

            # Update
            verbosity>0 && print("    updating... \r")

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

            if verbosity>0 
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

            update_single_loggers!(dom)

            # Update time t and dt
            t   += dt
            dom.env.t = t

            # Check for saving output file
            if abs(t - T) < ttol && save_outs
                env.cout += 1
                iout = env.cout
                update_output_data!(dom)
                update_composed_loggers!(dom)
                save(dom, "$outdir/$filekey-$iout.vtu", printlog=true)
                T += dT # find the next output time
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
            # Restore counters
            env.cinc -= 1
            inc -= 1

            if autoinc
                verbosity>0 && println("    increment failed.")
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

    if !save_outs
        update_output_data!(dom)
        update_composed_loggers!(dom)
    end

    # time spent
    verbosity==1 && println("  time spent: ", see(sw, format=:hms), "\e[K")

    return true
end
