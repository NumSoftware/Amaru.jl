# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export AcousticMechAnalysis

mutable struct AcousticMechAnalysisProps<:Analysis
    stressmodel::String # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    
    function AcousticMechAnalysisProps(;stressmodel="3d", thickness=1.0, g=0.0)
        @check stressmodel in ("plane-stress", "plane-strain", "axisymmetric", "3d")
        @check thickness>0
        @check g>=0
        return new(stressmodel, thickness, g)
    end
end

AcousticMechAnalysis = AcousticMechAnalysisProps

function am_mount_M(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_mass, (ty,)) || continue

            Me, rmap, cmap = elem_acoustic_mass(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end


function am_mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_stiffness, (ty,)) || continue

            Ke, rmap, cmap = elem_acoustic_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


# function complete_ut_T(model::Model)
#     haskey(model.node_data, "ut") || return
#     Ut = model.node_data["ut"]
#     T0 = get(model.env.mat, :T0, 0.0)

#     for elem in model.elems
#         elem.shape.family==BULKCELL || continue
#         elem.shape==elem.shape.basic_shape && continue
#         npoints  = elem.shape.npoints
#         nbpoints = elem.shape.basic_shape.npoints
#         map = [ elem.nodes[i].id for i in 1:nbpoints ]
#         Ute = Ut[map]
#         C = elem.shape.nat_coords
#         for i in nbpoints+1:npoints
#             id = elem.nodes[i].id
#             R = C[i,:]
#             N = elem.shape.basic_shape.func(R)
#             Ut[id] = dot(N,Ute)
#         end
#     end

#     model.node_data["T"] = Ut .+ T0
# end


"""
    am_solve!(D, bcs, options...) -> Bool

Performs one stage finite element thermo-mechanical analysis of a `domain`
subjected to a list of boundary conditions `bcs`.

# Arguments

`model` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`time_span = NaN` : Time lapse for the transient analysis in the current stage

`end_time = NaN` : Final time for the transient analysis in the current stage

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`


`tol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`Tmin     = 1e-8` : Pseudo-time tolerance

`scheme  = :FE` : Predictor-corrector scheme at each increment. Available scheme is `:FE`

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`outkey = ""` : File key for output files

`quiet = false` : verbosity level from 0 (silent) to 2 (verbose)

"""
function solve!(model::Model, ana::AcousticMechAnalysis; args...)
    name = "Solver for acoustic mechanical analyses"
    status = stage_iterator!(name, am_stage_solver!, model; args...)
    return status
end


function am_stage_solver!(model::Model, stage::Stage; 
    tol     :: Number  = 1e-2,
    rtol    :: Number  = 1e-2,
    Tmin    :: Number  = 1e-9,
    rspan   :: Number  = 1e-2,
    maxits  :: Int     = 5,
    autoinc :: Bool    = false,
    outdir  :: String  = ".",
    outkey  :: String  = "out",
    quiet  :: Bool    = false
    )

    env = model.env

    println(env.log, "Acoustic FE analysis: Stage $(stage.id)")
    stage.status = :solving

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    tspan     = stage.tspan
    env       = model.env
    save_outs = stage.nouts > 0
    ftol      = tol

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    # global DOFS = dofs
    # global CC = stage.bcs
    # error()


    ndofs    = length(dofs)
    umap     = 1:nu         # map for unknown bcs
    pmap     = nu+1:ndofs   # map for prescribed bcs
    model.ndofs = length(dofs)
    println(env.log, "unknown dofs: $nu")
    # println(env.alerts, "  unknown dofs: $nu") # TODO

    # quiet || nu==ndofs && println(env.alerts, "solve_system!: No essential boundary conditions", Base.warn_color) #TODO

    # Dictionary of data keys related with a dof
    components_dict = Dict( 
                            :up => (:up, :fq, :vp, :ap),
                            :ux => (:ux, :fx, :vx, :ax),
                            :uy => (:uy, :fy, :vy, :ay),
                            :uz => (:uz, :fz, :vz, :az),
                            :rx => (:rx, :mx, :vrx, :arx),
                            :ry => (:ry, :my, :vry, :ary),
                            :rz => (:rz, :mz, :vrz, :arz))
    
    model.env.transient = true

    # Setup quantities at dofs
    if stage.id == 1
        for dof in dofs
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end

        # Initial accelerations
        A = zeros(ndofs)
        V = zeros(ndofs)
        Uex, Fexa = get_bc_vals(model, bcs) # get values at time t
        
        M = am_mount_M(model.elems, ndofs)

        # initial acceleration
        sysstatus = solve_system!(M, A, Fexa, nu)
        # failed(sysstatus) && return failure("solver: Reduced mass matrix is singular")


        # Initial values at nodes
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[vs] = V[i]
            dof.vals[as] = A[i]
            dof.vals[fs] = Fexa[i]
        end

        # Save initial file and loggers
        update_output_data!(model)
        update_single_loggers!(model)
        update_multiloggers!(model)
        update_monitors!(model)
        # complete_ut_T(model)
        save_outs && save(model, "$outdir/$outkey-0.vtu", quiet=true)
    end
    lastflush = time()
    flushinterval = 5.0

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    
    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, 0.01))

    ΔTbk    = 0.0
    ΔTcheck = save_outs ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    t    = env.t     # current time
    
    inc  = 0             # increment counter
    iout = env.out       # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUi  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUk  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    Fina = zeros(ndofs)  # current internal forces
    
    
    ΔFexi = zeros(ndofs)  # increment of external natural bc
    ΔUexi = zeros(ndofs)  # increment of external essential bc
    Fexi! = zeros(ndofs)  # last value of external natural bc
    Uexi! = zeros(ndofs)  # last value of external essential bc

    sysstatus = ReturnStatus()

    # Get boundary conditions
    # Uex, Fexa = get_bc_vals(model, bcs, t) # get values at time t  #TODO pick internal forces and displacements instead!

    # for elem in model.elems
    #     elem_internal_forces(elem, Fin)
    # end
    # Fexa .-= Fin # add negative forces to external forces vector

    # Get global vectors from values at dofs
    for (i,dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64,Int64}
    local RHS::Array{Float64,1}

    while T < 1.0-Tmin
        env.ΔT = ΔT
        Δt     = tspan*ΔT
        env.t  = t + Δt

        # Update counters
        inc += 1
        env.inc = inc 

        println(env.log, "  inc $inc")

        # Get forces and displacements from boundary conditions
        Uexi, Fexi = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt
        ΔUexi = Uexi - Uexi!
        ΔFexi = Fexi - Fexi!

        ΔUexi[umap] .= 0.0
        ΔFexi[pmap] .= 0.0

        ΔF = ΔFexi
        ΔF[pmap] .= 0.0  # Zero at prescribed positions

        ΔUi .= 0.0
        ΔUk .= ΔUexi    # essential values at iteration i

        # Newton Rapshon iterations
        res   = 0.0
        # maxfails  = 3  # maximum number of it. fails with residual change less than 90%
        # nfails    = 0  # counter for iteration fails
        nits  = 0
        res1  = 0.0
        converged = false
        syserror   = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUk.=0.0) # essential values are applied only at first iteration
            lastres = res # res from last iteration

            K = am_mount_K(model.elems, ndofs)
            M = am_mount_M(model.elems, ndofs)

            K′ = K + 4/Δt^2*M   # pseudo-stiffness matrix
            ΔF′ = ΔF + M*(A + 4*V/Δt - 4*ΔUi/Δt^2)

            # Solve
            solve_system!(K′, ΔUk, ΔF′, nu)

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUit   = ΔUi + ΔUk
            ΔFin, sysstatus = update_state!(model.elems, ΔUit, Δt)
            failed(sysstatus) && (syserror=true; break)

            Vt = -V + 2/Δt*ΔUit
            At = -A + 4/Δt^2*(ΔUit - Vt*Δt)
            
            R = ΔFexi - (ΔFin + M*At)
            res = maximum(abs, R[umap] )

            # Update accumulated displacement
            ΔUi .+= ΔUk

            # Residual vector for next iteration
            # ΔF[pmap] .= 0.0  # Zero at prescribed positions

            @printf(env.log, "    it %d  res: %-10.4e\n", it, res)

            # @show res

            it==1 && (res1=res)
            res < tol  && (converged=true; break)
            isnan(res) && break
            it>maxits  && break
            it>1 && res>lastres && break
        end

        q = 0.0 # increment size factor for autoinc

        if syserror
            println(env.log, sysstatus.message)
            converged = false
        end
        quiet || sysstatus.message!="" && println(env.alerts, sysstatus.message)

        if converged
            # Update nodal natural and essential values for the current stage
            U  .+= ΔUi
            F  .+= ΔFin

            Uexi! = Uexi
            Fexi! = Fexi

            # @showm ΔUi

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                # @show dof.name, ΔUi[i]
                dof.vals[dof.name]    += ΔUi[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            # Update time
            t += Δt
            T += ΔT
            env.t = t
            env.T = T
            env.residue = res

            # Check for saving output file
            if T>Tcheck-Tmin && save_outs
                env.out += 1
                iout = env.out
                rm.(glob("*conflicted*.dat", "$outdir/"), force=true)
                
                update_output_data!(model)

                update_multiloggers!(model)
                # @show "$outdir/$outkey-$iout.vtu"
                save(model, "$outdir/$outkey-$iout.vtu", quiet=true)

                Tcheck += ΔTcheck # find the next output time
            end

            update_single_loggers!(model)
            update_monitors!(model)
            flush(env.log)

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        q = 1.0 + tanh(log10(tol/res1))
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
                    if T+ΔTtr>Tcheck-Tmin
                        ΔTbk = ΔT
                        ΔT = Tcheck-T
                        @assert ΔT>=0.0
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            # Restore counters
            inc -= 1
            env.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(env.log, "      increment failed")
                q = (1+tanh(log10(tol/res1)))
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < Tmin
                    solstatus = failure("solver did not converge")
                    break
                end
            else
                solstatus = failure("solver did not converge")
                break
            end
        end
    end

    if !save_outs
        update_output_data!(model)
        complete_ut_T(model)
        update_multiloggers!(model)
    end

    return solstatus

end
