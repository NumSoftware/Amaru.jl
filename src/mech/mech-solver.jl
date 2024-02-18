# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechAnalysis

mutable struct MechAnalysisProps<:Analysis
    stressmodel::String # plane stress, plane strain, etc.
    thickness::Float64  # thickness for 2d analyses
    g::Float64 # gravity acceleration
    
    function MechAnalysisProps(;stressmodel="3d", thickness=1.0, g=0.0)
        @check stressmodel in ("plane-stress", "plane-strain", "axisymmetric", "3d")
        @check thickness>0
        @check g>=0
        return new(stressmodel, thickness, g)
    end
end

MechAnalysis = MechAnalysisProps


# Assemble the global stiffness matrix
function mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            Ke, rmap, cmap = elem_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i,j]
                    abs(val) < eps() && continue
                    # @redude R = append!(Int64[], rmap[i])
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    # nelems = length(elems)
    # pool = Array{Any}(undef, nelems)

    # Ri = Array{Array{Int,1}}(undef, nelems)
    # Ci = Array{Array{Int,1}}(undef, nelems)
    # Vi = Array{Array{Float64,1}}(undef, nelems)

    # # Threads.nthreads()>1 && LinearAlgebra.BLAS.set_num_threads(1)


    # @time begin
    #     # ccall(:jl_enter_threaded_region, Cvoid, ())
    #     @sync for elem in elems
    #     # Threads.@threads for elem in elems


    #         Threads.@spawn begin
    #             Ke, rmap, cmap = elem_stiffness(elem)

    #             I, J, V = findnz(sparse(Ke))
    #             Ri[elem.id] = rmap[I]
    #             Ci[elem.id] = cmap[J]
    #             Vi[elem.id] = V

    #         end

            
    #     end
    #     # ccall(:jl_exit_threaded_region, Cvoid, ())

    #     println("100")

            
    #     R = reduce(vcat, Ri)
    #     C = reduce(vcat, Ci)
    #     V = reduce(vcat, Vi)
    #     # Threads.nthreads()==1 && LinearAlgebra.BLAS.set_num_threads(8)
        
    #     # R = reduce(vcat, p[1] for p in pool)
    #     # println("200")
    #     # C = reduce(vcat, p[2] for p in pool)
    #     # println("200")
    #     # V = reduce(vcat, p[3] for p in pool)
    #     # R, C, V = Int64[], Int64[], Float64[]



        
    # end
    
    
    # println("xxxxxxx")
    
    # error()
    
    local K
    try
        K = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show err
    end

    yield()

    return K
end


function mount_RHS(model::Model, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in active_elems
        F, map = elem_RHS(elem) # RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = F
    end
    return RHS
end


# Get internal forces and update data at integration points
function update_state!(active_elems::Array{<:Element,1}, ΔUt::Vect, t::Float64)
    ndofs = length(ΔUt)
    @withthreads begin 
        ΔFin     = zeros(ndofs)
        statuses = ReturnStatus[]
        for elem in active_elems
            ΔF, map, status = update_elem!(elem, ΔUt, t)
            if failed(status) 
                push!(statuses, status)
                break
            end
            ΔFin[map] .+= ΔF
        end
    end

    yield()
    
    length(statuses)>0 && return ΔFin, statuses[1]
    any(isnan.(ΔFin)) && return ΔFin, failure("solve_system!: NaN values in internal forces vector")
    return ΔFin, success()
end


function update_embedded_disps!(active_elems::Array{<:Element,1}, U::Matx)
    for elem in active_elems.embeddeds
        Ue, nodemap, dimmap = elem_displacements(elem)
        U[nodemap, dimmap] .= Ue
    end
end


"""
    solve!(model,options...) :: Bool

Performs one stage static finite element analysis of a domain `model`
subjected to a set of boundary conditions `bcs`.

# Arguments

`model` : A finite element domain

`bcs` : Array of boundary conditions given as an array of pairs ( location => condition)

# Keyword arguments

`nincs   = 1` : Number of increments

`maxits  = 5` : Maximum number of Newton-Rapson iterations per increment

`autoinc = false` : Sets automatic increments size. The first increment size will be `1/nincs`

`rtol     = 1e-2` : Tolerance for the maximum absolute error in forces vector

`ΔTmin     = 1e-8` : Pseudo-time tolerance

`scheme  = "FE"` : Predictor-corrector scheme at each increment. Available schemes are "FE", "ME", "BE", "Ralston"

`nouts   = 0` : Number of output files per analysis

`outdir  = ""` : Output directory

`outkey = ""` : File key for output files

`quiet = false` : verbosity level from 0 (silent) to 2 (verbose)

"""
function solve!(model::Model, ana::MechAnalysis; args...)
    name = "Solver for mechanical analyses"
    status = stage_iterator!(name, mech_stage_solver!, model; args...)
    return status
end


mech_stage_solver_params = [
    FunInfo( :mech_stage_solver!, "Solves a load stage of a mechanical analysis.", "M::Model, S::Stage"),
    ArgInfo( :tol, "Force tolerance", 0.01, condition=:(tol>0)),
    ArgInfo( :rtol, "Relative tolerance in terms of displacements", 0.01, condition=:(rtol>0)),
    ArgInfo( :dTmin, "Relative minimum increment size", 1e-7, condition=:(0<dTmin<1) ),
    ArgInfo( :dTmax, "Relative maximum increment size", 0.1, condition=:(0<dTmax<1) ),
    ArgInfo( :rspan, "Relative span to residue reapplication", 0.01, condition=:(0<rspan<1) ),
    ArgInfo( :scheme, "Global solving scheme", :FE, values=(:FE, :ME, :BE, :Ralston) ),
    ArgInfo( :maxits, "Maximum number of NR iterations", 5, condition=:(1<=maxits<=10)),
    ArgInfo( :autoinc, "Flag to set auto-increments", false),
    ArgInfo( :quiet, "Flat to set silent mode", false),
]
@doc make_doc(mech_stage_solver_params) mech_stage_solver!()

function mech_stage_solver!(model::Model, stage::Stage; args...)
    args = checkargs(args, mech_stage_solver_params)
    
    tol     = args.tol      
    rtol    = args.rtol     
    ΔTmin   = args.dTmin    
    ΔTmax   = args.dTmax   
    rspan   = args.rspan    
    scheme  = args.scheme   
    maxits  = args.maxits 
    autoinc = args.autoinc  
    quiet   = args.quiet    

    env = model.env
    println(env.log, "Mechanical FE analysis: Stage $(stage.id)")

    solstatus = success()
    # scheme in ("FE", "ME", "BE", "Ralston") || error("solve! : invalid scheme \"$(scheme)\"")

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    env       = model.env
    saveouts = stage.nouts > 0
    ftol      = tol

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, stage.bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    model.ndofs = length(dofs)

    println(env.info,"unknown dofs: $nu")
    println(env.log, "unknown dofs: $nu")

    quiet || nu==ndofs && println(env.alerts, "No essential boundary conditions")

    if stage.id == 1
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end

        update_records!(model)
    end
    lastflush = time()
    flushinterval = 5.0

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, ΔTmax))

    ΔTbk    = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    inc  = 0             # increment counter
    iout = env.out       # file output counter
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    still_linear = true
    sysstatus = ReturnStatus()

    # Get boundary conditions
    Uex, Fex = get_bc_vals(model, bcs)

    # Get unbalanced forces from activated elements
    if length(stage.toactivate)>0
        Fin = zeros(ndofs)
        for elem in stage.toactivate
            elem_internal_forces(elem, Fin)
        end
        Fex .-= Fin # add negative forces to external forces vector
    end

    if scheme==:FE
        p1=1.0; q11=1.0
    elseif scheme==:ME
        p1=1.0; q11=1.0; a1=0.5; a2=0.5
    elseif scheme==:BE
        p1=1.0; q11=1.0; a1=0.0; a2=1.0
    elseif scheme==:Ralston
        p1=2/3; q11=2/3; a1=1/4; a2=3/4
    end

    local K::SparseMatrixCSC{Float64,Int64}

    while T < 1.0-ΔTmin
        env.ΔT = ΔT

        # Update counters
        inc += 1
        env.inc = inc

        println(env.log, "  inc $inc")

        ΔUex, ΔFex = ΔT*Uex, ΔT*Fex     # increment of external vectors

        ΔTcr = min(rspan, 1-T)    # time span to apply cumulated residues
        αcr  = min(ΔT/ΔTcr, 1.0)  # fraction of cumulated residues to apply
        T<1-rspan && (ΔFex .+= αcr.*Rc) # addition of residuals

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex  # essential values at iteration i

        # Newton Rapshon iterations
        nits      = 0
        err       = 0.0
        res       = 0.0
        res1      = 0.0
        converged = false
        syserror  = false

        for it=1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = err # residue from last iteration

            # Predictor step
            K = mount_K(active_elems, ndofs)
            ΔUitr = p1*ΔUi
            Rtr   = q11*R
            sysstatus = solve_system!(K, ΔUitr, Rtr, nu)   # Changes unknown positions in ΔUi and R
            failed(sysstatus) && (syserror=true; break)
            copyto!.(State, StateBk)
            ΔUt   = ΔUa + ΔUitr
            ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
            failed(sysstatus) && (syserror=true; break)

            # Corrector step
            if scheme==:FE
                ΔUi = ΔUitr
            else
                K2 = mount_K(active_elems, ndofs)
                K = a1*K + a2*K2
                sysstatus = solve_system!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (syserror=true; break)
                copyto!.(State, StateBk)
                ΔUt   = ΔUa + ΔUi
                ΔFin, sysstatus = update_state!(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (syserror=true; break)
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions
            res = norm(R, Inf)

            err = norm(ΔUi, Inf)/norm(ΔUa, Inf)

            @printf(env.log, "    it %d  residue: %-10.4e\n", it, res)

            it==1 && (res1=res)
            it>1  && (still_linear=false)
            res<ftol && (converged=true; break)
            err<rtol && (converged=true; break)

            isnan(res) && break
            it>maxits  && break

            it>1 && res>lastres && break
        end

        env.residue = res
        q = 0.0 # increment size factor for autoinc

        if syserror
            println(env.alerts, sysstatus.message)
            converged = false
        end

        if converged
            # Update forces and displacement for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Rc .= (1.0-αcr).*Rc .+ R  # update cumulated residue

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            # # Update nodal displacements
            # for node in model.nodes
            #     ux = ΔUa[node.dofdict[:ux].eq_id] 
            #     uy = ΔUa[node.dofdict[:uy].eq_id] 
            #     uz = ΔUa[node.dofdict[:uz].eq_id] 
            #     node.coord = node.coord + [ux, uy, uz]
            # end

            # Update time
            T += ΔT
            env.T = T

            # Check for saving output file
            checkpoint = T>Tcheck-ΔTmin
            # && saveouts
            if checkpoint
                env.out += 1
                update_embedded_disps!(active_elems, model.node_data["U"])
                Tcheck += ΔTcheck # find the next output time
            end

            flush = time()-lastflush>flushinterval || T >= 1.0-ΔTmin

            rstatus = update_records!(model, checkpoint=checkpoint, flush=flush)
            if failed(rstatus)
                println(env.alerts, rstatus.message)
                return rstatus
            end
            
            if flush
                lastflush = time()
                GC.gc()
            end

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    # if nits==1
                    #     # q = 1+tanh(log10(ftol/res))
                    #     q = 1.333
                    # else
                    #     q = 1+tanh(log10(rtol/err))
                    # end
                    q = 1+tanh(log10(ftol/(res1+eps())))
                    q = max(q, 1.1)

                    if still_linear
                        ΔTtr = min(q*ΔT, ΔTmax, 1-T)
                    else
                        ΔTtr = min(q*ΔT, ΔTmax, 1-T)
                    end

                    if T+ΔTtr>Tcheck-ΔTmin
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
                # if nits==1
                #     # q = 1+tanh(log10(ftol/res))
                #     q = 0.666
                # else
                #     q = 1+tanh(log10(rtol/err))
                # end

                q = 1+tanh(log10(ftol/(res1+eps())))
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < ΔTmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Solver did not converge. Try `autoinc=true`. ")
                break
            end
        end
    end

    failed(solstatus) && update_records!(model, solverfailed=true)

    return solstatus

end
