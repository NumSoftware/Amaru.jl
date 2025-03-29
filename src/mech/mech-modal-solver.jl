# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MechModalAnalysis

MechModalAnalysis_params = [
    FunInfo(:DynAnalysis, "Dynamic analysis"),
    ArgInfo(:model, "Finite element model"),
]
@doc docstring(MechModalAnalysis_params) MechModalAnalysis


mutable struct MechModalAnalysis<:StaticAnalysis
    model ::FEModel
    ctx   ::MechContext
    sctx  ::SolverContext

    stages  ::Array{Stage}
    loggers ::Array{AbstractLogger,1}
    monitors::Array{AbstractMonitor,1}

    freqs::Array{Float64,1} # frequencies
    modes::Array{Float64,2} # modes

    function MechModalAnalysis(model::FEModel; outdir=".", outkey="out")
        this = new(model, model.ctx)
        this.stages = []
        this.loggers = []
        this.monitors = []  

        this.freqs = zeros(0)
        this.modes = zeros(0,0)

        this.sctx = SolverContext()
        this.sctx.outdir = outdir
        this.sctx.outkey = outkey

        model.ctx.thickness = model.thickness
        if model.ctx.stressmodel==:none
            if model.ctx.ndim==2
                model.ctx.stressmodel = :planestrain
            else
                model.ctx.stressmodel = :d3
            end
        end

        return this
    end
end


mech_modal_solver_params = [
    FunInfo( :mech_modal_solver!, "Finds the frequencies and vibration modes of a mechanical system."),
    ArgInfo( :model, "FEModel object"),
    ArgInfo( :stage, "Stage object"),
    KwArgInfo( :nmodes, "Number of modes to be calculated", 5),
    KwArgInfo( :rayleigh, "Flag to use Rayleigh-Ritz method for damping", false),
    KwArgInfo( :quiet, "Flag to set silent mode", false),
]
@doc docstring(mech_modal_solver_params) solve!(::MechModalAnalysis; args...)

function solve!(ana::MechModalAnalysis; args...)
    args = checkargs(args, mech_modal_solver_params)
    if !args.quiet
        printstyled("Solver for mechanical modal analyses", "\n", bold=true, color=:cyan)
        println("  stress model: ", ana.ctx.stressmodel)
    end

    status = stage_iterator!(mech_modal_solver!, ana; args...)
    return status
end


function mech_modal_solver!(ana::MechModalAnalysis, stage::Stage; kwargs...)
    args = NamedTuple(kwargs)
    # args = checkargs(kwargs, mech_modal_solver_params)
    nmodes = args.nmodes
    rayleigh = args.rayleigh
    quiet = args.quiet 

    model = ana.model
    ctx = model.ctx
    sctx = ana.sctx

    quiet || println(sctx.log, "Modal analysis for mechanical systems")

    stressmodel = ctx.stressmodel
    ctx.ndim==3 && @check stressmodel==:d3

    # todo: check there are not force boundary conditions

    # get only bulk elements
    model = FEModel(model.elems.bulks)

    # check density
    for elem in model.elems
        elem.props.ρ==0 && error("mech_modal_solver: density should not be zero")
    end

    # get dofs organized according to boundary conditions
    dofs, nu    = configure_dofs!(model, stage.bcs)

    ndofs       = length(dofs)
    model.ndofs = length(dofs)

    if !quiet
        println(sctx.log, "unknown dofs: $nu")
        println(sctx.info, "unknown dofs: $nu")
    end

    # setup quantities at dofs
    for dof in dofs
        dof.vals[dof.name]    = 0.0
        dof.vals[dof.natname] = 0.0
    end

    K11 = mount_K(model.elems, ndofs)[1:nu, 1:nu]
    M11 = mount_M(model.elems, ndofs)[1:nu, 1:nu]

    M11 = 0.5*(M11 + M11')

    L = sum(M11, dims=1)

    # m11 = zeros(nu) #Vetor of inverse matrix mass
    # for i in 1:size(M11,2)
    #     m11[i] = 1.0/M11[i,i]
    # end

    # inverse of the lumped matrix in vector form
    # P = m11.*K11
    P = K11 ./ L

    # @show K11

    # m11 = inv(Matrix(M11))
    # P = m11*K11

    # eingenvalues and eingenvectors
    # Eig = eigs(P, nev=nmodes, which=:SR, tol=1e-6, maxiter=600, check=1) # SR: smallest real part
    nenv = size(P,2) - 2
    Eig = eigs(P, nev=nmodes, which=:SM, tol=1e-6, ncv=min(2*nmodes, nenv), check=1) # SM: smallest magnitude
    # Eig = eigs(P, nev=nmodes, which=:LM) # LM: largest magnitude
    w0  = Eig[1] # frequencies
    wi  = copy(w0)

    @show wi
    @show wi
    @show wi

    # select possible vals
    filter = [ i for i in eachindex(wi) if isreal(wi[i]) && real(wi[i])>0 ]
    filter = filter[ unique(i -> wi[i], filter) ]  # todo: do not remove eingenvalues
    perm   = sortperm(real(wi[filter]))[1:nmodes]
    filter = filter[perm]

    wi = wi[filter] # sorted
    v  = Eig[2][:, filter]

    w = wi.^0.5  # true frequencies

    update_output_data!(model)
    save(model, joinpath(sctx.outdir, "$(sctx.outkey)-0.vtu"), quiet=true)

    sctx.inc = 1
    sctx.ΔT  = 1.0
    
    # save modes
    for i in 1:nmodes
        U = v[:,i] # modal displacements
        
        for (k,dof) in enumerate(dofs[1:nu])
            dof.vals[dof.name] = U[k]
        end
        sctx.T = i/nmodes
        sctx.out = i
        update_output_data!(model)
        save(model, joinpath(sctx.outdir, "$(sctx.outkey)-$i.vtu"), quiet=true)
    end

    # reset displacement values
    for dof in dofs
        dof.vals[dof.name] = 0.0
    end

    # show modal frequencies
    # if !quiet
    #     println(sctx.log, "modal frequencies:")
    #     for i in 1:nmodes
    #         println(sctx.log, "ω$i = ", abs(w[i]))
    #     end
    # end

    # show("Modal Frequencies rad/s")
    # @show w
    #show("Modal Shapes rad/s")
    #@show v

    if rayleigh
        #Rayleigh-Ritz method for Damping

        println("What is the damping relationship 1 (xi1 = c/c_c)?")
        xi1 = parse(Float64,chomp(readline()))
        println("What is the damping relationship 2 (xi2 = c/c_c)?")
        xi2 = parse(Float64,chomp(readline()))


        #Interacting with the user for select the correct frequêncies based on modal deformated

        println("See the vtk files of the modal deformated and decide")
        print("What is the first frequency that you want?")
        w1 = w[parse(Int,chomp(readline()))]
        print("What is the second frequency that you want?")
        w2 = w[parse(Int,chomp(readline()))]

        alpha = 2*(w1^2*w2*xi2 - w1*w2^2*xi1) / (w1^2-w2^2)
        beta  = 2*(w1*xi1 - w2*xi2)/(w1^2 - w2^2)

        info("Rayleigh alpha: ", alpha)
        info("Rayleigh beta:", beta)
    end

    # save frequencies and modes
    ana.freqs = w
    ana.modes = v

    return success()

end

