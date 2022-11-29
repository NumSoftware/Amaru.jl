export mod_solve!
#Solids stiffness matrix for computing frequencies

# function mount_Kf(model::Model, ndofs::Int)

#     R, C, V = Int64[], Int64[], Float64[]
#     elems=model.elems.solids

#     for elem in elems
#        Ke, rmap, cmap = elem_stiffness(elem)
#         nr, nc = size(Ke)
#         for i in 1:nr
#             for j in 1:nc
#                 push!(R, rmap[i])
#                 push!(C, cmap[j])
#                 push!(V, Ke[i,j])
#             end
#         end
#     end


#     return sparse(R, C, V, ndofs, ndofs)

# end


# #Solids mass matrix for computing frequencies

# function mount_Mf(model::Model, ndofs::Int)

#     R, C, V = Int64[], Int64[], Float64[]
#     elems=model.elems.solids

#     for elem in elems
#         Me, rmap, cmap = elem_mass(elem)
#         nr, nc = size(Me)
#         for i in 1:nr
#             for j in 1:nc
#                 push!(R, rmap[i])
#                 push!(C, cmap[j])
#                 push!(V, Me[i,j])
#             end
#         end
#     end

#     return sparse(R, C, V, ndofs, ndofs)
# end


function mod_solve!(model::Model, bcs::Array; nmods::Int=5, rayleigh=false, outdir="", quiet=false)

    quiet || printstyled("FEM modal analysis:\n", bold=true, color=:cyan)

    # get only bulk elements
    model = Model(model.elems.bulks)

    # get dofs organized according to boundary conditions
    dofs, nu    = configure_dofs!(model, bcs)

    ndofs       = length(dofs)
    model.ndofs = length(dofs)
    quiet || println("  unknown dofs: $nu")

    # setup quantities at dofs
    for dof in dofs
        dof.vals[dof.name]    = 0.0
        dof.vals[dof.natname] = 0.0
    end

    K11 = mount_K(model.elems, ndofs)[1:nu, 1:nu]
    M11 = mount_M(model.elems, ndofs)[1:nu, 1:nu]

    m11 = zeros(nu) #Vetor of inverse matrix mass
    for i in 1:size(M11,2)
        m11[i] = 1.0/M11[i,i]
    end

    # inverse of the lumped matrix in vector form
    P = m11.*K11

    # eingenvalues and eingenvectors
    Eig = eigs(P, nev=20, which=:SM)
    w0  = Eig[1] # frequencies
    wi  = copy(w0)

    # select possible vals
    filter = [ i for i in eachindex(wi) if isreal(wi[i]) && real(wi[i])>0 ]
    filter = filter[ unique(i -> wi[i], filter) ]
    perm   = sortperm(real(wi[filter]))[1:nmods]
    filter = filter[perm]

    wi     = wi[filter] # sorted
    v      = Eig[2][:, filter]

    w = wi.^0.5

    update_output_data!(model)
    save(model, joinpath(outdir, "mod-0.vtu"))

    # save modes
    for i in 1:nmods
        U = v[:,i] # modal displacements

        for (k,dof) in enumerate(dofs[1:nu])
            dof.vals[dof.name] = U[k]
        end
        update_output_data!(model)
        save(model, joinpath(outdir, "mod-$i.vtu"))
    end

    # reset displacement values
    for dof in dofs
        dof.vals[dof.name] = 0.0
    end

    # show modal frequencies
    if !quiet
        message("modal frequencies:")
        for i in 1:nmods
            # @show i
            # println("ω$i = ", abs(w[i]))
            info("ω$i = ", abs(w[i]))
        end
    end

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

end

function modsolvex!(model::Model, bcs::Array; nmods::Int=5, rayleigh=false, save=true, quiet=false)

    quiet || printstyled("FEM modal analysis:\n", bold=true, color=:cyan)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs!(model, bcs)

    ndofs = length(dofs)
    model.ndofs = length(dofs)
    !quiet && println("  unknown dofs: $nu")

    K = mount_Kf(model, ndofs)
    M = mount_Mf(model, ndofs)

    if nu>0
        K11 = K[1:nu, 1:nu]
        M11 = M[1:nu, 1:nu]
    end

    #Delete columns and rows that corresponding at bars and joint elements

    nodesbars=getnodes(model.elems.lines) #nodes of bars

    idsd = Array{Int64,1}() #idsd ids for delete

    for node in nodesbars
        for dof in node.dofs
            push!(idsd,dof.eq_id)
        end
    end
    idsd = unique(idsd)
    sort!(idsd)
    nidsd = length(idsd)

    for i in 1:nidsd
       M11=M11[1:size(M11,1)     .!= idsd[nidsd+1-i,1],:];
       M11=M11[:,1:size(M11,2) .!= idsd[nidsd+1-i,1]];
       K11=K11[1:size(K11,1)     .!= idsd[nidsd+1-i,1],:];
       K11=K11[:,1:size(K11,2) .!= idsd[nidsd+1-i,1]];
    end

    m11 = zeros(size(M11,2)) #Vetor of inverse matrix mass

    for i in 1:size(M11,2)
        m11[i] = 1/M11[i,i]
    end

    # Inverse of the lumped matrix in vector form
    #m11 = inv.(sum(M11, 2))
    P = m11.*K11


    #Eingenvalues and Eingenvectors

    nevn=size(P,2)-2 #eigvals number

    #using Arpack

    Eig=eigs(P, nev=20, which=:SM)

    #using LinearAlgebra
    w0 = Eig[1]
    #w0 = eigvals(P)
    wi = copy(w0)

    #Select logical and possible vals

    wi=real(wi[isreal.(wi)]) #delete imaginarian numbers
    wi=unique(wi) #all vals diferents
    wi=wi[wi[1:size(wi)[1]] .>0] # delete zero and negative vals
    wi=sort!(wi) # sort the vector

    w = Array{Float64,1}()

    for i in 1:nmods
        push!(w,wi[i])
    end

    v = Array{Float64,2}(undef,size(P,2),nmods)

    for i in 1:nmods
        col = findfirst(isequal(w[i]),w0)
        col = col[1]
        #v[:,i] = eigvecs(P)[:,col]
        v[:,i] = Eig[2][:,col]
    end

    w = w.^0.5

    if save

        #Modals Deformed

        nodessolids = getnodes(model.elems.solids) #nodes of solids

        idss = Array{Int64,1}() #idss ids of dof solids

        for node in nodessolids
            for dof in node.dofs
                push!(idss,dof.eq_id)
            end
        end
        idss = unique(idss)
        sort!(idss)
        nidss = length(idss)
        idefmod = 1

        for j in 1:nmods

            Umod = zeros(ndofs)

            for i in 1:length(v[:,j])
                Umod[idss[i]] = v[i,j]
            end

            for (k,dof) in enumerate(dofs)
                dof.vals[dof.name] = Umod[k]
            end

            save(model, model.filekey * "-defmod$idefmod.vtk") # saves current domain state for modal deformated
            idefmod += 1
        end

    end


    # reset displacement values
    for (k,dof) in enumerate(dofs)
        dof.vals[dof.name] = 0.0
    end

    show("Modal Frequencies rad/s")

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

        alpha = 2 * ((w1 ^ 2) * w2 * xi2 - w1 * (w2 ^ 2) * xi1) / ((w1 ^ 2)-(w2 ^ 2))
        beta  = 2 * (w1 * xi1 - w2 * xi2)/((w1 ^ 2) - (w2 ^ 2))

        info("Rayleigh alpha: ", alpha)
        info("Rayleigh beta:", beta)

    end


end

