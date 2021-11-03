using Amaru, LinearAlgebra
using Test

printstyled("\nShape functions\n", color=:blue, bold=true)

for shape in ALL_ISO_SHAPES
    print("shape : ", shape.name)
    n  = shape.npoints
    ndim = shape.ndim

    In = Matrix{Float64}(I,n,n)

    # Check at nodes
    RR = [ shape.nat_coords[i,:] for i=1:n ]
    NN = shape.func.(RR)

    II = hcat(NN...) # should provida an identity matrix
    @test II ≈ In atol=1e-10

    # Check at default set of integration points
    Q  = shape.quadrature[0]
    nip, _ = size(Q)

    RR = [ Q[i,:] for i=1:nip ]
    NN = shape.func.(RR)
    TR = @test sum(sum(NN)) ≈ nip atol=1e-10
    println(TR)
end

printstyled("\nShape functions derivatives\n", color=:blue, bold=true)

for shape in ALL_ISO_SHAPES
    print("shape : ", shape.name)
    n    = shape.npoints
    ndim = shape.ndim
    RR = [ shape.nat_coords[i,:] for i=1:n ]
    f  = shape.func

    In = Matrix{Float64}(I,n,n)
    Id = Matrix{Float64}(I,ndim,ndim)

    # numerical derivative
    δ  = 1e-8
    shape_pass = true
    TR = Test.Pass(:test, nothing, nothing, true)
    for R in RR
        RI = R .+ Id*δ
        fR = f(R)
        D  = zeros(n, ndim)
        for i=1:ndim
            Di     = 1/δ*(f(RI[:,i]) - fR)
            D[:, i] = Di
        end

        tr = @test D ≈ shape.deriv(R) atol=1e-6
        tr isa Test.Pass || (TR = tr)
    end
    println(TR)
end
