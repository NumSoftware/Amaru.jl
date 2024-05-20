using Amaru

tr, J2, J3, eigen, dev, I2 = Amaru.tr, Amaru.J2, Amaru.J3, Amaru.eigen, Amaru.dev, Amaru.I2

σ = Amaru.Vec6(10,20,30,4,5,6)
s = dev(σ)

println("Eigenvalues of the stress tensor")
A, V = eigen(σ)
@test sum(A) ≈ tr(σ) atol=1e-10

println("Invariants of the deviatoric tensor")
@test s ≈ σ - 1/3*tr(σ)*I2
@test dev(s) ≈ s atol=1e-10

@test tr(σ) ≈ 60 atol=1e-10
@test tr(s) ≈ 0 atol=1e-10
@test J2(σ) ≈ J2(s) atol=1e-10
@test J3(σ) ≈ J3(s) atol=1e-10
@test dot(σ,inv(σ)) ≈ 3 atol=1e-10


# @test J3(σ) ≈ 1/J3(inv(s)) atol=1e-10
# @show σ*inv(σ)'
# @show s*inv(s)'