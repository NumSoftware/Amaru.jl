# From the book:
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru
using Test

# Mesh generation

block = Block( [0 0; 1 1], nx=1, ny=1, cellshape=QUAD4, tag="solid")
mesh = Mesh(block, quiet=true)

# Model definition

mats = [
    "solid" => MechSolid => LinearElastic => (E=1.0, nu=0.25),
]

ctx = MechContext(stressmodel=:planestrain)
model = FEModel(mesh, mats, ctx)
ana = MechAnalysis(model)

bcs = [
    :(x==0.) => SurfaceBC(ux=0.),
    :(y==0.) => SurfaceBC(uy=0),
    :(y==1.) => SurfaceBC(ty=-1.),
]
addstage!(ana, bcs, nincs=1)
solve!(ana)


dis = 
    [
        0.0       0.0000e+00 
        0.3125    0.0000e+00 
        0.0      -9.3750e-01 
        0.3125   -9.3750e-01 
    ]
           

println("Displacements:")

D = get_data(model.nodes)[[:ux, :uy]]
println(D)

@test dis ≈ Array(D) atol=1e-5

println("Stress:")
S = elems_ip_vals(model.elems[1])[[:σxx, :σyy, :σxy]]
println(S)

println("Support reactions:")
F = get_data(model.nodes)[[:fx, :fy]]
println(F)