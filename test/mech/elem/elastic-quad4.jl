# From the book:
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru
using Test

# Mesh generation

block = Block( [0 0; 1 1], nx=1, ny=1, cellshape=QUAD4, tag="solid")
mesh = Mesh(block, verbosity=0, reorder=false)

# Domain definition

mats = [
    "solid" => ElasticSolid(E=1.0, nu=0.25),
]

dom = Domain(mesh, mats)

bcs = [
    :(x==0.) => FaceBC(ux=0.),
    :(y==0.) => FaceBC(uy=0),
    :(y==1.) => FaceBC(ty=-1.),
]

solve!(dom, bcs, verbosity=1)


dis = 
    [
        0.0       0.0000e+00 
        0.3125    0.0000e+00 
        0.0      -9.3750e-01 
        0.3125   -9.3750e-01 
    ]
           

println("Displacements:")
D = get_data(dom.nodes)[[:ux, :uy]]
println(D)

@test dis ≈ Array(D) atol=1e-5

println("Stress:")
S = elems_ip_vals(dom.elems[1])[[:sxx, :syy, :sxy]]
println(S)

println("Support reactions:")
F = get_data(dom.nodes)[[:fx, :fy]]
println(F)