# From the book:
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Amaru
using Test

# Mesh generation

block = Block( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, cellshape=HEX8, tag="solid")
mesh = Mesh(block, quiet=true)

# Model definition

materials = [
    "solid" => MechSolid => LinearElastic => (E=1.0, nu=0.3),
]

# Load cases

bcs1 = [
    :(x==0 && y==0 && z==0) => NodeBC(ux=0, uy=0), 
    :(x==1 && y==0 && z==0) => NodeBC(uy=0),
    :(x==0 && y==1 && z==0) => NodeBC(ux=0),
    :(z==0)                 => NodeBC(uz=0),
    :(z==1)                 => NodeBC(fz=1),
]

bcs2 = [
    :(x==0 && y==0 && z==0) => NodeBC(ux=0, uy=0),
    :(x==1 && y==0 && z==0) => NodeBC(uy=0),
    :(x==0 && y==1 && z==0) => NodeBC(ux=0),
    :(z==0)                 => NodeBC(uz=0),
    :(y==1 && z==1)         => EdgeBC(qy=2),
]

bcs3 = [
    :(x==0 && y==0 && z==0) => NodeBC(ux=0, uy=0),
    :(x==1 && y==0 && z==0) => NodeBC(uy=0),
    :(x==0 && y==1 && z==0) => NodeBC(ux=0),
    :(z==0)                 => NodeBC(uz=0),
    :(x==1)                 => SurfaceBC(tx=:(3*z)),
]

bcs4 = [
    :(x==0 && y==0 && z==0) => NodeBC(ux=0, uy=0),
    :(x==1 && y==0 && z==0) => NodeBC(uy=0),
    :(x==0 && y==1 && z==0) => NodeBC(ux=0),
    :(z==0)                 => NodeBC(uz=0),
    :(x>=0)                 => BodyC(wz=-1),
]

lc_list = ["Nodal load", "Edge load", "Triangular face load", "Volume load"]
bcs_list = [bcs1, bcs2, bcs3, bcs4]
dis_list = [
            [ 0.0, 0.0, 0.0, 0.0, 4.0    ,    4.0,     4.0,      4.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 3.32088, 3.1998, -4.6002, -4.32047 ],
            [ 0.0, 0.0, 0.0, 0.0, 1.51044,-2.4501,  1.4499, -2.31023 ],
            [ 0.0, 0.0, 0.0, 0.0, -0.5   ,   -0.5,    -0.5,     -0.5 ] ]

for (lc, bcs, dis) in zip(lc_list, bcs_list, dis_list)

    println("\nLoad case: $lc \n")

    ctx = MechContext()
    model = FEModel(mesh, materials, ctx)
    ana = MechAnalysis(model)
    addstage!(ana, bcs, nouts=1)
    solve!(ana)

    println("Displacements:")
    D = get_data(model.nodes)[[:ux, :uy, :uz]]
    println(D)

    @test dis ≈ D[:uz] atol=1e-5

    println("Stress:")
    S = elems_ip_vals(model.elems[1])[[:σxx, :σyy, :σzz, :σyz, :σxz, :σxy]]
    println(S)

    println("Support reactions:")
    F = get_data(model.nodes[:(z==0)])[[:fx, :fy, :fz]]
    println(F)
end


