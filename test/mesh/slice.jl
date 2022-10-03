using Amaru

bl = Block([0 0 0; 1 1 1], nx=15, ny=15, nz=15, cellshape=HEX8)
mesh = Mesh(bl)

mat = [ :bulks => ElasticSolid(E=100, nu=0.2) ]
model = Model(mesh, mat)

bcs = [
    :(z==0) => NodeBC(ux=0, uy=0, uz=0)
    :(z==1) => SurfaceBC(tz=-1)
]

addstage!(model, bcs)

solve!(model)

m = slice(model, base=[0.5,0.55,0.55], axis=[1,1,1.5])
nothing
# @show m.elems