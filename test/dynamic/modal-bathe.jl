using Amaru

# Entrada de dados
# Concreto
fc = 25800000/1000                    # resistência à compressão
Ec  = (6100/0.14503773800722)*1000
Ea = (30000/0.14503773800722)*1000
Aa = 2*(0.0254*0.0254)
p  = √(4*Aa*π)
roc = 2626.96/1000#((0.217e-6)/(1.1228705601979e-6))
roa = 7850/1000


# Mesh generation


bl1  = Block( (1/1000)*[0   0      0;6*25.4 20*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl2  = Block( (1/1000)*[0 20*25.4  0;6*25.4 38*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl3  = Block( (1/1000)*[0 38*25.4  0;6*25.4 50*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl4  = Block( (1/1000)*[0 50*25.4  0;6*25.4 59*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl5  = Block( (1/1000)*[0 59*25.4  0;6*25.4 68*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl6  = Block( (1/1000)*[0 68*25.4  0;6*25.4 77*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl7  = Block( (1/1000)*[0 77*25.4  0;6*25.4 86*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl8  = Block( (1/1000)*[0 86*25.4  0;6*25.4 98*25.4   3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl9  = Block( (1/1000)*[0 98*25.4  0;6*25.4 116*25.4  3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl10 = Block( (1/1000)*[0 116*25.4 0;6*25.4 136*25.4  3*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)

bl11  = Block( (1/1000)*[0   0      3*25.4;6*25.4 20*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl12  = Block( (1/1000)*[0 20*25.4  3*25.4;6*25.4 38*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl13  = Block( (1/1000)*[0 38*25.4  3*25.4;6*25.4 50*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl14  = Block( (1/1000)*[0 50*25.4  3*25.4;6*25.4 59*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl15  = Block( (1/1000)*[0 59*25.4  3*25.4;6*25.4 68*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl16  = Block( (1/1000)*[0 68*25.4  3*25.4;6*25.4 77*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl17  = Block( (1/1000)*[0 77*25.4  3*25.4;6*25.4 86*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl18  = Block( (1/1000)*[0 86*25.4  3*25.4;6*25.4 98*25.4   11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl19  = Block( (1/1000)*[0 98*25.4  3*25.4;6*25.4 116*25.4  11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)
bl20  = Block( (1/1000)*[0 116*25.4 3*25.4;6*25.4 136*25.4  11*25.4], nx=1, ny=1, nz=1, cellshape=HEX8)

bl21 = BlockInset( (1/1000)*[3*25.4 0 3*25.4; 3*25.4 136*25.4 3*25.4], curvetype="polyline")


mesh=Mesh(bl1, bl2, bl3, bl4, bl5, bl6, bl7, bl8, bl9, bl10, bl21, bl11, bl12, bl13, bl14, bl15, bl16, bl17, bl18, bl19, bl20, quiet=true)


# Model definition
mats = [
    :solids => MechSolid => DruckerPrager => (E=Ec, nu=0.2, rho=roc, alpha=0.4312, kappa=3771.2),
    :lines  => MechBar => LinearElastic => (E=Ea, A=Aa, rho=roa),
    :joints1d => MechBondSlip => ElasticBondSlip => (ks=1.e7, kn=1.e7, p=p),
    # MaterialBind(:solids, DruckerPrager => (E=Ec, nu=0.2, rho=roc, alpha=0.4312, kappa=3771.2) ),
    # MaterialBind(:joints1D, ElasticJoint1D(ks=1.e7, kn=1.e7, A=Aa) ),
    # MaterialBind(:lines   , ElasticRod => (E=Ea, A=Aa, rho=roa) ),
]

matl = MaterialList()
matl.add_material(:solids, MechSolid, DruckerPrager, E=Ec, nu=0.2, rho=roc, alpha=0.4312, kappa=3771.2)
matl.add_material(:lines,  MechBar, LinearElastic, E=Ea, A=Aa, rho=roa)
matl.add_material(:joints1d, MechBondSlip, ElasticBondSlip, ks=1.e7, kn=1.e7, p=p)

ctx   = MechContext()
model = FEModel(mesh, mats, ctx)
ana   = MechModalAnalysis(model)


# model = FEModel(mesh, mats)

# Finite element modeling
P = (((-13.5)/0.00022480894387096)/1000)
bcs = [
    (y==0, z==3*0.0254) => NodeBC(uz = 0, uy=0, ux=0 ),
    (y==136*0.0254, z==3*0.0254) => NodeBC(uz=0 ),
]

addstage!(ana, bcs)
run!(ana, nmods=5, rayleigh=false)
