using Amaru

# Entrada de dados
# Concreto
fc = 25800000/1000                    # resistência à compressão
Ec  = (6100/0.14503773800722)*1000
Ea = (30000/0.14503773800722)*1000
Aa = 2*(0.0254*0.0254)
P = ((-13.5)/0.00022480894387096)/1000
roc = 2626.96/1000#((0.217e-6)/(1.1228705601979e-6))
roa = 7850/1000


# Mesh generation


bl1  = Block3D( (1/1000)*[0   0      0;6*25.4 20*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl2  = Block3D( (1/1000)*[0 20*25.4  0;6*25.4 38*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl3  = Block3D( (1/1000)*[0 38*25.4  0;6*25.4 50*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl4  = Block3D( (1/1000)*[0 50*25.4  0;6*25.4 59*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl5  = Block3D( (1/1000)*[0 59*25.4  0;6*25.4 68*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl6  = Block3D( (1/1000)*[0 68*25.4  0;6*25.4 77*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl7  = Block3D( (1/1000)*[0 77*25.4  0;6*25.4 86*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl8  = Block3D( (1/1000)*[0 86*25.4  0;6*25.4 98*25.4   3*25.4], nx=1, ny=1, shape=HEX8)
bl9  = Block3D( (1/1000)*[0 98*25.4  0;6*25.4 116*25.4  3*25.4], nx=1, ny=1, shape=HEX8)
bl10 = Block3D( (1/1000)*[0 116*25.4 0;6*25.4 136*25.4  3*25.4], nx=1, ny=1, shape=HEX8)

bl11  = Block3D( (1/1000)*[0   0      3*25.4;6*25.4 20*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl12  = Block3D( (1/1000)*[0 20*25.4  3*25.4;6*25.4 38*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl13  = Block3D( (1/1000)*[0 38*25.4  3*25.4;6*25.4 50*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl14  = Block3D( (1/1000)*[0 50*25.4  3*25.4;6*25.4 59*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl15  = Block3D( (1/1000)*[0 59*25.4  3*25.4;6*25.4 68*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl16  = Block3D( (1/1000)*[0 68*25.4  3*25.4;6*25.4 77*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl17  = Block3D( (1/1000)*[0 77*25.4  3*25.4;6*25.4 86*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl18  = Block3D( (1/1000)*[0 86*25.4  3*25.4;6*25.4 98*25.4   11*25.4], nx=1, ny=1, shape=HEX8)
bl19  = Block3D( (1/1000)*[0 98*25.4  3*25.4;6*25.4 116*25.4  11*25.4], nx=1, ny=1, shape=HEX8)
bl20  = Block3D( (1/1000)*[0 116*25.4 3*25.4;6*25.4 136*25.4  11*25.4], nx=1, ny=1, shape=HEX8)

bl21 = BlockInset( (1/1000)*[3*25.4 0 3*25.4; 3*25.4 136*25.4 3*25.4] , curvetype="polyline")


mesh=Mesh(bl1, bl2, bl3, bl4, bl5, bl6, bl7, bl8, bl9, bl10, bl21, bl11, bl12, bl13, bl14, bl15, bl16, bl17, bl18, bl19, bl20,
verbose=true)  

#tag!(mesh.points[1], 1)
#tag!(mesh.cells    , 10) # all cells
#tag!(mesh.ips[:(x>0)], 1000)


# Domain definition

materials = [
    MaterialBind(:solids, DruckerPrager(E=Ec, nu=0.2, rho=roc, alpha=0.4312, kappa=3771.2) ),
    MaterialBind(:joints1D, ElasticJoint1D(ks=1.e7, kn=1.e7, A=Aa) ),
    MaterialBind(:lines   , ElasticRod(E=Ea, A=Aa, rho=roa) ),

]

dom = Domain(mesh, materials)

# Finite element modeling
P = (((-13.5)/0.00022480894387096)/1000)
bcs = [
    NodeBC(:(y==0  && z==3*0.0254), :( uz = 0, uy=0, ux=0 )),
    NodeBC(:(y==136*0.0254  && z==3*0.0254), :(uz=0 )),
    NodeBC(:(x==0 && z==11*0.0254 && y==50*0.0254), fz=:( (((-13.5)/0.00022480894387096)/1000)/4) ),
    NodeBC(:(x==6*0.0254 && z==11*0.0254 && y==50*0.0254), fz=:( (((-13.5)/0.00022480894387096)/1000)/4) ),
    NodeBC(:(x==0 && z==11*0.0254 && y==86*0.0254), fz=:( (((-13.5)/0.00022480894387096)/1000)/4) ),
    NodeBC(:(x==6*0.0254 && z==11*0.0254 && y==86*0.0254), fz=:( (((-13.5)/0.00022480894387096)/1000)/4) ),

    ]


dynsolve!(dom, bcs, time_span=0.05, nincs=1000, verbose=true, filekey="dyn", nouts=100, autoinc=true, alpha=0.0, beta=0.0)
#dynsolve!(dom, bcs, time_span=0.1, nincs=50, verbose=true, filekey="dyn", nouts=100, autoinc=true)
#using Glob
#rm.(glob("*.vtk"))
