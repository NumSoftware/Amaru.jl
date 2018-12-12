using Amaru

#Unidades kN,m,kPa


# Entrada de dados
# Concreto
fc = 40                    # resistência à compressão

# Modulo de Elasticidade
Ec0 = 21.5e6
ae  = 0.9
fcm = fc+8
Eci = Ec0*ae*(fcm/10)^(1/3)
ai  = 0.8+0.2*fcm/88
Ec  = ai*Eci
#Ec = 30e6
roc = 2.6

#VIGA1H###################

bl1  = Block3D( [0 0 2.7;0.3 4 2.9], nx=1, ny=10, nz=1, cellshape=HEX20)
bl2  = Block3D( [0 0 2.9;0.3 4 3.0], nx=1, ny=10, nz=1, cellshape=HEX20)

#VIGA2H###################
bl3  = Block3D( [4.3 0 2.7;4.6 4 2.9], nx=1, ny=10, nz=1, cellshape=HEX20)
bl4  = Block3D( [4.3 0 2.9;4.6 4 3], nx=1, ny=10, nz=1, cellshape=HEX20)

#VIGA1V###################
bl5  = Block3D( [0.3 -0.3 2.7;4.3 0.0 2.9], nx=10, ny=1, nz=1, cellshape=HEX20)
bl6  = Block3D( [0.3 -0.3 2.9;4.3 0.0 3], nx=10, ny=1, nz=1, cellshape=HEX20)
#VIGA2V###################
bl7  = Block3D( [0.3 4.0 2.7;4.3 4.3 2.9], nx=10, ny=1, nz=1, cellshape=HEX20)
bl8  = Block3D( [0.3 4.0 2.9;4.3 4.3 3], nx=10, ny=1, nz=1, cellshape=HEX20)
#COLUNA 1
bl9  = Block3D( [0. -0.3 0.;0.30 0.0 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#COLUNA 2
bl10  = Block3D( [0 4.0 0;0.30 4.3 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#COLUNA 3
bl11  = Block3D( [4.3 -0.3 0.;4.6 0.0 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#COLUNA 4
bl12  = Block3D( [4.3 4.0 0;4.60 4.3 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#LAJE
bl13  = Block3D( [0.3 0 2.9;4.3 4 3], nx=10, ny=10, nz=1, cellshape=HEX20)


#Blocos Adicionais  Apoios

bl14  = Block3D( [0 -0.3 2.7;0.3 0.0 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl15  = Block3D( [0 -0.3 2.9;0.3 0.0 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl16  = Block3D( [0 4.0 2.7;0.3 4.3 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl17  = Block3D( [0 4.0 2.9;0.3 4.3 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl18  = Block3D( [4.3 -0.3 2.7;4.6 0.0 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl19  = Block3D( [4.3 -0.3 2.9;4.6 0.0 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl20  = Block3D( [4.3 4.0 2.7;4.6 4.3 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl21  = Block3D( [4.3 4.0 2.9;4.6 4.3 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)


mesh=Mesh(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10,bl11,bl12,bl13,bl14,bl15,bl16,bl17,bl18,bl19,bl20,bl21,verbose=true)


#Definicao de material

materials = [
    MaterialBind(:solids, ElasticSolid(E=Ec, nu=0.2, rho=roc))
]


dom = Domain(mesh,materials)


#Definicao condicoes de contorno
bcs = [
    NodeBC( :(z==0), :( uz = 0, uy=0, ux=0 ))
]



dynsolve!(dom, bcs, time_span=35.0, sism=true, tds=30.9, tss=0.0, nincs=1000, verbose=true, filekey="dyn", nouts=100, autoinc=true, alpha=0.0, beta=0.0)


