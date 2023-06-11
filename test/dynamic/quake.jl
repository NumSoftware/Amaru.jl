using Amaru

#Units kN,m,kPa


# Input data
fc = 40     
Ec0 = 21.5e6
ae  = 0.9
fcm = fc+8
Eci = Ec0*ae*(fcm/10)^(1/3)
ai  = 0.8+0.2*fcm/88
Ec  = ai*Eci
#Ec = 30e6
roc = 2.6

#Beam1

bl1  = Block( [0 0 2.7;0.3 4 2.9], nx=1, ny=10, nz=1, cellshape=HEX20)
bl2  = Block( [0 0 2.9;0.3 4 3.0], nx=1, ny=10, nz=1, cellshape=HEX20)

#Beam2
bl3  = Block( [4.3 0 2.7;4.6 4 2.9], nx=1, ny=10, nz=1, cellshape=HEX20)
bl4  = Block( [4.3 0 2.9;4.6 4 3], nx=1, ny=10, nz=1, cellshape=HEX20)

#Beam3
bl5  = Block( [0.3 -0.3 2.7;4.3 0.0 2.9], nx=10, ny=1, nz=1, cellshape=HEX20)
bl6  = Block( [0.3 -0.3 2.9;4.3 0.0 3], nx=10, ny=1, nz=1, cellshape=HEX20)
#Beam4
bl7  = Block( [0.3 4.0 2.7;4.3 4.3 2.9], nx=10, ny=1, nz=1, cellshape=HEX20)
bl8  = Block( [0.3 4.0 2.9;4.3 4.3 3], nx=10, ny=1, nz=1, cellshape=HEX20)
#Column1
bl9  = Block( [0. -0.3 0.;0.30 0.0 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#Column2
bl10  = Block( [0 4.0 0;0.30 4.3 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#Column3
bl11  = Block( [4.3 -0.3 0.;4.6 0.0 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#Column4
bl12  = Block( [4.3 4.0 0;4.60 4.3 2.7], nx=1, ny=1, nz=9, cellshape=HEX20)

#Plate
bl13  = Block( [0.3 0 2.9;4.3 4 3], nx=10, ny=10, nz=1, cellshape=HEX20)


#Adicional Blocks to suports

bl14  = Block( [0 -0.3 2.7;0.3 0.0 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl15  = Block( [0 -0.3 2.9;0.3 0.0 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl16  = Block( [0 4.0 2.7;0.3 4.3 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl17  = Block( [0 4.0 2.9;0.3 4.3 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl18  = Block( [4.3 -0.3 2.7;4.6 0.0 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl19  = Block( [4.3 -0.3 2.9;4.6 0.0 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)

bl20  = Block( [4.3 4.0 2.7;4.6 4.3 2.9], nx=1, ny=1, nz=1, cellshape=HEX20)
bl21  = Block( [4.3 4.0 2.9;4.6 4.3 3.0], nx=1, ny=1, nz=1, cellshape=HEX20)


mesh=Mesh(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10,bl11,bl12,bl13,bl14,bl15,bl16,bl17,bl18,bl19,bl20,bl21)


#Definition of the material

materials = [
    MaterialBind(:solids, ElasticSolid(E=Ec, nu=0.2, rho=roc))
]


model = FEModel(mesh,materials)


#Definition of the boundary conditions
bcs = [
    NodeBC( :(z==0), :( uz = 0, uy=0, ux=0 ))
]



dynsolve!(model,tspan=33.0, sism=true, tds=30.9, tss=0.0, nincs=50, nouts=10, autoinc=true, alpha=0.0, beta=0.0)

#@test model.nodes[25].dofdict[:uz].vals[:uz] ≈ -0.0123
