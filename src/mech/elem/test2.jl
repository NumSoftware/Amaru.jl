using Amaru

bl1 = Block( [-300 0 0; -300 600 0; -300 300 0], nx=10,nz=10, cellshape=LIN2, tag="shell")

msh = Mesh(bl1, verbose=true)
msh = revolve(msh, angle=180, n=40)

# Finite element model

mats = [ "shell" => ElasticShellQUAD4(E=300, nu=0.3, thick = 3) ]

dom = Domain(msh, mats)

log1 = NodeGroupLogger()
loggers = [:(y ==300 && z==300) => log1]
setloggers!(dom, loggers)


        bcs =
            [
            :(y==0) => NodeBC( ux=0, uz=0, rx=0),
            :(y==600) => NodeBC(ux=0, uz=0, rx =0),
            :(y==300 && z==300) => NodeBC(fz=-1),
            ]

solve!(dom, bcs, tol = 0.1, nincs = 1,verbose=true)
save(dom, "shell4node_coberta.vtu")
save(log1, "shell4node_coberta.dat")

