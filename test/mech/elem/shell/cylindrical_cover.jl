using Amaru

bl1 = Block( [-300 0 0; -300 600 0; -300 300 0], nx=50,nz=50, cellshape=LIN2, tag="shell")

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

solve!(dom, bcs, tol = 0.1, nincs = 1,verbose=false).success
save(dom, "shell4node_coberta.vtu")
save(log1, "shell4node_coberta.dat")

mplot(dom, "coberta_shell4node.pdf", axis=false, field="uz", fieldscale=100,
      lw=0.05, colorbarscale=0.5, colorbarlabel=raw"$u_z$ [cm]",
      figsize=(6.5,4), azim=-50, warpscale =1000, elev=20, dist = 8, colormap="coolwarm")
