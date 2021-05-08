using Amaru

coord1 = [ 5.000000000e+01     1.000000000e+02    0.000000000e+00;
           1.000000000e+02     1.000000000e+02     1.000000000e+01;
           1.000000000e+02     5.000000000e+01     0.000000000e+00 ;
           5.000000000e+01     5.000000000e+01     0.000000000e+00]

coord2 = [   0.000000000e+00     5.000000000e+01     0.000000000e+00;
             0.000000000e+00     1.000000000e+02    -1.000000000e+01;
             5.000000000e+01     1.000000000e+02     0.000000000e+00;
             5.000000000e+01     5.000000000e+01     0.000000000e+00]

coord3 = [5.000000000e+01     0.000000000e+00     0.000000000e+00 ;
          0.000000000e+00     0.000000000e+00     1.000000000e+01 ;
          0.000000000e+00     5.000000000e+01     0.000000000e+00  ;
          5.000000000e+01     5.000000000e+01     0.000000000e+00]

coord4 = [   1.000000000e+02     5.000000000e+01     0.000000000e+00;
             1.000000000e+02     0.000000000e+00    -1.000000000e+01;
             5.000000000e+01     0.000000000e+00     0.000000000e+00;
             5.000000000e+01     5.000000000e+01     0.000000000e+00]

blt = [Block(coord1, nx = 1, ny =1, cellshape = QUAD4, tag="shell"),
       Block(coord2, nx = 1, ny =1, cellshape = QUAD4, tag="shell"),
       Block(coord3, nx = 1, ny =1, cellshape = QUAD4, tag="shell"),
       Block(coord4, nx = 1, ny =1, cellshape = QUAD4, tag="shell"),
       ]

msh = Mesh(blt, silent=true)

# Finite element model

mats = [ "shell" => ElasticShellQUAD4(E=2.85E4, nu=0.4, thick = 0.8) ]

dom = Domain(msh, mats)


log1 = NodeGroupLogger()
loggers = [:(x ==50 && y==50) => log1]
setloggers!(dom, loggers)

        bcs =
            [
            :(x==0) => NodeBC(rx=0, ry=0, ux=0, uy=0, uz=0),
            :(x==100) => NodeBC(rx=0, ry=0, ux=0, uy=0, uz=0),
            :(y==0) => NodeBC(rx=0, ry=0, ux=0, uy=0, uz=0),
            :(y==100) => NodeBC(rx=0, ry=0, ux=0, uy=0, uz=0),
            :(x==50 && y==50) => NodeBC(fz=-100),
            ]

solve!(dom, bcs, tol = 0.1, nincs = 1,verbose=true)
save(dom, "shell4node_hiperbolic.vtu")
save(log1, "shell4node_hiperbolic.dat")

#=
mplot(dom, "hiperbolic_shell_.pdf", axis=false, field="uz", fieldscale=1,
      lw=0.1, colorbarscale=0.5, colorbarlabel=raw"$u_z$ [cm]",
      figsize=(6.5,4), azim=250, warpscale =500, elev=-10, dist = 7.2, colormap="coolwarm")
=#
