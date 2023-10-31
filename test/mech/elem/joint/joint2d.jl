using Amaru
using Test

# mesh generation

bl1  = Block( [0 0; 0.075 0; 0.125 0.1; 0 0.1], nx=1, ny=1, cellshape=QUAD4)
bl2  = Block( [0.075 0; 0.2 0; 0.2 0.1; 0.125 0.1], nx=1, ny=1, cellshape=QUAD4)

# straight
bl1  = Block( [0 0; 0.1 0; 0.1 0.1; 0 0.1], nx=1, ny=1, cellshape=QUAD4)
bl2  = Block( [0.1 0; 0.2 0; 0.2 0.1; 0.1 0.1], nx=1, ny=1, cellshape=QUAD4)

msh = Mesh(bl1, bl2)


generate_joints!(msh, tag="joints")

# msh.nodes[5].coord[1] = 0.15
# msh.nodes[6].coord[1] = 0.15

tag!(msh.elems[BULKCELL, :(x<=0.125)].nodes, "left")
tag!(msh.elems[BULKCELL, :(x>=0.075)].nodes, "right")

# finite element analysis

E = 27.e6

allmats = [
    :bulks << MechSolid << LinearElastic << (E=E, nu=0.2),
    # :joints << MCJoint << (E=E, nu=0.2, ft=2.4e3, mu=1.4, zeta=5.0, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk" ),
    # :joints << MCJoint << (E=E, nu=0.2, ft=2.4e3, mu=1.4, zeta=5.0, wc=1.7e-4, ws=1.85e-5, softcurve="soft" ),

    :joints << MechJoint << TCJoint << (E=E, nu=0.2, fc=-24e3, ft=2.4e3, zeta=5.0, wc=1.7e-4, alpha=1.5, gamma=0.1, theta=1.5, softcurve="hordijk" ),
    # :joints << TCJoint << (E=E, nu=0.2, fc=-24e3, ft=2.4e3, zeta=5.0, wc=1.7e-4, alpha=1.5, gamma=0.1, theta=1.5, softcurve="soft" ),

    :joints << MechJoint << TCFJoint << (E=E, nu=0.2, fc=-24e3, ft=2.4e3, zeta=5.0, wc=1.7e-4, alpha=0.3, mu=0.1, softcurve="hordijk" ),
    # :joints << TCFJoint << (E=E, nu=0.2, fc=-24e3, ft=2.4e3, zeta=5.0, wc=1.7e-4, alpha=0.3, mu=0.1, softcurve="soft" ),
    # :joints << ElasticJoint << (E=E, nu=0.2, zeta=5.0)
]

for i in (2,3)

    mats = allmats[[1,i]]

    ana = MechAnalysis(stressmodel="plane-stress", thickness=1.0)
    model = FEModel(msh, mats, ana)

    # Loggers
    tag!(model.elems["joints"].ips, "jips")
    log1 = IpLogger()
    loggers = [
               "jips" << log1
              ]
    setloggers!(model, loggers)

    # Boundary conditions
    bcs = [
        "left" << NodeBC(ux=0, uy=0),
        
        # "right" << NodeBC(ux=-2e-5, uy=2e-4),
        # "right" << NodeBC(ux=2e-5, uy=1e-5),
        "right" << NodeBC(ux=1e-4, uy=3e-4),

        # "right" << NodeBC(ux=2e-4),
        # "right" << NodeBC(ux=-5e-5),
        # "right" << NodeBC(uy=2e-4),
    ]
    addstage!(model, bcs, nincs=10, nouts=20)
    # @test solve!(model, autoinc=true, maxits=3, tol=0.01, scheme="Ralston").success

    quiet = false
    quiet = true
    solve!(model, autoinc=true, maxits=3, tol=0.01, rspan=0.01, scheme="Ralston", quiet=quiet)

    if Amaru.makeplots
        using PyPlot
        table = log1.table
        # plot(table[:jup], table[:js1], marker="o")
        # plot(table[:jup], table[:js2], marker="o")

        # plot(table[:jw1], table[:js1], marker="o")
        # plot(table[:jw2], table[:js2], marker="o")
        plot(table[:js1], table[:js2], marker="o")
    end


end
