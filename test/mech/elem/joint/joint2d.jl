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
ft = 2.4e3
fc = -24e3

# fc = -30e3
# ft = 2.5e3
alpha=0.1

allmats = [
    :bulks => MechSolid => LinearElastic => (E=E, nu=0.2),
    # :joints => MCJoint => (E=E, nu=0.2, ft=ft, mu=1.4, zeta=5.0, wc=1.7e-4, ws=1.85e-5 ),
    # :joints => MCJoint => (E=E, nu=0.2, ft=ft, mu=1.4, zeta=5.0, wc=1.7e-4, ws=1.85e-5, softcurve="soft" ),

    :joints => MechJoint => TCJoint => (E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, wc=1.7e-4, alpha=1.5, gamma=0.05, theta=1.5 ),
    # :joints => MechJoint => TC2Joint => (E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, wc=1.7e-4, alpha=0.5, gamma=0.3, theta=1 ),
    # :joints => MechJoint => TCJoint => (E=E, nu=0.2, fc=-24e3, ft=ft, zeta=5.0, softcurve=[0 ft; 0.15*1.5e-4 1.1*ft; 0.3*1.5e-4 0.7*ft; 1.5e-4 0.5ft], alpha=1.5, gamma=0.1, theta=1.5 ),
    # :joints => TCJoint => (E=E, nu=0.2, fc=-24e3, ft=ft, zeta=5.0, wc=1.7e-4, alpha=1.5, gamma=0.1, theta=1.5, softcurve="soft" ),

    # :joints => MechJoint => TCFJoint => (E=E, nu=0.2, fc=-24e3, ft=ft, zeta=5.0, wc=1.7e-4, alpha=0.3, mu=0.1 ),
    # :joints => TCFJoint => (E=E, nu=0.2, fc=-24e3, ft=ft, zeta=5.0, wc=1.7e-4, alpha=0.3, mu=0.1, softcurve="soft" ),
    # :joints => ElasticJoint => (E=E, nu=0.2, zeta=5.0)
]

# for i in (2,3)
for i in 2

    mats = allmats[[1,i]]

    ana = MechAnalysis(stressmodel=:planestress, thickness=1.0)
    model = FEModel(msh, mats, ana)

    # Loggers
    tag!(model.elems["joints"].ips, "jips")
    log1 = IpLogger()
    loggers = [
               "jips" => log1
              ]
    setloggers!(model, loggers)

    # Boundary conditions
    bcs = [
        "left" => NodeBC(ux=0, uy=0),
        
        "right" => NodeBC(ux=-1e-9, uy=8e-5),
        # "right" => NodeBC(ux=2e-5, uy=1e-5),
        # "right" => NodeBC(ux=1e-4, uy=1e-5),
        # "right" => NodeBC(ux=2e-4),

        # "right" => NodeBC(ux=2e-4),
        # "right" => NodeBC(ux=-5e-5),
        # "right" => NodeBC(uy=2e-4),
    ]
    addstage!(model, bcs, nincs=80, nouts=20)

    solve!(model, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, scheme=:Ralston)

    if Amaru.makeplots
        table = log1.table

        chart = Chart(; xlabel=L"$u_p$", ylabel=L"\sigma")
        addplot!(chart, LinePlot(table[:jup], table[:js1], marker=:circle))
        save(chart, "up-sn.pdf")
        
        chart = Chart(; xlabel=L"\sigma_n", ylabel=L"\tau")
        addplot!(chart, LinePlot(table[:js1], table[:js2], marker=:circle))
        save(chart, "sn-tau.pdf")

    end


end
