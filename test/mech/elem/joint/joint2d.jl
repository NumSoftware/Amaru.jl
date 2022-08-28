using Amaru
using Test

# mesh generation

bl  = Block( [0 0; 0.2 0.1], nx=2, ny=6, cellshape=QUAD4, tag="solids")
msh = Mesh(bl)
generate_joints!(msh, tag="joints")

# finite element analysis

E = 27.e6

for i in 1:2

    mats = [
        "solids" => ElasticSolid(E=E, nu=0.2),
        "joints" => MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, zeta=5.0, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk" ),
        "joints" => ElasticJoint(E=E, nu=0.2, zeta=5.0)
    ]

    if i==1
        mats = mats[[1,2]]
    else
        mats = mats[[1,3]]
    end

    model = Model(msh, mats, modeltype="plane-stress", thickness=1.0)

    # Loggers
    tag!(model.elems["joints"].ips, "jips")
    log1 = IpLogger()
    loggers = [
               "jips" => log1
              ]
    setloggers!(model, loggers)

    # Boundary conditions
    bcs = [
           :(x==0)   => SurfaceBC(ux=0, uy=0 ),
           :(x==0.2) => SurfaceBC(ux=2.0*1.7e-4),
          ]
    addstage!(model, bcs, nincs=20, nouts=10)
    @test solve!(model, autoinc=true, maxits=3, tol=0.01, report=true, scheme="ME").success

    if @isdefined(makeplots) && makeplots
        using PyPlot
        table = log1.table
        #plot(table[:up], table[:s1])
        plot(table[:jw1], table[:js1])
    end
end
