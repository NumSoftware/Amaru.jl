using Amaru


# Mesh generation


blocks = [
    Block( [0 0; 1 10], nx=1, ny=10, cellshape=QUAD4, tag="solids"),
]

msh = Mesh(blocks)
insert_cohesive_elements!(msh, layers=3, tag="joints")


# Finite element analysis


# Analysis data
load   = 10.0
E      = 1300.0  # Young modulus
kappa  = 4.5e-13 # Intrinsic permeability
n      = 0.2975  # porosity
alpha  = 1.0     # biot constant
nu     = 0.40    # Poisson
eta    = 1e-6    # water viscosity
Ks     = 1e9     # Bulk modulus of solid
Kw     = 2e6     # Bulk modulus of fluid
gw     = 10.0    # water specific weight
hd     = 10.0    # drainage height
mv     = (1+nu)*(1-2*nu)/(E*(1-nu))
k      = (kappa*gw)/eta # permeability
cv     = k/(mv*gw)   # consolidation coefficient
T      = [ 0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0 ]
times  = T*(hd^2/cv)
lapses = diff(times)
t1     = lapses[1]/10

S = (alpha - n)/Ks + n/Kw 

for i in 1:2
    materials = [
        # "solids" => HMSolid => LinearElasticSeep => (E=E, nu=nu, k=k, alpha=alpha, S=S),
        "solids" => HMSolid => HMCombined{ConstPermeability,LinearElastic} => (E=E, nu=nu, k=k, alpha=alpha, S=S),

        "joints" => HMJoint => ElasticJointSeep => (E=E, nu=nu, zeta=100, eta=eta, kt=1),
        # "joints" => HMJoint => HMCombined{LeakoffJoint,ElasticJoint} => (E=E, nu=nu, zeta=100, eta=eta, kt=1),

        "joints" => HMJoint => MCJointSeep => (E=E, nu=nu, ft=2.4e3, mu=1.4, zeta=100, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk", eta=eta, kt=1),
        # "joints" => HMJoint => HMCombined{LeakoffJoint,MCJoint} => (E=E, nu=nu, ft=2.4e3, mu=1.4, zeta=100, wc=1.7e-4, ws=1.85e-5, eta=eta, kt=1),
    ]

    if i==1
        materials = materials[[1,2]]
    else
        materials = materials[[1,3]]
    end

    ana = HydromechAnalysis(gammaw=gw)
    model = FEModel(msh, materials, ana)

    log1 = NodeGroupLogger()
    loggers = [
        :(x==0) => log1
    ]
    setloggers!(model, loggers)


    # Stage 1: hydrostatic pore-pressure
    

    tlong = 10000*hd^2/cv
    bcs = [
        :(y==0)  => NodeBC(ux=0, uy=0),
        :(x>=0)  => NodeBC(ux=0),
        :(y==10) => NodeBC(uw=0),
    ]
    addstage!(model, bcs, tspan=tlong, nincs=2, nouts=1)
    solve!(model, tol=1e-2)
    model.ctx.t = 0.0

    # Stage 2: loading
    bcs = [
        :(y==0)  => NodeBC(ux=0., uy=0.),
        :(x>=0)  => NodeBC(ux=0.),
        :(y==10) => SurfaceBC(ty=-load),
        :(y==10) => NodeBC(uw=0.),
    ]
    addstage!(model, bcs, tspan=t1, nincs=4, nouts=1)
    solve!(model, tol=1e-2)

    # Stage 3: draining
    bcs = [
        :(y==0)  => NodeBC(ux=0., uy=0.),
        :(x>=0)  => NodeBC(ux=0.),
        :(y==10) => SurfaceBC(ty=-load),
        :(y==10) => NodeBC(uw=0.),
    ]

    Uw_vals = [] # A list with porepressure vectors
    for tspan in lapses
        addstage!(model, bcs, tspan=tspan, nincs=20, nouts=1)
        solve!(model, tol=1e-2)
        push!( Uw_vals, log1.book[end][:uw] )
    end


    # Output
    if @isdefined(makeplots) && makeplots
        # Terzaghi's 1-d consolidation
        function calc_Ue(Z, T)
            sum = 0.0
            for i in 0:4
                M = pi/2*(2*i+1)
                sum = sum + 2/M*sin(M*Z)*exp(-M^2*T)
            end
            return sum
        end

        using PyPlot


        # numerical curves
        book  = log1.book
        Uwini = book.tables[2][:uw]       # hydrostatic porepressure
        Z     = 1 .- book.tables[2][:y]/hd # normalized depth
        for Uw in Uw_vals
            dUw = (Uw - Uwini)/load       # excess of porepressure
            plot(dUw, Z, "-o")
        end

        # analytical curves
        for Ti in T
            Ue = calc_Ue.(Z, Ti)
            plot(Ue, Z, "k")
        end

        ylim(1,0)
        show()
    end
end
