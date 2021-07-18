using Amaru


# Mesh generation
# ===============

blocks = [
    Block( [0 0; 1 10], nx=1, ny=10, cellshape=QUAD4, tag="solids"),
]

msh = Mesh(blocks, printlog=false)
generate_joints!(msh, layers=3, tag="joints")


# Finite element analysis
# =======================

# Analysis data
load  = 10.0   
E     = 1300.0  # Young modulus
kappa = 4.5e-13 # Intrinsic permeability
n     = 0.2975  # porodity 
alpha = 1.0     # biot constant 
nu    = 0.40    # Poisson
eta   = 1e-6    # water viscosity
Ks    = 1e9     # Bulk modulus of solid
Kw    = 2e6     # Bulk modulus of fluid
gw    = 10.0    # water specific weight
hd    = 10.0    # drainage height
mv    = (1+nu)*(1-2*nu)/(E*(1-nu))
k = (kappa*gw)/eta # permeability
cv    = k/(mv*gw)   # consolidation coefficient
T = [ 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0 ]
times = T*(hd^2/cv)
t1 = times[1]/10
times .+= t1

S = (alpha - n)/Ks + n/Kw 

for i=1:2
    materials = [
        "solids" => ElasticSolidLinSeep(E=E, nu=nu, k=k, gammaw=gw, alpha=alpha, S=S),
        "joints" => ElasticJointSeep(E=E, nu=nu, zeta=100, gammaw=gw,eta=eta, kt=1),
        "joints" => MCJointSeep(E=E, nu=nu, ft=2.4e3, mu=1.4, zeta=100, wc=1.7e-4, ws=1.85e-5, softcurve="hordijk", gammaw=gw, eta=eta, kt=1),
    ]

    if i==1
        materials = materials[[1,2]]
    else
        materials = materials[[1,3]]
    end

    dom = Domain(msh, materials, gammaw=10)

    log1 = NodeGroupLogger()
    loggers = [
        :(x==0) => log1
    ]
    setloggers!(dom, loggers)


    # Stage 1: hydrostatic pore-pressure
    # ==================================

    tlong = 10000*hd^2/cv

    bcs = [
        :(y==0) => NodeBC(ux=0., uy=0.),
        :(x>=0) => NodeBC(ux=0.),
        :(y==10) => NodeBC(uw=0.),
    ]

    hm_solve!(dom, bcs, end_time=tlong, nincs=2, tol=1e-2, nouts=1, printlog=false)

    dom.env.t = 0.0


    # Stage 2: loading
    # ================

    pt(t) = -load
    #pt(t) = t>=t1? -load : -load/t1*t

    bcs = [
        :(y==0)  => NodeBC(ux=0., uy=0.),
        :(x>=0)  => NodeBC(ux=0.),
        :(y==10) => FaceBC(ty=-load),
        :(y==10) => NodeBC(uw=0.),
    ]

    hm_solve!(dom, bcs, end_time=t1, nincs=4, tol=1e-2, nouts=1, printlog=false)


    # Stage 3: draining
    # =================

    bcs = [
        :(y==0)  => NodeBC(ux=0., uy=0.),
        :(x>=0)  => NodeBC(ux=0.),
        :(y==10) => FaceBC(ty=-load),
        :(y==10) => NodeBC(uw=0.),
    ]

    Uw_vals = [] # A list with porepressure vectors
    for t in times[1:end]
        hm_solve!(dom, bcs, end_time=t, nincs=20, tol=1e-2, nouts=1, printlog=false)
        push!( Uw_vals, log1.book[end][:uw] )
    end


    # Output
    # ======

    if Amaru.config.makeplots
        # Terzaghi's 1-d consolidation
        function calc_Ue(Z, T)
            sum = 0.0
            for i=0:4
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
