# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PlateMZC

mutable struct PlateMZC<:Mechanical
    id    ::Int
    shape ::ShapeType

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function PlateMZC()
        return new()
    end
end

matching_shape_family(::Type{PlateMZC}) = LINE_SHAPE

function plate_B_matrix(elem::PlateMZC)
        C  = get_coords(elem)
        a  = C[2,1]-C[1,1] # element length in X direction
        b  = C[2,2]-C[2,1] # element length in Y direction

        gauss_x = zeros(4,1)
        gauss_y = zeros(4,1)

        gauss_x[1] = -1/sqrt(3);
        gauss_y[1] = -1/sqrt(3);

        gauss_x[2] =  1/sqrt(3);
        gauss_y[2] =-1/sqrt(3);

        gauss_x[3] = 1/sqrt(3);
        gauss_y[3] = 1/sqrt(3);

        gauss_x[4] =-1/sqrt(3);
        gauss_y[4] = 1/sqrt(3);


    for igaus = 1 : 4
        x = gauss_x(igaus) # x = Local X coordinate of the Gauss point
        y = gauss_y(igaus) # y = Local Y coordinate of the Gauss point

        d2N = zeros(4,3)
        d2N[1,1] = 3*( x - x*y )/(4*a^2);
        d2N[2,1] = 3*(-x + x*y )/(4*a^2);
        d2N[3,1] = 3*(-x - x*y )/(4*a^2);
        d2N[4,1] = 3*( x + x*y )/(4*a^2);

        d2N[1,2] = 3*( y - x*y )/(4*b^2);
        d2N[2,2] = 3*( y + x*y )/(4*b^2);
        d2N[3,2] = 3*(-y - x*y )/(4*b^2);
        d2N[4,2] = 3*(-y + x*y )/(4*b^2);

        d2N[1,3] = 2*(  1/2 - 3*x^2/8 - 3*y^2/8)/(a*b);
        d2N[2,3] = 2*( -1/2 + 3*x^2/8 + 3*y^2/8)/(a*b);
        d2N[3,3] = 2*(  1/2 - 3*x^2/8 - 3*y^2/8)/(a*b);
        d2N[4,3] = 2*( -1/2 + 3*x^2/8 + 3*y^2/8)/(a*b);
        #------------------------------------------------------
        d2NN = zeros(4,3)
        d2NN[1,1] = ( (3*a*x - 3*a*x*y - a + a*y)/4 )/a^2;
        d2NN[2,1] = ( (3*a*x - 3*a*x*y + a - a*y)/4 )/a^2;
        d2NN[3,1] = ( (3*a*x + 3*a*x*y + a + a*y)/4 )/a^2;
        d2NN[4,1] = ( (3*a*x + 3*a*x*y - a - a*y)/4 )/a^2;

        d2NN[1,2] = 0;
        d2NN[2,2] = 0;
        d2NN[3,2] = 0;
        d2NN[4,2] = 0;

        d2NN[1,3] = 2*( -3/8*a*x^2 + a*x/4 + a/8 )/(a*b);
        d2NN[2,3] = 2*( -3/8*a*x^2 - a*x/4 + a/8 )/(a*b);
        d2NN[3,3] = 2*(  3/8*a*x^2 + a*x/4 - a/8 )/(a*b);
        d2NN[4,3] = 2*(  3/8*a*x^2 - a*x/4 - a/8 )/(a*b);
        #------------------------------------------------------
        d2NNN = zeros(4,3)
        d2NNN[1,1] = 0;
        d2NNN[2,1] = 0;
        d2NNN[3,1] = 0;
        d2NNN[4,1] = 0;

        d2NNN[1,2] = ( (3*b*y - 3*b*x*y - b + b*x)/4 )/b^2;
        d2NNN[2,2] = ( (3*b*y + 3*b*x*y - b - b*x)/4 )/b^2;
        d2NNN[3,2] = ( (3*b*y + 3*b*x*y + b + b*x)/4 )/b^2;
        d2NNN[4,2] = ( (3*b*y - 3*b*x*y + b - b*x)/4 )/b^2;

        d2NNN[1,3] = 2*( -3/8*b*y^2 + b*y/4 + b/8 )/(a*b);
        d2NNN[2,3] = 2*(  3/8*b*y^2 - b*y/4 - b/8 )/(a*b);
        d2NNN[3,3] = 2*(  3/8*b*y^2 + b*y/4 - b/8 )/(a*b);
        d2NNN[4,3] = 2*( -3/8*b*y^2 - b*y/4 + b/8 )/(a*b);
        #-----------------------------------------------------
        bmat_1  = [ -d2N[1,1] -d2NN[1,1] -d2NNN[1,1]
                    -d2N[1,2] -d2NN[1,2] -d2NNN[1,2]
                    -d2N[1,3] -d2NN[1,3] -d2NNN[1,3]];

        bmat_2  = [ -d2N[2,1] -d2NN[2,1] -d2NNN[2,1]
                    -d2N[2,2] -d2NN[2,2] -d2NNN[2,2]
                    -d2N[2,3] -d2NN[2,3] -d2NNN[2,3]];

        bmat_3  = [ -d2N[3,1] -d2NN[3,1] -d2NNN[3,1]
                    -d2N[3,2] -d2NN[3,2] -d2NNN[3,2]
                    -d2N[3,3] -d2NN[3,3] -d2NNN[3,3]];

        bmat_4  = [ -d2N[4,1] -d2NN[4,1] -d2NNN[4,1]
                    -d2N[4,2] -d2NN[4,2] -d2NNN[4,2]
                    -d2N[4,3] -d2NN[4,3] -d2NNN[4,3]];

        Bb = [bmat_1 bmat_2 bmat_3 bmat_4]; # strain-displacement matrix
    end
   return Bb
end

function D_matrix(elem::PlateMZC)
    th     = elem.env.thickness
    coef = elem.mat.E*th^3/(12*(1-elem.mat.nu^2));

    D_mat = coef*[1 elam.mat.nu 0
                  elam.mat. 1 0
                  0  0 (1-elem.mat.nu)/2];
    return D_mat
end

function elem_config_dofs(elem::PlateMZC)
    ndim = elem.env.ndim
    ndim == 1 && error("PlateMZC: Plate elements do not work in 1d analyses")
    if ndim==2
        for node in elem.nodes
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :uz, :fz)
        end
    else
        error("PlateMZC: Plate elements do not work in this analyses")
        #=
        for node in elem.nodes

            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
        =#
    end
end

function elem_map(elem::PlateMZC)::Array{Int,1}
    #=
    if elem.env.ndim==2
        dof_keys = (:ux, :uy, :rz)
    else
        dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    end
        =#

    dof_keys = (:uz, :rx, :ry)
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

# Return the class of element where this material can be used
#client_shape_class(mat::PlateMZC) = LINE_SHAPE
#=
function calcT(elem::PlateMZC, C)
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,1])/L
    return

end
=#

function elem_stiffness(elem::PlateMZC)
#=
    C  = get_coords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    L3 = L*L*L
    mat = elem.mat
    EA = mat.E*mat.A
    EI = mat.E*mat.I

    K0 = [ EA/L     0         0         -EA/L    0         0
           0       12*EI/L3   6*EI/L2    0     -12*EI/L3   6*EI/L2
           0        6*EI/L2   4*EI/L     0      -6*EI/L2   2*EI/L
          -EA/L     0          0         EA/L     0        0
           0      -12*EI/L3  -6*EI/L2    0      12*EI/L3  -6*EI/L2
           0        6*EI/L2   2*EI/L     0      -6*EI/L2   4*EI/L  ]

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L

    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]
=#
    nnodes = length(elem.nodes)
    th     = elem.env.thickness
    C  = get_coords(elem)
    a  = C[2,1]-C[1,1] # element length in X direction
    b  = C[2,2]-C[2,1] # element length in Y direction

    # K_elem = zeros( nnodes*3 , nnodes*3 )

    Bb = plate_B_matrix(elem)
    D_mat = D_matrix(elem)

    # K_elem = K_elem + bmat'*D_mat*bmat_b*a*b
    K_elem = bmat'*D_mat*bmat_b*a*b
    map = elem_map(elem)

    return K_elem, map, map
end

#=
function elem_mass(elem::PlateMZC)
    C  = get_coords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    mat = elem.mat
    EA = mat.E*mat.A
    EI = mat.E*mat.I


    M0 = mat.ρ*L/420.0*[ 140   0      0      70    0      0
                         0     156    22*L   0     54    -13*L
                         0     22*L   4*L2   0     13*L  -3*L2
                         70    0      0      140   0      0
                         0     54     13*L   0     156   -22*L
                         0    -13*L  -3*L2   0    -22*L   4*L2 ]

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L
    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]

    map = elem_map(elem)
    return T'*M0*T, map, map
end
=#

function elem_update!(elem::PlateMZC, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
end
