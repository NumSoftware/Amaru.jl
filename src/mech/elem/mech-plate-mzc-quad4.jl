# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PlateMZC

mutable struct PlateMZCElem<:MechElem
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function PlateMZC()
        return new()
    end
end

matching_shape_family(::Type{PlateMZCElem}) = BULKCELL


function D_matrix(elem::PlateMZCElem)

    coef = elem.matparams.E*elem.matparams.thick^3/(12*(1-elem.matparams.nu^2));

    D_mat = coef*[1 elem.matparams.nu 0
                  elem.matparams.nu 1 0
                  0  0 (1-elem.matparams.nu)/2];
    return D_mat
end


function elem_config_dofs(elem::PlateMZCElem)
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


function elem_map(elem::PlateMZCElem)::Array{Int,1}

    #if elem.env.ndim==2
        dof_keys = (:uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end


function elem_stiffness(elem::PlateMZCElem)

    nnodes = length(elem.nodes)
    C  = getcoords(elem)
    a  = abs(C[2,1]-C[1,1])/2 # element length in X direction
    b  = abs(C[2,1]-C[1,1])/2# element length in Y direction

    d2N = zeros(4,3)
    d2NN = zeros(4,3)
    d2NNN = zeros(4,3)
    K_elem = zeros( nnodes*3 , nnodes*3 )

    D_mat = D_matrix(elem)


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
        x = gauss_x[igaus] # x = Local X coordinate of the Gauss point
        y = gauss_y[igaus] # y = Local Y coordinate of the Gauss point

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

        Bb = [bmat_1 bmat_2 bmat_3 bmat_4];

        K_elem += Bb'*D_mat*Bb*a*b

    end
    map = elem_map(elem)

    return K_elem, map, map
end


function update_elem!(elem::PlateMZCElem, U::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
end
