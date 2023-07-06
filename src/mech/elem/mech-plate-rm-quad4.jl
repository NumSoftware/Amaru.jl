# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PlateRM

mutable struct PlateRMElem<:MechElem
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    matparams::MatParams
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function PlateRM()
        return new()
    end
end

matching_shape_family(::Type{PlateRMElem}) = BULKCELL

function D_matrix(elem::PlateRMElem)

    coef = elem.matparams.E/(1-elem.matparams.nu^2);

    D_mat = coef*[1 elem.matparams.nu 0 0 0
                  elem.matparams.nu 1 0 0 0
                  0  0 (1/2)*(1-elem.matparams.nu) 0 0
                  0  0 0 (5/12)*(1-elem.matparams.nu) 0
                  0  0 0 0 (5/12)*(1-elem.matparams.nu)];
    return D_mat
end


function elem_config_dofs(elem::PlateRMElem)
    ndim = elem.env.ndim
    ndim == 1 && error("PlateRM: Plate elements do not work in 1d analyses")
    if ndim==2
        for node in elem.nodes
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :uz, :fz)
        end
    else
        error("PlateRM: Plate elements do not work in this analyses")
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

function elem_map(elem::PlateRMElem)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function elem_stiffness(elem::PlateRMElem)

    nnodes = length(elem.nodes)
    th     = 0.15 # COLOCAR AUTOMÁTICO
    C  = getcoords(elem)

    a  = abs(C[2,1]-C[1,1])  # element length in X direction
    b  = abs(C[2,1]-C[1,1]) # ALTERAR element length in Y direction

    K_elem = zeros(12,12)

    D_mat = D_matrix(elem)

    # Gauss_Legendre_Numerical Integration constants

    H=[0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170]
    W=[0.932469514203152;
       0.661209386466265;
       0.238619186083197;
      -0.238619186083197;
      -0.661209386466265;
      -0.932469514203152];

      # This region 6 segment integration
         for k in 1:6;
            for  i=1:6;
               for j in 1:6;
                   e=W[i];
                   n=W[j];
                   z=0.5*th*W[k];

      #Element  shape function
                    N1=1/4*(1-e)*(1-n);
                    N2=1/4*(1+e)*(1-n);
                    N3=1/4*(1+e)*(1+n);
                    N4=1/4*(1-e)*(1+n);

      # Element shape function partial derivative
                    Hx1=-(1-n)/a/2.0;
                    Hx2= (1-n)/a/2.0;
                    Hx3= (n+1)/a/2.0;
                    Hx4=-(n+1)/a/2.0;
                    Hy1 = -(1-e)/b/2.0;
                    Hy2 = -(e+1)/b/2.0;
                    Hy3 =  (e+1)/b/2.0;
                    Hy4 =  (1-e)/b/2.0;

                    Jac=[a/2   0  ;
                          0   b/2];

      #Element connection matrix

                    B=[0     0     z*Hx1   0     0     z*Hx2   0    0      z*Hx3    0      0      z*Hx4   ;
                0  -z*Hy1     0     0   -z*Hy2    0     0   -z*Hy3    0      0    -z*Hy4     0     ;
                0  -z*Hx1   z*Hy1   0   -z*Hx2  z*Hy2   0   -z*Hx3   z*Hy3   0    -z*Hx4   z*Hy4   ;
               Hy1   -N1      0    Hy2   -N2      0    Hy3  -N3       0     Hy4    -N4       0     ;
               Hx1    0       N1   Hx2    0      N2    Hx3    0       N3    Hx4     0       N4    ];


      #Local axis element stiffness matrix
                    K_elem += 0.5*th*det(Jac)*H[i]*H[j]*H[k]*B'*D_mat*B
                  end
              end
           end

    map = elem_map(elem)

  return K_elem, map, map
end


function update_elem!(elem::PlateRMElem, U::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
end
