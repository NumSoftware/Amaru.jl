# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PlateRM8node

mutable struct PlateRM8node<:Mechanical
    id    ::Int
    shape ::CellShape

    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function PlateRM8node()
        return new()
    end
end

matching_shape_family(::Type{PlateRM8node}) = SOLID_CELL

function D_matrix(elem::PlateRM8node)

    coef = elem.mat.E/(1-elem.mat.nu^2);

    D_mat = coef*[1 elem.mat.nu 0 0 0
                  elem.mat.nu 1 0 0 0
                  0  0 (1/2)*(1-elem.mat.nu) 0 0
                  0  0 0 (5/12)*(1-elem.mat.nu) 0
                  0  0 0 0 (5/12)*(1-elem.mat.nu)];
    return D_mat
end


function elem_config_dofs(elem::PlateRM8node)
    ndim = elem.env.ndim
    ndim == 1 && error("PlateRM8node: Plate elements do not work in 1d analyses")
    if ndim==2
        for node in elem.nodes
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :uz, :fz)
        end
    else
        error("PlateRM8node: Plate elements do not work in this analyses")
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

function elem_map(elem::PlateRM8node)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function elem_stiffness(elem::PlateRM8node)

    nnodes = length(elem.nodes)
    th     = 0.15 # COLOCAR AUTOMÁTICO
    C  = getcoords(elem)
    Area = cell_extent(Area)

    a  = abs(C[2,1]-C[1,1]) # element length in X direction
    b  = abs(C[2,1]-C[1,1]) # VERIFICAR element length in Y direction

    K_elem = zeros(24,24)

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
                    N1=1/4*(1-(e))*(1-(n))*(-(e)-(n)-1);
                    N2=1/4*(1+(e))*(1-(n))*((e)-(n)-1);
                    N3=1/4*(1+(e))*(1+(n))*((e)+(n)-1);
                    N4=1/4*(1-(e))*(1+(n))*(-(e)+(n)-1);
                    N5=1/2*(1-(e)^2)*(1-(n));
                    N6=1/2*(1-(n)^2)*(1+(e));
                    N7=1/2*(1-(e)^2)*(1+(n));
                    N8=1/2*(1-(n)^2)*(1-(e));


      # Element shape function partial derivative Hx
                    Hx1 = 2*(-(1-n)*(-n-e-1)/4.0-(1-e)*(1-n)/4.0)/a;
                    Hx2 = 2*((1-n)*(-n+e-1)/4.0+(e+1)*(1-n)/4.0)/a;
                    Hx3 = 2*((n+1)*(n+e-1)/4.0+(e+1)*(n+1)/4.0)/a;
                    Hx4 = 2*(-(n+1)*(n-e-1)/4.0-(1-e)*(n+1)/4.0)/a;
                    Hx5 =-2*e*(1-n)/a;
                    Hx6 =(1-n^2)/a;
                    Hx7 =-2*e*(n+1)/a;
                    Hx8 =-(1-n^2)/a;

      # Element shape function partial derivative Hy

                    Hy1 = 2*(-(1-e)*(-n-e-1)/4.0-(1-e)*(1-n)/4.0)/b;
                    Hy2 = 2*(-(e+1)*(-n+e-1)/4.0-(e+1)*(1-n)/4.0)/b;
                    Hy3 = 2*((e+1)*(n+e-1)/4.0+(e+1)*(n+1)/4.0)/b;
                    Hy4 = 2*((1-e)*(n-e-1)/4.0+(1-e)*(n+1)/4.0)/b;
                    Hy5 = -(1-e^2)/b;
                    Hy6 = -2*(e+1)*n/b;
                    Hy7 = (1-e^2)/b;
                    Hy8 = -2*(1-e)*n/b;


                    Jac=[a/2   0  ;
                          0   b/2];

      #Element connection matrix

                    B=[0     0     z*Hx1   0     0     z*Hx2   0    0      z*Hx3    0      0      z*Hx4    0     0     z*Hx5   0     0     z*Hx6   0    0      z*Hx7    0      0      z*Hx8
                       0  -z*Hy1     0     0   -z*Hy2    0     0   -z*Hy3    0      0    -z*Hy4     0      0  -z*Hy5     0     0   -z*Hy6    0     0   -z*Hy7    0      0    -z*Hy8     0
                       0  -z*Hx1   z*Hy1   0   -z*Hx2  z*Hy2   0   -z*Hx3   z*Hy3   0    -z*Hx4   z*Hy4    0  -z*Hx5   z*Hy5   0   -z*Hx6  z*Hy6   0   -z*Hx7   z*Hy7   0    -z*Hx8   z*Hy8
                       Hy1   -N1      0    Hy2   -N2      0    Hy3  -N3       0     Hy4    -N4       0     Hy5   -N5      0    Hy6   -N6      0    Hy7  -N7       0     Hy8    -N8       0
                       Hx1    0       N1   Hx2    0      N2    Hx3    0       N3    Hx4     0       N4     Hx5    0       N5   Hx6    0      N6    Hx7    0       N7    Hx8     0       N8    ];


      #Local axis element stiffness matrix
                    K_elem += 0.5*th*det(Jac)*H[i]*H[j]*H[k]*B'*D_mat*B
                  end
              end
           end

    map = elem_map(elem)
   #println(K_elem)
  return K_elem, map, map
end


function elem_update!(elem::PlateRM8node, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    #println(K)
    dU  = U[map]
    #println(dU)
    F[map] += K*dU
    #println(F)
end
