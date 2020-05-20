# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ShellQuad4node

mutable struct ShellQuad4node<:Mechanical
    id    ::Int
    shape ::ShapeType
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function ShellQuad4node()
        return new()
    end
end

matching_shape_family(::Type{ShellQuad4node}) = SOLID_SHAPE

# the strain-displacement matrix for membrane forces
function D_matrixm(elem::ShellQuad4node)

    coef1 = elem.mat.E*elem.mat.thick /(1-elem.mat.nu^2)
    coef2 = elem.mat.nu*coef1
    coef3 = elem.mat.E*elem.mat.thick /2/(1+elem.mat.nu)

        D_matm = [coef1  coef2 0
                  coef2  coef1 0
                  0      0     coef3];
    return D_matm
end

# the strain-displacement matrix for shear forces

function D_matrixs(elem::ShellQuad4node)

    coef = (5/6)*elem.mat.E*elem.mat.thick/2/(1+elem.mat.nu)

            D_mats = [coef    0
                        0     coef];
    return D_mats
end

# the strain-displacement matrix for bending moments
function D_matrixb(elem::ShellQuad4node)

    D_matm = D_matrixm(elem)

    D_matb = D_matm*(elem.mat.thick^2/12)

    return D_matb
end

# Rotation Matrix

function RotMatrix(elem::ShellQuad4node)

    cxyz  = get_coords(elem) # Como chamar as coordenadas em Z?

    v12 = zeros(3,1)
    v13 = zeros(3,1)
    vxe = zeros(3,1)
    vye = zeros(3,1)
    vze = zeros(3,1)

    v12[1] = cxyz[2,1] - cxyz[1,1]
    v12[2] = cxyz[2,2] - cxyz[1,2]
    v12[3] = 0 #cxyz[2,3] - cxyz[1,3] # ARRUMAR

    v13[1] = cxyz[3,1] - cxyz[1,1]
    v13[2] = cxyz[3,2] - cxyz[1,2]
    v13[3] = 0 #cxyz[3,3] - cxyz[1,3] # ARRUMAR

    vze[1] = v12[2]*v13[3] - v12[3]*v13[2]
    vze[2] = v12[3]*v13[1] - v12[1]*v13[3]
    vze[3] = v12[1]*v13[2] - v12[2]*v13[1]

    dz = sqrt(vze[1]^2 + vze[2]^2 + vze[3]^2);

  # Unit vector normal to element surface
    vze[1] = vze[1]/dz
    vze[2] = vze[2]/dz
    vze[3] = vze[3]/dz

  # XZ plane intesection with element surface
    vxe[1] =  1/sqrt(1+(vze[1]/vze[3])^2)
    vxe[2] =  0
    vxe[3] = -1/sqrt(1+(vze[3]/vze[1])^2)

    dd = vxe[1]*vze[1] + vxe[3]*vze[3];
    if (abs(dd) > 1e-8)
      vxe[3] = -vxe[3]
    end

    if ((vze[3] == 0) && (vze[1] == 0))
      vxe[1] =  1
      vxe[2] =  0
      vxe[3] =  0
    end

  # Vector product
    vye[1] = vze[2]*vxe[3] - vxe[2]*vze[3]
    vye[2] = vze[3]*vxe[1] - vxe[3]*vze[1]
    vye[3] = vze[1]*vxe[2] - vxe[1]*vze[2]

    dy = sqrt(vye[1]^2 + vye[2]^2 + vye[3]^2)
    vye[1] = vye[1]/dy
    vye[2] = vye[2]/dy
    vye[3] = vye[3]/dy

    if (vye[2] < 0 )
      vye[1] = -vye[1]
      vye[2] = -vye[2]
      vye[3] = -vye[3]
      vxe[1] = -vxe[1]
      vxe[2] = -vxe[2]
      vxe[3] = -vxe[3]
    end

    Rot = [ vxe[1] vxe[2] vxe[3]
           vye[1] vye[2] vye[3]
           vze[1] vze[2] vze[3] ]

    return Rot
end


function elem_config_dofs(elem::ShellQuad4node)
    ndim = elem.env.ndim
    ndim == 1 && error("ShellQuad4node: Plate elements do not work in 1d analyses")
    if ndim==2
        for node in elem.nodes
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
        end
    else
        error("ShellQuad4node: Plate elements do not work in this analyses")
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

function elem_map(elem::ShellQuad4node)::Array{Int,1}

    #if elem.env.ndim==2
        dof_keys = (:ux, :uy, :uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:ux, :uy, :uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function elem_stiffness(elem::ShellQuad4node)

    nnodes = length(elem.nodes)
    C  = get_coords(elem)
    #a  = abs(C[2,1]-C[1,1])/2 # element length in X direction
    #b  = abs(C[2,1]-C[1,1])/2# element length in Y direction

    D_matm = D_matrixm(elem)
    D_mats = D_matrixs(elem)
    D_matb = D_matrixb(elem)

    Rot = RotMatrix(elem)

    gauss_x = zeros(4,1)
    gauss_y = zeros(4,1)
    gauss_w = zeros(4,1)

    gauss_x[1] = -1/sqrt(3);
    gauss_y[1] = -1/sqrt(3);
    gauss_w[1] =  1.0;

    gauss_x[2] =  1/sqrt(3);
    gauss_y[2] = -1/sqrt(3);
    gauss_w[2] =  1.0;

    gauss_x[3] =  1/sqrt(3);
    gauss_y[3] =  1/sqrt(3);
    gauss_w[3] =  1.0;

    gauss_x[4] = -1/sqrt(3);
    gauss_y[4] =  1/sqrt(3);
    gauss_w[4] =  1.0;

    xjacm = zeros(2,2)

    dxNl = zeros(4,1)
    dyNl = zeros(4,1)
    dxN = zeros(4,1)
    dyN = zeros(4,1)

    K_elem = zeros( nnodes*5 , nnodes*5 )

    cxyz  = zeros(4,3)
    cxyz[:,1:2]  = get_coords(elem)

    ctxy = cxyz*Rot';    # Rotate coordinates to element mid plane

    x = ctxy[1:4,1]; # Local X coordinate of the element
    y = ctxy[1:4,2];  # Local Y coordinate of the element

    for igaus = 1 : 4
        #-----------------------------
      xgs = gauss_x[igaus] # Local X coordinate of the Gauss point
      ygs = gauss_y[igaus] # Local Y coordinate of the Gauss point

      dxNl[1] = (-1+ygs)/4;
      dxNl[2] = ( 1-ygs)/4;
      dxNl[3] = ( 1+ygs)/4;
      dxNl[4] = (-1-ygs)/4;

      dyNl[1] = (-1+xgs)/4;
      dyNl[2] = (-1-xgs)/4;
      dyNl[3] = ( 1+xgs)/4;
      dyNl[4] = ( 1-xgs)/4;

      xjacm[1,1] = x[1]*dxNl[1] + x[2]*dxNl[2] + x[3]*dxNl[3] + x[4]*dxNl[4]
      xjacm[1,2] = y[1]*dxNl[1] + y[2]*dxNl[2] + y[3]*dxNl[3] + y[4]*dxNl[4]
      xjacm[2,1] = x[1]*dyNl[1] + x[2]*dyNl[2] + x[3]*dyNl[3] + x[4]*dyNl[4]
      xjacm[2,2] = y[1]*dyNl[1] + y[2]*dyNl[2] + y[3]*dyNl[3] + y[4]*dyNl[4]

      xjaci = inv(xjacm);

      area = abs(xjacm[1,1]*xjacm[2,2] - xjacm[2,1]*xjacm[1,2]);

      dxN[1] = xjaci[1,1]*dxNl[1]+xjaci[1,2]*dyNl[1]
      dxN[2] = xjaci[1,1]*dxNl[2]+xjaci[1,2]*dyNl[2]
      dxN[3] = xjaci[1,1]*dxNl[3]+xjaci[1,2]*dyNl[3]
      dxN[4] = xjaci[1,1]*dxNl[4]+xjaci[1,2]*dyNl[4]

      dyN[1] = xjaci[2,1]*dxNl[1]+xjaci[2,2]*dyNl[1]
      dyN[2] = xjaci[2,1]*dxNl[2]+xjaci[2,2]*dyNl[2]
      dyN[3] = xjaci[2,1]*dxNl[3]+xjaci[2,2]*dyNl[3]
      dyN[4] = xjaci[2,1]*dxNl[4]+xjaci[2,2]*dyNl[4]

      #-----------------------------
      bmat_b1  = [ 0 0 0 -dxN[1] 0
             0 0 0      0  -dyN[1]
             0 0 0 -dyN[1] -dxN[1]];

             bmat_b2  = [ 0 0 0 -dxN[2]     0
             0 0 0      0 -dyN[2]
             0 0 0 -dyN[2] -dxN[2]];

             bmat_b3  = [ 0 0 0 -dxN[3]     0
             0 0 0      0 -dyN[3]
             0 0 0 -dyN[3] -dxN[3]];

             bmat_b4  = [ 0 0 0 -dxN[4]     0
             0 0 0      0 -dyN[4]
             0 0 0 -dyN[4] -dxN[4]];


             bmat_b = [bmat_b1 bmat_b2 bmat_b3 bmat_b4];
             #-----------------------------
             bmat_m1d  = [ dxN[1]     0  0
                   0 dyN[1] 0
              dyN[1] dxN[1] 0];

              bmat_m2d  = [ dxN[2]      0 0
                   0 dyN[2] 0
              dyN[2] dxN[2] 0];

              bmat_m3d  = [ dxN[3]      0 0
                   0 dyN[3] 0
              dyN[3] dxN[3] 0];

              bmat_m4d  = [ dxN[4]      0 0
                   0 dyN[4] 0
              dyN[4] dxN[4] 0];

              bmat_mir  = [ 0 0
                            0 0
                            0 0];

              bmat_m1 = [bmat_m1d*Rot bmat_mir];
              bmat_m2 = [bmat_m2d*Rot bmat_mir];
              bmat_m3 = [bmat_m3d*Rot bmat_mir];
              bmat_m4 = [bmat_m4d*Rot bmat_mir];

              bmat_m = [bmat_m1 bmat_m2 bmat_m3 bmat_m4];

              #-----------------------------
              cx = [ 0 1 0 -1]
              cy = [-1 0 1 0 ]

              c     = zeros(8,8);

              b_bar = zeros(8,12)

              N= zeros(4,1)

              for i = 1 : 4
                N[1] = (1-cx[i])*(1-cy[i])/4 ;
                N[2] = (1+cx[i])*(1-cy[i])/4 ;
                N[3] = (1+cx[i])*(1+cy[i])/4 ;
                N[4] = (1-cx[i])*(1+cy[i])/4 ;

                dxNl[1] = (-1+cy[i])/4;
                dxNl[2] = ( 1-cy[i])/4;
                dxNl[3] = ( 1+cy[i])/4;
                dxNl[4] = (-1-cy[i])/4;

                dyNl[1] = (-1+cx[i])/4;
                dyNl[2] = (-1-cx[i])/4;
                dyNl[3] = ( 1+cx[i])/4;
                dyNl[4] = ( 1-cx[i])/4;

                xjacm[1,1] = x[1]*dxNl[1] + x[2]*dxNl[2] + x[3]*dxNl[3] + x[4]*dxNl[4];
                xjacm[1,2] = y[1]*dxNl[1] + y[2]*dxNl[2] + y[3]*dxNl[3] + y[4]*dxNl[4];
                xjacm[2,1] = x[1]*dyNl[1] + x[2]*dyNl[2] + x[3]*dyNl[3] + x[4]*dyNl[4];
                xjacm[2,2] = y[1]*dyNl[1] + y[2]*dyNl[2] + y[3]*dyNl[3] + y[4]*dyNl[4];

                jpos = [ i*2-1  i*2 ];

                c[jpos,jpos] = xjacm;

                bmat_s1  = [ dxN[1] -N[1]    0
                    dyN[1]     0 -N[1]];

                    bmat_s2  = [ dxN[2] -N[2]    0
                    dyN[2]     0 -N[2]];

                    bmat_s3  = [ dxN[3] -N[3]    0
                    dyN[3]     0 -N[3]];

                    bmat_s4  = [ dxN[4] -N[4]    0
                    dyN[4]     0 -N[4]];

                    bmat_s = [bmat_s1 bmat_s2 bmat_s3 bmat_s4];

                    b_bar[2*i-1,:] = bmat_s[1,:];
                    b_bar[2*i,:] = bmat_s[2,:]
                end
                #-----------------------------
                    T_mat = [ 1  0  0  0  0  0  0  0
                    0  0  0  1  0  0  0  0
                    0  0  0  0  1  0  0  0
                    0  0  0  0  0  0  0  1 ];

                    P_mat = [ 1  -1   0   0
                    0   0   1   1
                    1   1   0   0
                    0   0   1  -1 ];

                    A_mat = [ 1  ygs  0    0
                    0    0  1  xgs ];

                    bmat_ss = xjaci * A_mat * inv(P_mat) * T_mat * c * b_bar

                    bmat_s1 = [0  0 bmat_ss[1, 1]
                    0  0 bmat_ss[2, 1]];

                    bmat_s2 = [0  0 bmat_ss[1, 4]
                    0  0 bmat_ss[2, 4]];

                    bmat_s3 = [0  0 bmat_ss[1, 7]
                    0  0 bmat_ss[2, 7]];

                    bmat_s4 = [0  0 bmat_ss[1,10]
                    0  0 bmat_ss[2,10]];

                    bmat_s1 = [bmat_s1*Rot bmat_ss[:,2:3]];

                    bmat_s2 = [bmat_s2*Rot bmat_ss[:,5:6]];

                    bmat_s3 = [bmat_s3*Rot bmat_ss[:,8:9]];

                    bmat_s4 = [bmat_s4*Rot bmat_ss[:,11:12]];

                    bmat_s = [bmat_s1 bmat_s2 bmat_s3 bmat_s4]

                    #-----------------------------
                    K_b = bmat_b'*D_matb*bmat_b*area*gauss_w[igaus];
                    K_m = bmat_m'*D_matm*bmat_m*area*gauss_w[igaus];
                    K_s = bmat_s'*D_mats*bmat_s*area*gauss_w[igaus];

                    K_elem += K_b + K_m + K_s


            end
        map = elem_map(elem)

    return K_elem, map, map
end


function elem_update!(elem::ShellQuad4node, U::Array{Float64,1}, F::Array{Float64,1}, dt::Float64)
    K, map, map = elem_stiffness(elem)
    dU  = U[map]
    F[map] += K*dU
end
