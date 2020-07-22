# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct TMShell<:Thermomechanical
    id    ::Int
    shape ::ShapeType
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    active::Bool
    linked_elems::Array{Element,1}
    env::ModelEnv

    function TMShell();
        return new()
    end
end

matching_shape_family(::Type{TMShell}) = SOLID_SHAPE

#=
function elem_config_dofs(elem::TMShell)
    nbnodes = elem.shape.basic_shape.npoints
    for (i, node) in enumerate(elem.nodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.env.ndim==3 && add_dof(node, :uz, :fz)
        if  i<=(nbnodes)
            add_dof(node, :ut, :ft)
        end
    end
end
=#

function elem_config_dofs(elem::TMShell)
    ndim = elem.env.ndim
    ndim == 1 && error("ShellQUAD4: Shell elements do not work in 1d analyses")
    #if ndim==2
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :ut, :ft)
        end
    #else
        #error("ShellQUAD4: Shell elements do not work in this analyses")
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
    #end
end


function elem_init(elem::TMShell)
    nothing
end


function distributed_bc(elem::TMShell, facet::Union{Facet,Nothing}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.env.ndim
    th    = elem.env.thickness
    suitable_keys = (:tx, :ty, :tz, :tn, :tq)

    # Check keys
    key in suitable_keys || error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")
    (key == :tz && ndim==2) && error("distributed_bc: boundary condition $key is not applicable in a 2D analysis")

    target = facet!=nothing ? facet : elem
    nodes  = target.nodes
    nnodes = length(nodes)
    t      = elem.env.t

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = get_coords(nodes, ndim)

    # Vector with values to apply
    Q = zeros(ndim)

    # Calculate the nodal values
    F     = zeros(nnodes, ndim)
    shape = target.shape
    ips   = get_ip_coords(shape)

    if key == :tq # energy per area
        for i=1:size(ips,1)
            R = vec(ips[i,:])
            w = R[end]
            N = shape.func(R)
            D = shape.deriv(R)
            J = D*C
            nJ = norm2(J)
            X = C'*N
            if ndim==2
                x, y = X
                vip = eval_arith_expr(val, t=t, x=x, y=y)
            else
                x, y, z = X
                vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            end
            coef = vip*nJ*w
            F .+= N*coef # F is a vector
        end

        # generate a map
        map  = [ node.dofdict[:ut].eq_id for node in target.nodes ]

        return F, map
    end

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
        X = C'*N
        if ndim==2
            x, y = X
            vip = eval_arith_expr(val, t=t, x=x, y=y)
            if key == :tx
                Q = [vip, 0.0]
            elseif key == :ty
                Q = [0.0, vip]
            elseif key == :tn
                n = [J[1,2], -J[1,1]]
                Q = vip*normalize(n)
            end
            if elem.env.modeltype=="axisymmetric"
                th = 2*pi*X[1]
            end
        else
            x, y, z = X
            vip = eval_arith_expr(val, t=t, x=x, y=y, z=z)
            if key == :tx
                Q = [vip, 0.0, 0.0]
            elseif key == :ty
                Q = [0.0, vip, 0.0]
            elseif key == :tz
                Q = [0.0, 0.0, vip]
            elseif key == :tn && ndim==3
                n = cross(J[1,:], J[2,:])
                Q = vip*normalize(n)
            end
        end
        coef = norm2(J)*w*th
        @gemm F += coef*N*Q' # F is a matrix
    end

    # generate a map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in target.nodes for key in keys ]

    return reshape(F', nnodes*ndim), map
end

#=
@inline function set_Bu(elem::Element, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B) # using function setB from mechanical analysis
end
=#

# the strain-displacement matrix for membrane forces
function D_matrixm(elem::TMShell)

    coef1 = elem.mat.thick*elem.mat.E/(1-elem.mat.nu^2)
    coef2 = elem.mat.nu*coef1
    coef3 = coef1*(1-elem.mat.nu)/2

        D_matm = [coef1  coef2 0
                  coef2  coef1 0
                  0      0     coef3];
    return D_matm
end

# the strain-displacement matrix for bending moments
function D_matrixb(elem::TMShell)

    D_matm = D_matrixm(elem)

    D_matb = D_matm*(elem.mat.thick^2/12)

    return D_matb
end

# the strain-displacement matrix for shear forces

function D_matrixs(elem::TMShell)

    coef = elem.mat.thick*(5/6)*elem.mat.E/(2*(1+elem.mat.nu))

            D_mats = [coef    0
                        0     coef];
    return D_mats
end

# Rotation Matrix

function RotMatrix(elem::TMShell)

    v12 = zeros(3,1)
    v13 = zeros(3,1)
    vxe = zeros(3,1)
    vye = zeros(3,1)
    vze = zeros(3,1)

    cxyz = get_coords(elem)

        if size(cxyz,2)==2
            v12[3] = 0
            v13[3] = 0
        else
            v12[3] = cxyz[2,3] - cxyz[1,3]
            v13[3] = cxyz[3,3] - cxyz[1,3]
        end

    v12[1] = cxyz[2,1] - cxyz[1,1]
    v12[2] = cxyz[2,2] - cxyz[1,2]

    v13[1] = cxyz[3,1] - cxyz[1,1]
    v13[2] = cxyz[3,2] - cxyz[1,2]

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

function elem_map(elem::TMShell)::Array{Int,1}

    #if elem.env.ndim==2
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry)
    #else
    #    dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz) # VERIFICAR
    #end

    dof_keys = (:ux, :uy, :uz, :rx, :ry)

    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)

end

function elem_stiffness(elem::TMShell)

    nnodes = length(elem.nodes)

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

    C = get_coords(elem)

        if size(C,2)==2
            cxyz  = zeros(4,3)
            cxyz[:,1:2]  = get_coords(elem)
        else

            cxyz = C
        end


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

                jpos = [i*2-1  i*2];

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


# matrix C
function elem_coupling_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = get_coords(elem)
    Bu  = zeros(6, nnodes*ndim)
    Cut = zeros(nnodes*ndim, nbnodes) # u-t coupling matrix

    J    = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    m    = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]
    β    = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu) # thermal stress modulus

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        # compute Cut
        Nt    = elem.shape.basic_shape.func(ip.R)
        coef  = β
        coef *= detJ*ip.w*th
        mNt   = m*Nt'
        @gemm Cut -= coef*Bu'*mNt
    end
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return Cut, map_u, mat_t
end

# thermal conductivity
function elem_conductivity_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)
    H      = zeros(nnodes, nnodes)
    Bt     = zeros(ndim, nnodes)
    KBt    = zeros(ndim, nnodes)
    J    = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        dNdR  = elem.shape.deriv(ip.R)
        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm J  = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")
        @gemm Bt = inv(J)*dNtdR

        # compute H
        K = calcK(elem.mat, ip.state)
        coef = detJ*ip.w*th
        @gemm KBt = K*Bt
        @gemm H  -= coef*Bt'*KBt
    end

    # map
    map = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    return H, map, map
end

function elem_mass_matrix(elem::TMShell)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)
    M      = zeros(nnodes, nnodes)

    J  = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        Nt   = elem.shape.basic_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute Cut
        coef  = elem.mat.ρ*elem.mat.cv
        coef *= detJ*ip.w*th
        M    -= coef*Nt*Nt'
    end

    # map
    map = [  node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes]  ]

    return M, map, map
end

#=
function elem_internal_forces(elem::TMShell, F::Array{Float64,1}, DU::Array{Float64,1})
    ndim   = elem.env.ndim
    th     = elem.env.thickness # VERIFICAR ESPESSURA
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C   = get_coords(elem)
    T0     = elem.env.T0 + 273.15

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nbnodes)
    Bt  = zeros(ndim, nbnodes)

    m = [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ] # = tI

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, ndim, nnodes)
    Jp  = Array{Float64}(undef, ndim, nbnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    dUt = DU[mat_t] # nodal temperature increments
    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu matrix and Bt
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        @gemm dNdX = inv(J)*dNdR
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        Jp = dNtdR*Ct
        @gemm dNtdX = inv(Jp)*dNtdR
        Bt = dNtdX
        # compute N

        # internal force
        ut   = ip.state.ut + 273
        β   = elem.mat.E*elem.mat.α/(1-2*elem.mat.nu)
        σ    = ip.state.σ - β*ut*m # get total stress
        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*σ

        # internal volumes dFt
        ε    = ip.state.ε
        εvol = dot(m, ε)
        coef = β*detJ*ip.w*th
        dFt  -= coef*Nt*εvol

        coef = detJ*ip.w*elem.mat.ρ*elem.mat.cv*th/T0
        dFt -= coef*Nt*ut

        QQ   = ip.state.QQ
        coef = detJ*ip.w*th/T0
        @gemv dFt += coef*Bt'*QQ
    end

    F[map_u] += dF
    F[mat_t] += dFt
end
=#


function elem_update!(elem::TMShell, DU::Array{Float64,1}, DF::Array{Float64,1}, Δt::Float64)
    ndim   = elem.env.ndim
    th     = elem.env.thickness
    T0     = elem.env.T0 + 273.15
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.basic_shape.npoints
    C      = get_coords(elem)

    E = elem.mat.E
    α = elem.mat.α
    ρ = elem.mat.ρ
    nu = elem.mat.nu
    cv = elem.mat.cv
    β = E*α/(1-2*nu)

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    mat_t  = [ node.dofdict[:ut].eq_id for node in elem.nodes[1:nbnodes] ]

    dU  = DU[map_u] # nodal displacement increments
    dUt = DU[mat_t] # nodal temperature increments
    Ut  = [ node.dofdict[:ut].vals[:ut] for node in elem.nodes]
    Ut += dUt # nodal tempeture at step n+1
    m   = tI  # [ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 ]  #

    dF  = zeros(nnodes*ndim)
    Bu  = zeros(6, nnodes*ndim)
    dFt = zeros(nnodes)
    Bt  = zeros(ndim, nnodes)

    J     = Array{Float64}(undef, ndim, ndim)
    dNdX  = Array{Float64}(undef, ndim, nnodes)
    dNtdX = Array{Float64}(undef, ndim, nbnodes)
    Δε = zeros(6)

    for ip in elem.ips
        elem.env.modeltype=="axisymmetric" && (th = 2*pi*ip.coord.x)

        # compute Bu and Bt matrices
        dNdR = elem.shape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        invJ = inv(J)
        @gemm dNdX = invJ*dNdR
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.basic_shape.deriv(ip.R)
        @gemm dNtdX = invJ*dNtdR

        # compute Nt
        Nt = elem.shape.basic_shape.func(ip.R)
        # compute Δε

        @gemv Δε = Bu*dU

        # compute Δut
        Δut = Nt'*dUt # interpolation to the integ. point

        # compute thermal gradient G
        Bt = dNtdX
        G  = Bt*Ut

        # internal force dF
        Δσ, q = stress_update(elem.mat, ip.state, Δε, Δut, G, Δt)
        Δσ -= β*Δut*m # get total stress

        coef = detJ*ip.w*th
        @gemv dF += coef*Bu'*Δσ

        # internal volumes dFt
        Δεvol = dot(m, Δε)
        coef  = β*Δεvol*T0
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt

        coef  = ρ*cv
        coef *= detJ*ip.w*th
        dFt  -= coef*Nt*Δut

        coef  = Δt
        coef *= detJ*ip.w*th
        @gemv dFt += coef*Bt'*q
    end

    DF[map_u] += dF
    DF[mat_t] += dFt
end
