function EfieldatpointSpherical = calculateEfieldEFIE(Const, r, theta_degrees, phi_degrees, ...
    Solver_setup, Isol)

  %   Usage:
    %       [EfieldAtPointSpherical] =  calculateEfieldAtPointRWG(r, theta, phi, nodes, ...
    %           sharedEdgesList, triangleDataList, Isol, k )
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       r, theta_degrees, phi_degrees
    %           Spherical co-ordinate of the observation point. Note the angles are in degrees.
    %       Solver_setup (struct)
    %           Struct containing the frequency range, nodal co-ordinates of the 
    %           triangles, as well as information on the shared edges, RWG basis functions, 
    %           etc. 
    %       Isol
    %           The solution vector, i.e. the MoM expansion coefficients
    %   Output Arguments:
    %       EfieldAtPointSpherical
    %           The E-field strength calculated at the point (r,theta,phi)
    %           EfieldAtPoint = A r^ + B theta^ + C phi^
    %
    %   Description:
    %       Calculates the E-field value at a certain point
    %       (r,theta,phi). Makes use of surface integrals and RWG elements
    EfieldAtPointCartesian = zeros(Solver_setup.frequencies.freq_num, 3);
    EfieldAtPointSpherical = zeros(Solver_setup.frequencies.freq_num, 3);
    % Populate now the variables as used by [1] - esp. replace the use of
    % global variables.
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    elements = Solver_setup.triangle_vertices;        % Replacing global ELEMENTS
    node_coord = Solver_setup.nodes_xyz;              % Replacing global NODE_COORD
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    
    quad_pts = Const.QUAD_PTS;
    sing     = Const.SING;
    eps_0    = Const.EPS_0;
    mu_0     = Const.MU_0;
    
    syms x y z 

    % Extract the triangle midpoints
    r_c = Solver_setup.triangle_centre_point;
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;
    
    [Px,Py,Pz] = transformSpericalCoordinateToCartesian(r,theta_degrees*Const.DEG2RAD,phi_degrees*Const.DEG2RAD);
    
    for freq_index = 1:Solver_setup.frequencies.freq_num

        freq = Solver_setup.frequencies.samples(freq_index);
        % Calculate some frequency dependent parameters required below
        omega = 2*pi*freq;       % Radial frequency
        lambda = Const.C0/freq;  % Wavelength in m
        k  = 2*pi/lambda;        % Wavenumber in rad/m
        jk = 1i*k;
       

        % The E-field at the point is calculated as the superposition of the individual 
        % RWG contributions.
        for mm = 1:num_dofs
            
            %pp_pls = EDGECONXELEMS(mm,1);            
            field_pt = Solver_setup.rwg_basis_functions_trianglePlus(mm);
            %pp_mns = EDGECONXELEMS(mm,2);
           % pp_mns = Solver_setup.rwg_basis_functions_triangleMinus(mm);
          
            for nn = 1:num_dofs
                
                  %qq_pls = EDGECONXELEMS(nn,1);
                   qq_pls = Solver_setup.rwg_basis_functions_trianglePlus(nn);
                  %qq_mns = EDGECONXELEMS(nn,2);
                   qq_mns = Solver_setup.rwg_basis_functions_triangleMinus(nn);
                
                   triangle_tn_plus_free_vertex = Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(nn);
                   triangle_tn_minus_free_vertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(nn);
            
                   [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,field_pt,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    
                    Amn_pls_source_pls = -MagVecPot;
                    Phi_mn_pls_source_pls = +ScalPot;
                    
                    
                   [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,field_pt,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    
                    Amn_pls_source_mns = + MagVecPot;
                    Phi_mn_pls_source_mns = -ScalPot;
                                    
                    Amn = Amn_pls_source_pls + Amn_pls_source_mns;
                    Phi_mn = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
                    
                   % Phix = gradient(Phi_mn);
                    
                    Hfield = curl(Amn, [x,y,z]);
                   
                    Efield_x = - 1i*omega*Amn(1);
                    Efield_y = -1i*omega*Amn(2);
                    Efield_z = -1i*omega*Amn(3);
                    
                    
                    EfieldAtPointCartesian(freq_index, 1) = EfieldAtPointCartesian(freq_index, 1) + Efield_x;
                    EfieldAtPointCartesian(freq_index, 2) = EfieldAtPointCartesian(freq_index, 2) + Efield_y;
                    EfieldAtPointCartesian(freq_index, 3) = EfieldAtPointCartesian(freq_index, 3) + Efield_z;
                    
                    
                    
            
        end
        end
    
    end
    EfieldAtPointSpherical(freq_index,:) = transformCartesianVectorToSpherical(EfieldAtPointCartesian(freq_index,:), Px, Py, Pz);
end
            
            


function [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, field_pt,source_pt,field_edge, source_edge, ii, k,r_c,quad_pts,sing,eps_0,mu_0,omega)
        
    % This subfunction computes the magnetic vector and scalar potentials for 
    % field point face field_pt and source point face source_pt. 

    % Code for debugging singularity scheme below:
    %[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq_debug(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    
    % Extract the nodes of the source triangle (q)
    qnodes = elements(source_pt,:);    
    r = zeros(3,3); % Store all three vertices of triangle q here in r()
    r(1,:) = node_coord(qnodes(1),:);
    r(2,:) = node_coord(qnodes(2),:);
    r(3,:) = node_coord(qnodes(3),:);
    
    % Extract the position of the ith free vertex (i.e. associated with the
    % ith RWG)
    %ii_nodes = elements(ii,:);
    rii = zeros(1,3);
    rii(1,1) = node_coord(ii,1);
    rii(1,2) = node_coord(ii,2);
    rii(1,3) = node_coord(ii,3);
        
    
    %ii = DOFLOCALNUM(source_edge,source_tri); % This is the free vertex
    %associated with the source_edge - which is now passed here as an
    %argument by the calling routine.    

    % [RWG82, Eq. (32) - without sign
    MagVecPot = mu_0*ell(source_edge)/(4*pi)*...
        ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - rii(1,:)*Ipq);

    % [RWG82, Eq. (33) - without sign
    ScalPot = ell(source_edge)/(1i*2*pi*omega*eps_0) * Ipq;   
        
        
        
        
        
end


function vectorAspherical = transformCartesianVectorToSpherical(vectorAcartesian, x, y, z)
    %vectorAspherical
    %   Usage:
    %           vectorAspherical = transformCartesianVectorToSpherical(vectorAcartesian, x, y, z)
    %
    %   Input Arguments:
    %       vectorAcartesian
    %           The x, y and z components of the spatial vector A = [Ax, Ay, Az]
    %       x,y,z
    %           The cartesian point where the vector is defined
    %
    %   Output Arguments:
    %       vectorAspherical
    %           The spherical components of the vector A = [Ar Atheta Aphi]
    %
    %   Description:
    %       Transforms A vector given in Cartesian Coordinates, to one
    %       in Spherical Coordinates
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   =======================

    r = sqrt(x^2 + y^2 + z^2);
    theta = acos(x/r);
    phi = atan(y/x); 

    Ax = vectorAcartesian(1);
    Ay = vectorAcartesian(2);
    Az = vectorAcartesian(3);

    Ar = Ax*sin(theta)*cos(phi) + Ay*sin(theta)*sin(phi) + Az*cos(theta);
    Atheta = Ax*cos(theta)*cos(phi) + Ay*cos(theta)*sin(phi) - Az*sin(theta);
    Aphi = -Ax*sin(phi) + Ay*cos(phi);

    vectorAspherical = [Ar Atheta Aphi];
    
end