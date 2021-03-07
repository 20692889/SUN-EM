function [Solver_setup] = AddRWGatObservationPt(Const, Solver_setup)

 %AddRWGatObservationPt
    %   Usage:
    %        [Solver_setup] = AddRWGatObservationPt(Const, Solver_setup)

    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       x, y, z
    %           
    %       Solver_setup (struct)
    %           Struct containing the frequency range, nodal co-ordinates of the 
    %           triangles, as well as information on the shared edges, RWG basis functions, 
    %           etc. Also contains the observation point vectors.
    %    
    %   Output Arguments:
    %       Solver_setup (struct)
    %           The updated struct Solver_setup containing the RWG
    %           information at the observation points
 
    lambda = Const.C0/Solver_setup.frequencies.samples(Solver_setup.frequencies.freq_num);
    lp = lambda/10; %calculate edge_length of RWG at observation point
    Solver_setup.lengthP = lp;
    rc_obs = Solver_setup.r_c_obs;
    %vector setup of x-directed RWG element (used to ultimately calculate rho vectors for observation RWGs)
    vert1_plusX = [rc_obs(:,1), rc_obs(:,2) - lp/2, rc_obs(:,3)];
    vert2_plusX = [rc_obs(:,1), rc_obs(:,2) + lp/2, rc_obs(:,3)];
    vert3_plusX = [rc_obs(:,1) + (1.5*lp), rc_obs(:,2), rc_obs(:,3)];
    
    Solver_setup.r_ct_plusx = (1/3) .* (vert1_plusX + vert2_plusX + vert3_plusX);
    
    vert1_minX = [rc_obs(:,1), rc_obs(:,2) - lp/2, rc_obs(:,3)];
    vert2_minX = [rc_obs(:,1), rc_obs(:,2) + lp/2, rc_obs(:,3)];
    vert3_minX = [rc_obs(:,1) - (1.5*lp), rc_obs(:,2), rc_obs(:,3)];
    
    Solver_setup.r_ct_minx = (1/3) .* (vert1_minX + vert2_minX + vert3_minX);
    
    Solver_setup.rho_pc_plusx = Solver_setup.r_ct_plusx - vert3_plusX; %rho vector directed from vertex
    Solver_setup.rho_pc_minx = -(Solver_setup.r_ct_minx - vert3_minX);%rho vector directed to vertex  
    
    %y-direction 
    vert1_plusY = [rc_obs(:,1)- lp/2, rc_obs(:,2), rc_obs(:,3)];
    vert2_plusY = [rc_obs(:,1)+ lp/2, rc_obs(:,2) , rc_obs(:,3)];
    vert3_plusY = [rc_obs(:,1) , rc_obs(:,2)+ (1.5*lp), rc_obs(:,3)];
    
    Solver_setup.r_ct_plusy = (1/3) .* (vert1_plusY + vert2_plusY + vert3_plusY);
    
    vert1_minY = [rc_obs(:,1) - lp/2, rc_obs(:,2), rc_obs(:,3)];
    vert2_minY = [rc_obs(:,1)+ lp/2, rc_obs(:,2) , rc_obs(:,3)];
    vert3_minY = [rc_obs(:,1) , rc_obs(:,2)- (1.5*lp), rc_obs(:,3)];
    
    Solver_setup.r_ct_miny = (1/3) .* (vert1_minY + vert2_minY + vert3_minY);
    
    Solver_setup.rho_pc_plusy = Solver_setup.r_ct_plusy - vert3_plusY; %rho vector directed from vertex
    Solver_setup.rho_pc_miny = -(Solver_setup.r_ct_miny - vert3_minY);%rho vector directed to vertex
    
    %z-direction
    vert1_plusZ = [rc_obs(:,1) - lp/2, rc_obs(:,2), rc_obs(:,3)];
    vert2_plusZ = [rc_obs(:,1) + lp/2, rc_obs(:,2), rc_obs(:,3)];
    vert3_plusZ = [rc_obs(:,1) , rc_obs(:,2), rc_obs(:,3)+ (1.5*lp)];
    
    Solver_setup.r_ct_plusz = (1/3) .* (vert1_plusZ + vert2_plusZ + vert3_plusZ);
    
    vert1_minZ = [rc_obs(:,1)- lp/2, rc_obs(:,2) , rc_obs(:,3)];
    vert2_minZ = [rc_obs(:,1)+ lp/2, rc_obs(:,2) , rc_obs(:,3)];
    vert3_minZ = [rc_obs(:,1), rc_obs(:,2), rc_obs(:,3) - (1.5*lp)];
    
    Solver_setup.r_ct_minz = (1/3) .* (vert1_minZ + vert2_minZ + vert3_minZ);
    
    Solver_setup.rho_pc_plusz = Solver_setup.r_ct_plusz - vert3_plusZ; %rho vector directed from vertex
    Solver_setup.rho_pc_minz = -(Solver_setup.r_ct_minz - vert3_minZ);%rho vector directed to vertex
    
    %Solver_setup.r_ct_plusx = [rc_obs(:,1)-(1/3)*(lp/2)*tand(60), rc_obs(:,2), rc_obs(:,3)];
   % Solver_setup.r_ct_minx = [rc_obs(:,1)+(1/3)*(lp/2)*tand(60), rc_obs(:,2), rc_obs(:,3)];
   % fvex_plusx = [rc_obs(:,1)-(lp/2)*tand(60), rc_obs(:,2), rc_obs(:,3)]; %free vertex Tp+ oriented with free vertices on x-axis
   % fvex_minx = [rc_obs(:,1)+(lp/2)*tand(60), rc_obs(:,2), rc_obs(:,3)]; %free vertex Tp- oriented with free vertices on x-axis
    
  %  Solver_setup.rho_pc_plusx = Solver_setup.r_ct_plusx - fvex_plusx; %rho vector directed from vertex
   % Solver_setup.rho_pc_minx = -(Solver_setup.r_ct_minx - fvex_minx);%rho vector directed to vertex
    
    %vector setup of y-directed RWG element (used to ultimately calculate rho vectors for observation RWGs)
   % Solver_setup.r_ct_plusy = [rc_obs(:,1), rc_obs(:,2)-(1/3)*(lp/2)*tand(60), rc_obs(:,3)];
   % Solver_setup.r_ct_miny = [rc_obs(:,1), rc_obs(:,2)+(1/3)*(lp/2)*tand(60), rc_obs(:,3)];
   % fvex_plusy = [rc_obs(:,1), rc_obs(:,2)-(lp/2)*tand(60), rc_obs(:,3)];%free vertex Tp+ oriented with free vertices on y-axis
   % fvex_miny = [rc_obs(:,1), rc_obs(:,2)+(lp/2)*tand(60), rc_obs(:,3)];%free vertex Tp- oriented with free vertices on y-axis
    
   % Solver_setup.rho_pc_plusy = Solver_setup.r_ct_plusy - fvex_plusy; %rho vector directed from vertex
   % Solver_setup.rho_pc_miny = -(Solver_setup.r_ct_miny - fvex_miny);%rho vector directed to vertex
    
    %vector setup of z-directed RWG element (used to ultimately calculate rho vectors for observation RWGs)
   % Solver_setup.r_ct_plusz = [rc_obs(:,1), rc_obs(:,2), rc_obs(:,3)-(1/3)*(lp/2)*tand(60)];
   % Solver_setup.r_ct_minz = [rc_obs(:,1), rc_obs(:,2), rc_obs(:,3)+(1/3)*(lp/2)*tand(60)];
   % fvex_plusz = [rc_obs(:,1), rc_obs(:,2), rc_obs(:,3)-(lp/2)*tand(60)];%free vertex Tp+ oriented with free vertices on z-axis
   % fvex_minz = [rc_obs(:,1), rc_obs(:,2), rc_obs(:,3)+(lp/2)*tand(60)];%free vertex Tp- oriented with free vertices on z-axis
    
   % Solver_setup.rho_pc_plusz = Solver_setup.r_ct_plusz - fvex_plusz;%rho vector directed from vertex
   % Solver_setup.rho_pc_minz = -(Solver_setup.r_ct_minz - fvex_minz);%rho vector directed to vertex
    

    end