% Author: Danie Ludick (dludick@sun.ac.za)
% Project: Strip dipole antenna simulation
%
% Note: Each project directory / example directory needs to have a sunem_initialise.m
% script, that is used to setup the correct environment for that example.
%
% Refer to the /doc folder for more information

% --------------------------------------------------------------------------------------------------
% Initialise the environment
% --------------------------------------------------------------------------------------------------
% Project output directory: './dipoles/'
% Debug: True/False
Const = sunem_initialise('dipoles',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings--------------------------------------------------------------------------------------------------
% 

% Choose the solvers that will be executed
Const.runMoMsolver              = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'strip_dipole.mat'; %Z-matrix by FEKO
Const.FEKOstrfilename          = 'strip_dipole.str'; %I-vector by FEKO
Const.FEKOrhsfilename          = 'strip_dipole.rhs'; %V-vector by FEKO
Const.FEKOoutfilename          = 'strip_dipole.out';
Const.FEKOefefilename          = 'strip_dipole.efe'; %electric nearfield by FEKO
Const.FEKOffefilename          = 'strip_dipole.ffe'; %electric farfield by FEKO
Const.FEKOhfefilename          = 'strip_dipole(1).hfe';%magnetic nearfield by FEKO

% The Following file is used to port solutions to FEKO 
% (for post-processing in POSTFEKO).
% TO-DO: [DL] Add this.
% Const.output_strfilename    = '';
% Const.writeFEKOstrfile = [0 0 0 0];


% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
% Const.useACA = 1;
% %Observation grid cuts
% z=100;
% y_grid=10:1:39;
% x_grid=100:1:100;
% num_y_samples = length(y_grid);
% x_rep = repelem(x_grid, num_y_samples).';
% z_rep = repelem(z, num_y_samples).';
% %vector of observation points
% Solver_setup.r_c_obs = [x_rep, y_grid.', z_rep];
% %function to add RWG elements at observation point
% [Solver_setup] = AddRWGatObservationPt(Const, Solver_setup);
% %extract SUNEM mom matrix (for now it calculates the Zmn matrix between
% %antenna surface RWGs and observation point RWGs whos vertices all lie on the
% %y-axis
% tic
% [Const, zMatricesSUNEMx, zMatricesSUNEMy, zMatricesSUNEMz, yVectorsSUNEM] = extractSUNEMMoMmatrixEq(Const, Solver_setup);
% TotZmnTime = toc;
% % --------------------------------------------------------------------------------------------------
% % Run the EM solver
% % --------------------------------------------------------------------------------------------------
% SourceVecInd = 1:Solver_setup.num_metallic_edges;
% ObservVecInd = 1:length(Solver_setup.r_c_obs);
% tic
% [Zmn,Ux,Vx] = calcZmn(Const,1, Solver_setup, zMatricesSUNEMx, Solver_setup.frequencies.samples, zMatricesSUNEMx.mBasis, zMatricesSUNEMx.nBasis, ObservVecInd, SourceVecInd);
% [Zmn,Uy,Vy] = calcZmn(Const,2, Solver_setup, zMatricesSUNEMy, Solver_setup.frequencies.samples, zMatricesSUNEMy.mBasis, zMatricesSUNEMy.nBasis, ObservVecInd, SourceVecInd);
% [Zmn,Uz,Vz] = calcZmn(Const,3, Solver_setup, zMatricesSUNEMz, Solver_setup.frequencies.samples, zMatricesSUNEMz.mBasis, zMatricesSUNEMz.nBasis, ObservVecInd, SourceVecInd);
% TotUVtime = toc; 
% Const.useACA = 0;
% [Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
% 
% % --------------------------------------------------------------------------------------------------
% % Postprocess the results, e.g. calculate the Electric field
% % --------------------------------------------------------------------------------------------------
% 
% %[Uvector, Vvector] = AdapadtiveCrossApproximation(zMatricesSUNEMy.values, zMatricesSUNEMy.mBasis, zMatricesSUNEMy.nBasis);
% %
% %[Uapprox, Vapprox] = AdapadtiveCrossApproximation(zMatricesSUNEMx.values, zMatricesSUNEMx.mBasis, zMatricesSUNEMx.nBasis);
% %Eqn (4) in "Field Computations Through the ACA Algorithm" by Maaskant, R. ; Lancelotti, V. (2015) 
% 
% lp = Solver_setup.lengthP;
% tic
% ExACA = -(1/(lp*lp))*Ux*(Vx*Solution.mom.Isol);
% CalcExACAtime = toc;
% tic
% EyACA = -(1/(lp*lp))*Uy*(Vy*Solution.mom.Isol);
% CalcEyACAtime = toc;
% tic
% EzACA = -(1/(lp*lp))*Uz*(Vz*Solution.mom.Isol);
% CalcEzACAtime = toc;
% TotTimeACA = CalcExACAtime + CalcEyACAtime + CalcEzACAtime;
% tic
% Ex = -(1/(lp*lp))*zMatricesSUNEMx.values*Solution.mom.Isol;
% CalcExMaastime = toc;
% tic
% Ey = -(1/(lp*lp))*zMatricesSUNEMy.values*Solution.mom.Isol;
% CalcEyMaastime = toc;
% tic
% Ez = -(1/(lp*lp))*zMatricesSUNEMz.values*Solution.mom.Isol;
% CalcEzMaastime = toc;
% TotTimeMaas = CalcExMaastime + CalcEyMaastime + CalcEzMaastime;
% %(1/(lp*lp))*
% Efield_ACA = [ExACA, EyACA, EzACA];
% Efield_Maaskant = [Ex, Ey, Ez];
% 
% 
% r = sqrt((x_rep).^2 + (y_grid.').^2 + (z_rep).^2);
% 
% 
% %r = 100;%100;
% Efield_ACA_mag = sqrt(abs(r.*Efield_ACA(:,1)).^2 + ...
%             abs(r.*Efield_ACA(:,2)).^2 + ...
%             abs(r.*Efield_ACA(:,3)).^2);
% Efield_Maaskant_mag= sqrt(abs(r.*Efield_Maaskant(:,1)).^2 + ...
%             abs(r.*Efield_Maaskant(:,2)).^2 + ...
%             abs(r.*Efield_Maaskant(:,3)).^2);
% % Loop over a few theta and phi points and compare the results with that of FEKO
% %theta_grid = 90:1:90;
% %phi_grid = 1:1:91;
% 
% num_x_samples = length(x_grid);
% num_y_samples = length(y_grid);
% num_z_samples = length(z);
% %num_theta_samples = length(theta_grid);
% %num_phi_samples = length(phi_grid);
% total_efield_samples = num_x_samples*num_y_samples*num_z_samples;
% Efield_magnitude = zeros(total_efield_samples,1);
% Hfield_magnitude = zeros(total_efield_samples,1);
% 
% difference = zeros(total_efield_samples,3);
% 
% Hfield_vectors = zeros(total_efield_samples,3);
% Efield_vectors = zeros(total_efield_samples,3);
% FEKO_Efield_vectors = zeros(total_efield_samples,3);
% FEKO_Hfield_vectors = zeros(total_efield_samples,3);
% 
% %Read the FEKO data from a *.hfe file for comparison
% FEKO_HfieldAtPoint = parseFEKOhfefile(Const, Const.FEKOhfefilename);
% FEKO_total_hfield_samples = FEKO_HfieldAtPoint.number_x_samples * ...
%     FEKO_HfieldAtPoint.number_y_samples * FEKO_HfieldAtPoint.number_z_samples;
% FEKO_Hfield_magnitude = sqrt(abs(r.*FEKO_HfieldAtPoint.Hx).^2 + ...
%     abs(r.*FEKO_HfieldAtPoint.Hy).^2 + ...
%     abs(r.*FEKO_HfieldAtPoint.Hz).^2);
% 
% FEKO_Hfield_vectors(:, 1) = FEKO_HfieldAtPoint.Hx;
% FEKO_Hfield_vectors(:, 2) = FEKO_HfieldAtPoint.Hy;
% FEKO_Hfield_vectors(:, 3) = FEKO_HfieldAtPoint.Hz;
% 
% % Read the FEKO data from a *.efe file for comparison
% FEKO_EfieldAtPoint = parseFEKOefefile(Const, Const.FEKOefefilename);
% FEKO_total_efield_samples = FEKO_EfieldAtPoint.number_x_samples * ...
%     FEKO_EfieldAtPoint.number_y_samples * FEKO_EfieldAtPoint.number_z_samples;
% FEKO_Efield_magnitude = sqrt(abs(r.*FEKO_EfieldAtPoint.Ex).^2 + ...
%     abs(r.*FEKO_EfieldAtPoint.Ey).^2 + ...
%     abs(r.*FEKO_EfieldAtPoint.Ez).^2);
% 
% %vectors with Efield values from FEKO
% FEKO_Efield_vectors(:, 1) = FEKO_EfieldAtPoint.Ex;
% FEKO_Efield_vectors(:, 2) = FEKO_EfieldAtPoint.Ey;
% FEKO_Efield_vectors(:, 3) = FEKO_EfieldAtPoint.Ez;
% 
% %[x, y, z] = transformSpericalCoordinateToCartesian(FEKO_EfieldAtPoint.Er(1),FEKO_EfieldAtPoint.Etheta(1),FEKO_EfieldAtPoint.Ephi(1));
% 
% %reading of far field values done by setting .efe in far field
% % Read also now FEKO's far field values here. 
% %FEKO_farfield = parseFEKOffefile(Const, Const.FEKOffefilename);
% %FEKO_total_farfield_samples =  FEKO_farfield.number_theta_samples * FEKO_farfield.number_phi_samples;
% %FEKO_farfield_magnitude = sqrt(abs(FEKO_farfield.Etheta).^2 + abs(FEKO_farfield.Ephi).^2);
% %FEKO_eff_etheta = FEKO_farfield.Etheta;
% %FEKO_eff_Ephi = FEKO_farfield.Ephi;
% %FEKO_Efield_vectors(:, 2) = FEKO_farfield.Etheta;
% %FEKO_Efield_vectors(:, 3) = FEKO_farfield.Ephi;
% %FEKO_Hfield_vectors(:, 2) = -(1/376.73)*FEKO_farfield.Ephi;
% %FEKO_Hfield_vectors(:, 3) = (1/376.73)*FEKO_farfield.Etheta;
% %FEKO_hff_theta = -(1/376.73)*FEKO_farfield.Ephi;
% %FEKO_hff_phi = (1/376.73)*FEKO_farfield.Etheta;
% %FEKO_farfield_magnitude_H = sqrt(abs(FEKO_farfield.Etheta*(1/376.73)).^2 + abs(FEKO_farfield.Ephi*(1/376.73)).^2);
% 
% 
% % Calculate now the E-field value here internal(Makarov Method)
% index = 0;
% for x_steps = x_grid
%     for y_steps = y_grid 
%         index = index + 1;
% 
%        
%         r = sqrt((x_steps)^2 + (y_steps)^2 + (z)^2);
%         % Uncomment below for additional debug output
%         %fprintf(sprintf('Calculating now the E-field at: (r,theta,phi) = (%2.f,%2.f,%2.f)\n',r,theta_degrees,phi_degrees));
%         
%         % Calculate the Electric field in spherical co-ordinates
%        [EfieldAtPointCartesian, HfieldAtPointCartesian] =  calculateEfieldAtPointRWG(Const, z, x_steps, y_steps, Solver_setup, Solution.mom.Isol);
%         Efield_vectors(index, :) = EfieldAtPointCartesian;
%         Hfield_vectors(index, :) = HfieldAtPointCartesian;  
%         
%         
%         % Calculate now the magnitude of the E-field vector. 
%         % Note: Change the unit of the E-field now fom V/m to V by 
%         % multiplying with the distance [r], at which
%         % the E-field was calculated.
%         Efield_magnitude(index) = sqrt(abs(r.*EfieldAtPointCartesian(1))^2 + ...
%             abs(r.*EfieldAtPointCartesian(2))^2 + ...
%             abs(r.*EfieldAtPointCartesian(3))^2);
%         Hfield_magnitude(index) = sqrt(abs(r.*HfieldAtPointCartesian(1))^2 + ...
%             abs(r.*HfieldAtPointCartesian(2))^2 + ...
%             abs(r.*HfieldAtPointCartesian(3))^2);
%         
% 
%         
%     end%for
% end%for
% 
% %error in difference between FEKO-Makarov & FEKO-Maaskant
% Efield_MAASKANT_errorx = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_Maaskant(:,1));
% Efield_MAASKANT_errory = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_Maaskant(:,2));
% Efield_MAASKANT_errorz = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_Maaskant(:,3));
% Efield_errorNormPercentagex = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_vectors(:,1));
% Efield_errorNormPercentagey = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_vectors(:,2));
% Efield_errorNormPercentagez = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_vectors(:,3));
% Hfield_errorNormPercentage = calculateErrorNormPercentage(FEKO_Hfield_vectors, Hfield_vectors);
% Efield_ACA_errorx = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_ACA(:,1));
% Efield_ACA_errory = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_ACA(:,2));
% Efield_ACA_errorz = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_ACA(:,3));
% 
% % Plot now the total E-field gri0d
% figure;
% hold on;
% grid on;
% 
% % Plot the normalised values
% plot_normalised_field = false;
% if (plot_normalised_field)
%     max_Efield_magnitude = max(Efield_magnitude);
%     max_FEKO_Efield_magnitude = max(FEKO_Efield_magnitude);
%     max_Efield_Maaskant_mag = max(Efield_Maaskant_mag);
%     max_Efield_ACA_mag = max(Efield_ACA_mag);
%     %max_FEKO_farfield_magnitude = max(FEKO_farfield_magnitude);
% else
%     % No normalisation applied, i.e. set the factor to 1.
%     max_Efield_magnitude = 1.0;
%     max_FEKO_Efield_magnitude = 1.0;
%     max_Efield_Maaskant_mag = 1.0;
%     max_Efield_ACA_mag = 1.0;
%    % max_FEKO_farfield_magnitude = 1.0;
% end%if
% 
% plot(1:total_efield_samples,Efield_magnitude./max_Efield_magnitude,'LineWidth',3);
% plot(1:FEKO_total_efield_samples,FEKO_Efield_magnitude./max_FEKO_Efield_magnitude,'x','LineWidth',3);
% plot(1:num_y_samples,Efield_Maaskant_mag./max_Efield_Maaskant_mag,'o','LineWidth',3);
% plot(1:num_y_samples,Efield_ACA_mag./max_Efield_ACA_mag,'+','LineWidth',3);
% %plot(1:FEKO_total_farfield_samples,FEKO_farfield_magnitude./max_FEKO_farfield_magnitude,'o','LineWidth',3);
% legend('SUN-EM-Makarov','FEKO (*.efe file)','SUN-EM-Maaskant','SUN-EM-ACA');
% set(get(gca, 'XLabel'), 'String', ('Sample index'));
% set(get(gca, 'YLabel'), 'String', ('|E-field| [V]'));
% 
% % There is a constant offset here. Try and figure out why we have this. Ideally, this should be around 1.
% % The following code was just required to extract the 2*pi constant that our code differs from FEKO
% % for the near-field calculation. Uncomment again if observing a large difference in the field values.
% % figure;
% % hold on;
% % grid on;
% % Efield_scaling_factor = Efield_magnitude./FEKO_Efield_magnitude;
% % plot(1:total_efield_samples,Efield_scaling_factor);
