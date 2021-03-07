% Author: Danie Ludick (dludick@sun.ac.za)
% Project: PEC plate array example
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
Const = sunem_initialise('bow_tie_antenna',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver       = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'bow_tie_antenna.mat'; % Z-matrix calculated by FEKO
Const.FEKOstrfilename          = 'bow_tie_antenna.str'; % I-vector calculated by FEKO
Const.FEKOrhsfilename          = 'bow_tie_antenna.rhs'; % V-vector calculated by FEKO
Const.FEKOoutfilename          = 'bow_tie_antenna.out';
Const.FEKOefefilename          = 'bow_tie_antenna.efe';
Const.FEKOffefilename          = 'bow_tie_antenna.ffe';
Const.FEKOhfefilename          = 'bow_tie_antenna(1).hfe';%magnetic nearfield by FEKO
% Const.FEKO_mag_far = load('FEKO_mag_Far.mat');
% Const.FEKO_mag_mid = load('FEKO_mag_mid.mat');
% Const.FEKO_mag_near = load('FEKO_mag_near.mat');
% Const.FM_mag_Far = load('FM_mag_Far.mat');
% Const.FM_mag_mid = load('FM_mag_mid.mat');
% Const.FM_mag_near = load('FM_mag_near.mat');
% Const.DM_mag_Far = load('DM_mag_Far.mat');
% Const.DM_mag_mid = load('DM_mag_mid.mat');
% Const.DM_mag_near = load('DM_mag_near.mat');



% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMifbmomstrfilename   = '';%'ifbmom_bow_tie_array_6by1.str';
Const.SUNEMcbfmstrfilename     = '';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% TO-DO: Setup some documentation for this - also assign default values if they are not
% defined explicitely here.
%Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.


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

%Const.useACA = 1;
%Observation grid cuts
z=1000;
%x_grid = -0.125:0.025:0.1;
%y_grid = 0.1:0.025:0.325;
y_grid=0.25:1:4.25;
%y_grid=-149.775:10:149.225;
%y_grid=-1499.775:100:1499.225;
%y_grid=-7499.775:500:7499.225;
%y_grid=-59.775:1:59.225;
%y_grid=-29.775:1:29.225;
%x_grid=-50:10:40;
x_grid = -25:10:15;
%x_grid=-100:10:90;
num_y_samples = length(y_grid);
num_x_samples = length(x_grid);
num_samples = num_y_samples*num_x_samples;
x_rep = repmat(x_grid, [num_y_samples,1]).';
x_rep = reshape(x_rep, [],1);
y_rep = reshape(repmat(reshape(y_grid',[],1),1,5)',[],size(y_grid,1));
z_rep = repelem(z, num_samples).';
%vector of observation points    
Solver_setup.r_c_obs = [x_rep, y_rep, z_rep];
%function to add RWG elements at observation point
[Solver_setup] = AddRWGatObservationPt(Const, Solver_setup);
%extract SUNEM mom matrix (for now it calculates the Zmn matrix between
%antenna surface RWGs and observation point RWGs whos vertices all lie on the
%y-axis
tic
[Const, zMatricesSUNEMx, zMatricesSUNEMy, zMatricesSUNEMz, yVectorsSUNEM] = extractSUNEMMoMmatrixEq(Const, Solver_setup);
TotZmnTime = toc;

SourceVecInd = 1:Solver_setup.num_metallic_edges;
ObservVecInd = 1:length(Solver_setup.r_c_obs);
tic
%[Zmn,Ux,Vx] = calcZmn(Const,1, Solver_setup, zMatricesSUNEMx, Solver_setup.frequencies.samples, zMatricesSUNEMx.mBasis, zMatricesSUNEMx.nBasis, ObservVecInd, SourceVecInd);
[Zmn,Uy,Vy] = calcZmn(Const,2, Solver_setup, zMatricesSUNEMy, Solver_setup.frequencies.samples, zMatricesSUNEMy.mBasis, zMatricesSUNEMy.nBasis, ObservVecInd, SourceVecInd);
%[Zmn,Uz,Vz] = calcZmn(Const,3, Solver_setup, zMatricesSUNEMz, Solver_setup.frequencies.samples, zMatricesSUNEMz.mBasis, zMatricesSUNEMz.nBasis, ObservVecInd, SourceVecInd);
TotUVtime = toc; 
Const.useACA = false;


% --------------------------------------------------------------------------------------------------
% Run the EM solver
% --------------------------------------------------------------------------------------------------
 [Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
% 
% Postprocess the results, e.g. calculate the Electric field
% --------------------------------------------------------------------------------------------------

%[Uvector, Vvector] = AdapadtiveCrossApproximation(zMatricesSUNEMy.values, zMatricesSUNEMy.mBasis, zMatricesSUNEMy.nBasis);
%
%[Uapprox, Vapprox] = AdapadtiveCrossApproximation(zMatricesSUNEMx.values, zMatricesSUNEMx.mBasis, zMatricesSUNEMx.nBasis);
%Eqn (4) in "Field Computations Through the ACA Algorithm" by Maaskant, R. ; Lancelotti, V. (2015) 

lp = Solver_setup.lengthP;
TotTimeACA = 0;
TotTimeMaas = 0;
CalcExACAtime = 0;
CalcEyACAtime = 0;
CalcEzACAtime = 0;
%iter = 0;
for i=1:100
    tic
 %   mul = (Vx*Solution.mom.Isol);
 %   ExACA = -(1/(lp*lp))*Ux*mul;
    %CalcExACAtime = CalcExACAtime + toc;
    %tic
    mul = (Vy*Solution.mom.Isol);
    EyACA = -(1/(lp*lp))*Uy*mul;
    %CalcEyACAtime = CalcEyACAtime + toc;
    %tic
 %   mul = (Vz*Solution.mom.Isol);
 %   EzACA = -(1/(lp*lp))*Uz*mul;
    %CalcEzACAtime = CalcEzACAtime + toc;
   % aTotTimeACA = toc;
    TotTimeACA = TotTimeACA + toc;
    %TotTimeACA = TotTimeACA + CalcExACAtime + CalcEyACAtime + CalcEzACAtime;
    tic
   % Ex = -(1/(lp*lp))*zMatricesSUNEMx.values*Solution.mom.Isol;
    %CalcExMaastime = toc;
   % tic
    Ey = -(1/(lp*lp))*zMatricesSUNEMy.values*Solution.mom.Isol;
   % CalcEyMaastime = toc;
   % tic
   % Ez = -(1/(lp*lp))*zMatricesSUNEMz.values*Solution.mom.Isol;
   % CalcEzMaastime = toc;
   % aTotTimeMaas = toc;
    TotTimeMaas = TotTimeMaas + toc;
    %iter = iter + 1;
   % TotTimeMaas = TotTimeMaas + CalcExMaastime + CalcEyMaastime + CalcEzMaastime;
end

%avgXtime = CalcExACAtime/100;
%avgYtime = CalcEyACAtime/100;
%avgZtime = CalcEzACAtime/100;
avgACAtime = TotTimeACA/100;
avgMAAStime = TotTimeMaas/100;
%(1/(lp*lp))*

a = zeros(1,300);
ExACA = complex(a,a).';
EzACA = complex(a,a).';
Ex = complex(a,a).';
Ez = complex(a,a).';
Efield_ACA = [ExACA, EyACA, EzACA];
Efield_Maaskant = [Ex, Ey, Ez];


r = sqrt((x_rep).^2 + (y_rep).^2 + (z_rep).^2);


%r = 100;%100;
Efield_ACA_mag = sqrt(abs(r.*Efield_ACA(:,1)).^2 + ...
            abs(r.*Efield_ACA(:,2)).^2 + ...
            abs(r.*Efield_ACA(:,3)).^2);
Efield_Maaskant_mag= sqrt(abs(r.*Efield_Maaskant(:,1)).^2 + ...
            abs(r.*Efield_Maaskant(:,2)).^2 + ...
            abs(r.*Efield_Maaskant(:,3)).^2);
% Loop over a few theta and phi points and compare the results with that of FEKO
%theta_grid = 90:1:90;
%phi_grid = 1:1:91;

num_x_samples = length(x_grid);
num_y_samples = length(y_grid);
num_z_samples = length(z);
%num_theta_samples = length(theta_grid);
%num_phi_samples = length(phi_grid);
total_efield_samples = num_x_samples*num_y_samples*num_z_samples;
Efield_magnitude = zeros(total_efield_samples,1);
%Hfield_magnitude = zeros(total_efield_samples,1);

%difference = zeros(total_efield_samples,3);

%Hfield_vectors = zeros(total_efield_samples,3);
Efield_vectors = zeros(total_efield_samples,3);
FEKO_Efield_vectors = zeros(total_efield_samples,3);
%FEKO_Hfield_vectors = zeros(total_efield_samples,3);

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

% Read the FEKO data from a *.efe file for comparison
FEKO_EfieldAtPoint = parseFEKOefefile(Const, Const.FEKOefefilename);
FEKO_total_efield_samples = FEKO_EfieldAtPoint.number_x_samples * ...
    FEKO_EfieldAtPoint.number_y_samples * FEKO_EfieldAtPoint.number_z_samples;
FEKO_Efield_magnitude = sqrt(abs(r.*FEKO_EfieldAtPoint.Ex).^2 + ...
    abs(r.*FEKO_EfieldAtPoint.Ey).^2 + ...
    abs(r.*FEKO_EfieldAtPoint.Ez).^2);

%vectors with Efield values from FEKO
FEKO_Efield_vectors(:, 1) = FEKO_EfieldAtPoint.Ex;
FEKO_Efield_vectors(:, 2) = FEKO_EfieldAtPoint.Ey;
FEKO_Efield_vectors(:, 3) = FEKO_EfieldAtPoint.Ez;

% --------------------------------------------------------------------------------------------------
% Postprocess the results, e.g. calculate the Electric field
% --------------------------------------------------------------------------------------------------

%if (false)

   % r = 100;%100;
    % Loop over a few theta and phi points and compare the results with that of FEKO
    theta_grid = 1;
    phi_grid = 1;

    num_theta_samples = 1;
    num_phi_samples = 1;
   % total_efield_samples = num_theta_samples*num_phi_samples;
   % Efield_magnitude = zeros(total_efield_samples,1);

    % Read the FEKO data from a *.efe file for comparison
    %FEKO_EfieldAtPoint = parseFEKOefefile(Const, Const.FEKOefefilename);
   % FEKO_total_efield_samples = FEKO_EfieldAtPoint.number_r_samples * ...
   %     FEKO_EfieldAtPoint.number_theta_samples * FEKO_EfieldAtPoint.number_phi_samples;
   % FEKO_Efield_magnitude = sqrt(abs(r.*FEKO_EfieldAtPoint.Er).^2 + ...
   %     abs(r.*FEKO_EfieldAtPoint.Etheta).^2 + ...
   %     abs(r.*FEKO_EfieldAtPoint.Ephi).^2);

    % Read also now FEKO's far field values here.
%     FEKO_farfield = parseFEKOffefile(Const, Const.FEKOffefilename);
%     FEKO_total_farfield_samples =  FEKO_farfield.number_theta_samples * FEKO_farfield.number_phi_samples;
%     FEKO_farfield_magnitude = sqrt(abs(FEKO_farfield.Etheta).^2 + abs(FEKO_farfield.Ephi).^2);

    % Calculate now the E-field value here internal
%     calcdmtime = 0;
%    % for i=1:100
%         tic
%     index = 0;
%     for y_steps = y_grid
%         for x_steps = x_grid
%             index = index + 1;
% 
%             r = sqrt((x_steps)^2 + (y_steps)^2 + (z)^2);
%             % Uncomment below for additional debug output
%             %fprintf(sprintf('Calculating now the E-field at: (r,theta,phi) = (%2.f,%2.f,%2.f)\n',r,theta_degrees,phi_degrees));
%             
%             % Calculate the Electric field in spherical co-ordinates
%             [EfieldAtPointCartesian] =  calculateEfieldAtPointRWG(Const, z, x_steps, y_steps, Solver_setup, Solution.mom.Isol);
%             Efield_vectors(index, :) = EfieldAtPointCartesian;
%           %  Hfield_vectors(index, :) = HfieldAtPointCartesian;
%             
%             
%             % Calculate now the magnitude of the E-field vector.
%             % Note: Change the unit of the E-field now fom V/m to V by
%             % multiplying with the distance [r], at which
%             % the E-field was calculated.
%             Efield_magnitude(index) = sqrt(abs(r.*EfieldAtPointCartesian(1))^2 + ...
%                 abs(r.*EfieldAtPointCartesian(2))^2 + ...
%                 abs(r.*EfieldAtPointCartesian(3))^2);
%             %Hfield_magnitude(index) = sqrt(abs(r.*HfieldAtPointCartesian(1))^2 + ...
%             %    abs(r.*HfieldAtPointCartesian(2))^2 + ...
%             %    abs(r.*HfieldAtPointCartesian(3))^2);
%         end%for
%     end%for
%     calcdmtime = calcdmtime + toc;
%     %end
%     CalcDMtime = calcdmtime/100;

    %error in difference between FEKO-Makarov & FEKO-Maaskant
%      Efield_FM_errorx = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_Maaskant(:,1));
%      Efield_FM_errory = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_Maaskant(:,2));
%      Efield_FM_errorz = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_Maaskant(:,3));
%      Efield_DM_errorx = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_vectors(:,1));
%      Efield_DM_errory = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_vectors(:,2));
%      Efield_DM_errorz = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_vectors(:,3));
%      %Hfield_errorNormPercentage = calculateErrorNormPercentage(FEKO_Hfield_vectors, Hfield_vectors);
%      Efield_ACA_errorx = calculateErrorNormPercentage(FEKO_Efield_vectors(:,1), Efield_ACA(:,1));
%      Efield_ACA_errory = calculateErrorNormPercentage(FEKO_Efield_vectors(:,2), Efield_ACA(:,2));
%      Efield_ACA_errorz = calculateErrorNormPercentage(FEKO_Efield_vectors(:,3), Efield_ACA(:,3));
    % Plot now the total E-field grid
    
   % E_veld_DM_error_total = calculateMatrixErrorNormPercentage(FEKO_Efield_vectors, Efield_vectors);
    E_veld_FM_error_total = calculateMatrixErrorNormPercentage(FEKO_Efield_vectors, Efield_Maaskant);
    E_veld_ACA_error_total = calculateMatrixErrorNormPercentage(FEKO_Efield_vectors, Efield_ACA);
%     %figure;
%     %hold on;
%     %grid on;
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
%     % max_FEKO_farfield_magnitude = 1.0;
% end%if
% % save('FEKO_mag_mid.mat','FEKO_Efield_magnitude');
% % save('FM_mag_mid.mat','Efield_Maaskant_mag');
% % save('DM_mag_mid.mat','Efield_magnitude');
% %   FEKO_mags_zsteps = [Const.FEKO_mag_near.FEKO_Efield_magnitude, Const.FEKO_mag_mid.FEKO_Efield_magnitude, Const.FEKO_mag_far.FEKO_Efield_magnitude];
% %   DM_mags_zsteps = [Const.DM_mag_near.Efield_magnitude, Const.DM_mag_mid.Efield_magnitude, Const.DM_mag_Far.Efield_magnitude];
% % % 
% %   FM_mags_zsteps = [Const.FM_mag_near.Efield_Maaskant_mag,  Const.FM_mag_mid.Efield_Maaskant_mag,  Const.FM_mag_Far.Efield_Maaskant_mag];
% % 
% %   figure(1);
% % % % %surf();
% %   xrep3 = [x_rep, x_rep, x_rep];
% %   yrep3 = [y_rep, y_rep, y_rep];
% %   error_vector_ACAdelY1 = [2.416, 2.422, 11.3265, 61.5907];
% %  x = logspace(-3,1,4);
% % % y(1) = error_vector_ACAdelY1(x(1))
% %  semilogx(x,y,'o','MarkerFaceColor',[0 0.447 0.741])
% %  grid on
%    i = 0;
%    MagMat = zeros(num_y_samples, num_x_samples);
%    MagMatDM = zeros(num_y_samples, num_x_samples);
%    MagMatFM = zeros(num_y_samples, num_x_samples);
% %   MagMatACA = zeros(num_y_samples, num_x_samples);
%    for y = 1:30
%        for x=1:10
%            i=i+1;
%            MagMat(y,x)= FEKO_Efield_magnitude(i);
%            MagMatDM(y,x)= Efield_magnitude(i);
%            MagMatFM(y,x)= Efield_Maaskant_mag(i);
%           % MagMatACA(y,x)= Efield_ACA_mag(i);
%            
%        end
%    end
%   [X,Y] = meshgrid(x_grid, y_grid);
% %  
% % %  %surfaceDM
%   figure(1);
%   surf(X, Y, MagMat,'FaceAlpha',0.2,'FaceColor','r');
%   hold on
% % % % %surf(X, Y, MagMatACA,'FaceAlpha',0.2,'FaceColor','b');
%   surf(X, Y, MagMatDM,'FaceAlpha',0.2,'FaceColor','y');
% % %  title('Surface DM v FEKO rad near');
% % % % 
%   hold off
% % %  
% % %  %surfaceFM
%   figure(2);
%   surf(X, Y, MagMat,'FaceAlpha',0.2,'FaceColor','r');
%   hold on
% % % % 
%  surf(X, Y, MagMatFM,'FaceAlpha',0.2,'FaceColor','b');
% % %  title('Surface FM v FEKO rad near');
% % % % 
%  hold off
%  %contour3(xrep3, yrep3, FEKO_mags_zsteps);
% figure(1);
% % surf(X, Y, MagMat,'FaceAlpha',0.6,'FaceColor','r');
% % hold on
% % 
%  surf(X, Y, MagMatACA,'FaceAlpha',0.2,'FaceColor','b');
%  title('Surface ACA tol0 v FEKO rad near');
% 
% hold off
%[X,Y] = meshgrid(x_grid, y_grid);
%Z = griddata(x,y,z,X,Y);
%contour3(X,Y,Z);
% 
%  figure(2);
%  plot3(x_rep,y_rep, FEKO_mags_zsteps(:,1),'-k','lineWidth',0.5);
%    hold on
%    grid on
%  plot3(x_rep,y_rep, FEKO_mags_zsteps(:,2),'-k','lineWidth',0.5);
%  plot3(x_rep,y_rep, FEKO_mags_zsteps(:,3),'-k','lineWidth',0.5);
%  hold off
% %  plot3(x_rep, y_rep,Efield_magnitude,':y','lineWidth',2);
%   xlabel('X-position [m]')
%   ylabel('Y-position [m]')
%   zlabel('|E-Field| [V]')
%   legend('FEKO','Dipole Moment')
%   hold off
% figure(2);
%   plot3(x_rep,y_rep, FEKO_Efield_magnitude,':r','lineWidth',2);
%  hold on
%  grid on
%  plot3(x_rep,y_rep,Efield_magnitude,':y','lineWidth',2);
%  plot3(x_rep, y_rep,Efield_Maaskant_mag,':b','lineWidth',2);
%  xlabel('X-position [m]')
%  ylabel('Y-position [m]')
%  zlabel('|E-Field| [V]')
%  legend('FEKO','Dipole Moment','Direct Full-MoM')
%  hold off
 
% b = FEKO_Efield_vectors(:,2);
% FEKO_y_mags = [b(6),b(16),b(26),b(36),b(46),b(56),b(66),b(76),b(86),b(96),b(106),b(116),b(126),b(136),b(146),b(156),b(166),b(176),b(186),b(196),b(206),b(216),b(226),b(236),b(246),b(256),b(266),b(276),b(286),b(296)];
% b=Efield_Maaskant(:,2);
% FM_y_mags = [b(6),b(16),b(26),b(36),b(46),b(56),b(66),b(76),b(86),b(96),b(106),b(116),b(126),b(136),b(146),b(156),b(166),b(176),b(186),b(196),b(206),b(216),b(226),b(236),b(246),b(256),b(266),b(276),b(286),b(296)];
% b=Efield_vectors(:,2);
% DM_y_mags = [b(6),b(16),b(26),b(36),b(46),b(56),b(66),b(76),b(86),b(96),b(106),b(116),b(126),b(136),b(146),b(156),b(166),b(176),b(186),b(196),b(206),b(216),b(226),b(236),b(246),b(256),b(266),b(276),b(286),b(296)];
 
%complex plots in x
% b = FEKO_Efield_vectors(:,1);
%  FEKO_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_Maaskant(:,1);
%  FM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_vectors(:,1);
%  DM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_ACA(:,1);
%  ACA_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
% figure(2);
% % plot3(y_grid, real(FEKO_y_comp),imag(FEKO_y_comp),'rx','lineWidth',2);
% % hold on
% % grid on
% %  plot3(y_grid, real(FM_y_comp),imag(FM_y_comp),'b*','lineWidth',2);
% %  plot3(y_grid, real(DM_y_comp),imag(DM_y_comp),'y+','lineWidth',2);
%  plot3(y_grid, real(ACA_y_comp),imag(ACA_y_comp),'bo','lineWidth',2);
%  grid on   
%  xlabel('Y-position [m]')
%  ylabel('Im\{Ex [V/m]\}')
%  zlabel('Re\{Ex [V/m]\}')
% % legend('FEKO','ACA')
%  title('Complex Ex ACA rad near')
% % hold off
% %complex plots in y 
%  b = FEKO_Efield_vectors(:,2);
%  FEKO_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_Maaskant(:,2);
%  FM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_vectors(:,2);
%  DM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_ACA(:,2);
%  ACA_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  figure(3);
% % plot3(y_grid, real(FEKO_y_comp),imag(FEKO_y_comp),'rx','lineWidth',2);
% % hold on
% % grid on
%  plot3(y_grid, real(ACA_y_comp),imag(ACA_y_comp),'bo','lineWidth',2);
%  grid on
% % plot3(y_grid, real(FM_y_comp),imag(FM_y_comp),'b*','lineWidth',2);
% % plot3(y_grid, real(DM_y_comp),imag(DM_y_comp),'y+','lineWidth',2);
% 
%  xlabel('Y-position [m]')
%  ylabel('Im\{Ey [V/m]\}')
%  zlabel('Re\{Ey [V/m]\}')
%  %legend('FEKO','ACA')
%  title('Complex Ey ACA rad near')
% % hold off
%  %complex plots in z
% b = FEKO_Efield_vectors(:,3);
%  FEKO_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_Maaskant(:,3);
%  FM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  b=Efield_vectors(:,3);
%  DM_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)]; 
%  b=Efield_ACA(:,3);
%  ACA_y_comp = [b(5),b(15),b(25),b(35),b(45),b(55),b(65),b(75),b(85),b(95)];
%  figure(4);
%  %plot3(y_grid, real(FEKO_y_comp),imag(FEKO_y_comp),'rx','lineWidth',2);
%  %hold on
% % grid on
% % plot3(y_grid, real(FM_y_comp),imag(FM_y_comp),'b*','lineWidth',2);
% % plot3(y_grid, real(DM_y_comp),imag(DM_y_comp),'y+','lineWidth',2);
%  plot3(y_grid, real(ACA_y_comp),imag(ACA_y_comp),'bo','lineWidth',2);
%  grid on
%  xlabel('Y-position [m]')
%  ylabel('Im\{Ez [V/m]\}')
%  zlabel('Re\{Ez [V/m]\}')
%  title('Complex Ez ACA rad near')
 %legend('FEKO','ACA')
 %hold off
 
 %ra = abs(FEKO_Efield_vectors(:,1));
%thet = angle(FEKO_Efield_vectors(:,1));
% figure(1);
% polarplot(FEKO_Efield_vectors(30:50,1),'k+');
% hold on
% polarplot(Efield_vectors(30:50,1),'bo','MarkerSize',10);
% legend('FEKO','Dipole Moment');
% hold off
% figure(2);
% polarplot(FEKO_Efield_vectors(30:50,2),'k+');
% hold on
% polarplot(Efield_vectors(30:50,2),'bo','MarkerSize',10);
% legend('FEKO','Dipole Moment');
% hold off
% figure(3);
% polarplot(FEKO_Efield_vectors(30:50,3),'k+');
% hold on
% polarplot(Efield_vectors(30:50,3),'bo','MarkerSize',10);
% legend('FEKO','Dipole Moment');
% hold off
% figure(4);
% polarplot(FEKO_Efield_vectors(30:50,1),'k+');
% hold on
% polarplot(Efield_Maaskant(30:50,1),'rx','MarkerSize',20);
% polarplot(Efield_vectors(30:50,1),'bo');
% legend('FEKO','Direct Full-MoM','Dipole Moment');
% hold off
% figure(5);
% polarplot(FEKO_Efield_vectors(30:50,2),'k+');
% hold on
% polarplot(Efield_Maaskant(30:50,2),'rx','MarkerSize',20);
% polarplot(Efield_vectors(30:50,2),'bo');
% legend('FEKO','Direct Full-MoM','Dipole Moment');
% hold off
% figure(6);
% polarplot(FEKO_Efield_vectors(30:50,3),'k+');
% hold on
% polarplot(Efield_Maaskant(30:50,3),'rx','MarkerSize',20);
% polarplot(Efield_vectors(30:50,3),'bo');
% legend('FEKO','Direct Full-MoM','Dipole Moment');
% hold off
% figure(1);
% plot(y_grid,FEKO_y_mags,'-x','LineWidth',2);
% hold on
% plot(y_grid, ACA_y_mags,'-o','LineWidth',4);
% set(get(gca, 'XLabel'), 'String', ('Y-position'));
% set(get(gca, 'YLabel'), 'String', ('|E-field| [V]'));
% hold off
%plot(1:total_efield_samples,Efield_magnitude./max_Efield_magnitude,'LineWidth',3);
%plot(1:FEKO_total_efield_samples,FEKO_Efield_magnitude./max_FEKO_Efield_magnitude,'x','LineWidth',4);
%plot(1:num_samples,Efield_Maaskant_mag./max_Efield_Maaskant_mag,'o','LineWidth',2);
%plot(1:num_samples,Efield_ACA_mag./max_Efield_ACA_mag,'+','LineWidth',1);
%plot(1:FEKO_total_farfield_samples,FEKO_farfield_magnitude./max_FEKO_farfield_magnitude,'o','LineWidth',3);
%legend('Dipole Moment','FEKO (*.efe file)','Direct Full MoM','Direct ACA');
%set(get(gca, 'XLabel'), 'String', ('Sample index'));
%set(get(gca, 'YLabel'), 'String', ('|E-field| [V]'));
%figure(1);
%plot(real(Efield_vectors(1:50,1)),imag(Efield_vectors(1:50,1)),'b+','LineWidth',2,'MarkerSize',10);
%hold on
%plot(real(FEKO_Efield_vectors(1:50,1)),imag(FEKO_Efield_vectors(1:50,1)),'ko','lineWidth',2);
%plot(real(Efield_ACA(1:50,1)),imag(Efield_ACA(1:50,1)),'rx','lineWidth',3);
%plot(real(Efield_Maaskant(1:50,1)),imag(Efield_Maaskant(1:50,1)),'g*','lineWidth',1);
%legend('Dipole Moment Method','FEKO Baseline');
%legend('Dipole Moment','FEKO','Direct ACA','Direct Full MoM');
%set(get(gca, 'XLabel'), 'String', ('Re\{E-field strength [V/m]\}'));
%set(get(gca, 'YLabel'), 'String', ('Im\{E-field strength [V/m]\}'));
%hold off
%figure(2);
%plot(real(Efield_vectors(1:50,2)),imag(Efield_vectors(1:50,2)),'b+','LineWidth',2,'MarkerSize',10);
%hold on
%plot(real(FEKO_Efield_vectors(1:50,2)),imag(FEKO_Efield_vectors(1:50,2)),'ko','lineWidth',2);
%plot(real(Efield_ACA(1:50,2)),imag(Efield_ACA(1:50,2)),'rx','lineWidth',3);
%plot(real(Efield_Maaskant(1:50,2)),imag(Efield_Maaskant(1:50,2)),'g*','lineWidth',1);
%legend('Dipole Moment','FEKO','Direct ACA','Direct Full MoM');
%legend('Dipole Moment Method','FEKO Baseline');
%set(get(gca, 'XLabel'), 'String', ('Re\{E-field strength [V/m]\}'));
%set(get(gca, 'YLabel'), 'String', ('Im\{E-field strength [V/m]\}'));
%hold off
%figure(3);
%plot(real(Efield_vectors(1:50,3)),imag(Efield_vectors(1:50,3)),'b+','LineWidth',2,'MarkerSize',10);
%hold on
%plot(real(FEKO_Efield_vectors(1:50,3)),imag(FEKO_Efield_vectors(1:50,3)),'ko','lineWidth',2);
%plot(real(Efield_ACA(1:50,3)),imag(Efield_ACA(1:50,3)),'rx','lineWidth',3);
%plot(real(Efield_Maaskant(1:50,3)),imag(Efield_Maaskant(1:50,3)),'g*','lineWidth',1);
%legend('Dipole Moment','FEKO','Direct ACA','Direct Full MoM');
%legend('Dipole Moment Method','FEKO Baseline');
%set(get(gca, 'XLabel'), 'String', ('Re\{E-field strength [V/m]\}'));
%set(get(gca, 'YLabel'), 'String', ('Im\{E-field strength [V/m]\}'));
%hold off
    % There is a constant offset here. Try and figure out why we have this. Ideally, this should be around 1.
    % The following code was just required to extract the 2*pi constant that our code differs from FEKO
    % for the near-field calculation. Uncomment again if observing a large difference in the field values.
    % figure;
    % hold on;
    % grid on;
    % Efield_scaling_factor = Efield_magnitude./FEKO_Efield_magnitude;
    % plot(1:total_efield_samples,Efield_scaling_factor);

%end%if