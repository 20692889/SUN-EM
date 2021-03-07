% Author: Willem de la Bat (20692889@sun.ac.za)
% Project: Feature collection for hybrid solver classification
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
Const = sunem_initialise('PF_COEW',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver              = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
%Const.FEKOmatfilename          = 'PF_COEW.mat';
%Const.FEKOstrfilename          = 'PF_COEW.str';
%Const.FEKOrhsfilename          = 'PF_COEW.rhs';
Const.FEKOoutfilename          = 'PF_COEW.out';
Const.FEKOefefilename          = 'PF_COEW_NearField1.efe';
Const.FEKOffefilename          = 'PF_COEW_FF.ffe';
Const.FEKOosfilename           = 'PF_COEW_Currents1.os';


[Const, Solver_setup] = parseGEOout(Const);
[Idx] = splitNormVecsByLabel(Solver_setup);

norms_feed = [Solver_setup.triangle_normal_vector(1:431, :); Solver_setup.triangle_normal_vector(2556:2986, :); Solver_setup.triangle_normal_vector(5111:5541, :); Solver_setup.triangle_normal_vector(7666:8096, :)];
norms_MR = [Solver_setup.triangle_normal_vector(432:2555, :); Solver_setup.triangle_normal_vector(2987:5110, :); Solver_setup.triangle_normal_vector(5542:7665, :); Solver_setup.triangle_normal_vector(8097:10220, :)];
%[linVectMR, indMR] = GeoSurface(norms_MR);
%[linVectfeed, indfeed] = GeoSurface(norms_feed);
area_feed = [Solver_setup.triangle_area_m2(1:431, :); Solver_setup.triangle_area_m2(2556:2986, :); Solver_setup.triangle_area_m2(5111:5541, :); Solver_setup.triangle_area_m2(7666:8096, :)];
area_feed_tot = sum(area_feed, 'all');

area_MR = [Solver_setup.triangle_area_m2(432:2555, :); Solver_setup.triangle_area_m2(2987:5110, :); Solver_setup.triangle_area_m2(5542:7665, :); Solver_setup.triangle_area_m2(8097:10220, :)];
area_MR_tot = sum(area_MR, 'all');
num_tri_feed = size(norms_feed,1);
num_tri_MR = size(norms_MR,1);

std_feed_norms = std(norms_feed);
std_MR_norms = std(norms_MR);
CV_feed = std_feed_norms/mean(norms_feed);
CV_MR = std_MR_norms/mean(norms_MR);
%NormVecs = Solver_setup.triangle_normal_vector;
%[linVect, ind] = GeoSurface(Solver_setup.triangle_normal_vector);
%ind(end+1) = length(Solver_setup.triangle_normal_vector);
%[area_vec] = CalcFaceAreas(Solver_setup,ind);
% The Following file is used to port solutions to FEKO 
% (for post-processing in POSTFEKO).
% TO-DO: [DL] Add this.
% Const.output_strfilename    = '';
% Const.writeFEKOstrfile = [0 0 0 0];
r = 100;%100;

% Read the FEKO data from a *.efe file for comparison
FEKO_EfieldAtPoint = parseFEKOefefile(Const, Const.FEKOefefilename);
FEKO_total_efield_samples = FEKO_EfieldAtPoint.number_r_samples * ...
    FEKO_EfieldAtPoint.number_theta_samples * FEKO_EfieldAtPoint.number_phi_samples;
FEKO_Efield_magnitude = sqrt(abs(r.*FEKO_EfieldAtPoint.Er).^2 + ...
    abs(r.*FEKO_EfieldAtPoint.Etheta).^2 + ...
    abs(r.*FEKO_EfieldAtPoint.Ephi).^2);

% Read also now FEKO's far field values here.
FEKO_farfield = parseFEKOffefile(Const, Const.FEKOffefilename);
FEKO_total_farfield_samples =  FEKO_farfield.number_theta_samples * FEKO_farfield.number_phi_samples;
FEKO_farfield_magnitude = sqrt(abs(FEKO_farfield.Etheta).^2 + abs(FEKO_farfield.Ephi).^2);

%Read also now FEKO's current values here.
FEKO_currents = parseFEKOosfile (Const, Const.FEKOosfilename);

% save('FEKO_farfield_data','FEKO_farfield')
% save('FEKO_nearfield_data', 'FEKO_EfieldAtPoint')
% save('FEKO_current_data','FEKO_currents')
FEKO_farfield_data = load('FEKO_farfield_data.mat');
FEKO_nearfield_data = load('FEKO_nearfield_data.mat');
FEKO_current_data = load('FEKO_current_data.mat');
% 

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%Calculate Error Norms between MLFMM and MOM/PO 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

Efields_MLFMM_ff = [FEKO_farfield_data.FEKO_farfield.Etheta FEKO_farfield_data.FEKO_farfield.Ephi];
Efields_MOMPO_ff = [FEKO_farfield.Etheta FEKO_farfield.Ephi];
Efield_MLFMM_nf = [FEKO_nearfield_data.FEKO_EfieldAtPoint.Er FEKO_nearfield_data.FEKO_EfieldAtPoint.Etheta FEKO_nearfield_data.FEKO_EfieldAtPoint.Ephi];
Efield_MOMPO_nf = [FEKO_EfieldAtPoint.Er FEKO_EfieldAtPoint.Etheta FEKO_EfieldAtPoint.Ephi];
Current_MLFMM = [FEKO_current_data.FEKO_currents.Jx FEKO_current_data.FEKO_currents.Jy FEKO_current_data.FEKO_currents.Jz];
Current_MOMPO = [FEKO_currents.Jx FEKO_currents.Jy FEKO_currents.Jz];
% 
Err_ff_phi = calculateErrorNormPercentage(FEKO_farfield_data.FEKO_farfield.Ephi, FEKO_farfield.Ephi);
Err_ff_theta = calculateErrorNormPercentage(FEKO_farfield_data.FEKO_farfield.Etheta, FEKO_farfield.Etheta);
Err_ff_tot = calculateMatrixErrorNormPercentage(Efields_MLFMM_ff, Efields_MOMPO_ff);
% 
Err_nf_x = calculateErrorNormPercentage(FEKO_nearfield_data.FEKO_EfieldAtPoint.Er, FEKO_EfieldAtPoint.Er);
Err_nf_y = calculateErrorNormPercentage(FEKO_nearfield_data.FEKO_EfieldAtPoint.Ephi, FEKO_EfieldAtPoint.Ephi);
Err_nf_z = calculateErrorNormPercentage(FEKO_nearfield_data.FEKO_EfieldAtPoint.Etheta, FEKO_EfieldAtPoint.Etheta);
Err_nf_tot = calculateMatrixErrorNormPercentage(Efield_MLFMM_nf, Efield_MOMPO_nf);

Err_cur_x =  calculateErrorNormPercentage(FEKO_current_data.FEKO_currents.Jx, FEKO_currents.Jx);
Err_cur_y =  calculateErrorNormPercentage(FEKO_current_data.FEKO_currents.Jy, FEKO_currents.Jy);
Err_cur_z =  calculateErrorNormPercentage(FEKO_current_data.FEKO_currents.Jz, FEKO_currents.Jz);
Err_cur_tot = calculateMatrixErrorNormPercentage(Current_MLFMM, Current_MOMPO);


