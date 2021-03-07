% Author: Willem de la Bat (20692889@sun.ac.za)
% Project: MoM/PO feature extraction
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
Const = sunem_initialise('LbandHorn',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver       = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'LbandHorn.mat'; % Z-matrix calculated by FEKO
Const.FEKOstrfilename          = 'LbandHorn.str'; % I-vector calculated by FEKO
Const.FEKOrhsfilename          = 'LbandHorn.rhs'; % V-vector calculated by FEKO
Const.FEKOoutfilename          = 'LbandHorn.out';
%Const.FEKOefefilename          = 'bow_tie_antenna.efe';
%Const.FEKOffefilename          = 'bow_tie_antenna.ffe';
%Const.FEKOhfefilename          = 'bow_tie_antenna(1).hfe';%magnetic nearfield by FEKO




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

