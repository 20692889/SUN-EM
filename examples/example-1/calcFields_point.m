
function [EfieldAtPointSpherical, HfieldAtPointSpherical] = calcFields_point(Const, r, theta_degrees, phi_degrees, ...
    Solver_setup, Isol)

    EfieldAtPointCartesian = zeros(Solver_setup.frequencies.freq_num, 3);
    EfieldAtPointSpherical = zeros(Solver_setup.frequencies.freq_num, 3); 
    HfieldAtPointCartesian = zeros(Solver_setup.frequencies.freq_num, 3);
    HfieldAtPointSpherical = zeros(Solver_setup.frequencies.freq_num, 3);
    [Px,Py,Pz] = transformSpericalCoordinateToCartesian(r,theta_degrees*Const.DEG2RAD,phi_degrees*Const.DEG2RAD);

        for freq_index = 1:Solver_setup.frequencies.freq_num

            lambda = Const.C0/Solver_setup.frequencies.samples(freq_index);
            k = 2*pi/lambda;
            K = Const.ETA_0/(4*pi);
            A = 1j*k/(4*pi);
            U = 1j*k;
            
            for m=1:Solver_setup.num_metallic_edges
                
                  lengthM = Solver_setup.rwg_basis_functions_length_m(m);
                  Im = Isol(m);
                  triangle_plus = Solver_setup.rwg_basis_functions_trianglePlus(m);
                  triangle_minus = Solver_setup.rwg_basis_functions_triangleMinus(m);
                  Point1 = Solver_setup.triangle_centre_point(triangle_plus,:);
                  Point2 = Solver_setup.triangle_centre_point(triangle_minus,:);
                  DipoleCenter(:,m)=0.5*(Point1+Point2);
                  DipoleMoment(:,m)=lengthM*Im*(-Point1+Point2);
            end
            
            ObservationPoint=[Px; Py; Pz];
            [E,H]=point(ObservationPoint,K,A,U,DipoleMoment,DipoleCenter);

            %find the sum of all dipole contributions
            EfieldAtPointCartesian(freq_index,:)=sum(E,2).'; 
            HfieldAtPointCartesian(freq_index,:)=sum(H,2).';
           % EfieldAtPointCartesian(freq_index,:) = EField.';
            %HfieldAtPointCartesian(freq_index,:) = Hfield.';
            EfieldAtPointSpherical(freq_index,:) = transformCartesianVectorToSpherical(EfieldAtPointCartesian(freq_index,:), Px, Py, Pz);
            HfieldAtPointSpherical(freq_index,:) = transformCartesianVectorToSpherical(HfieldAtPointCartesian(freq_index,:), Px, Py, Pz);
%k=omega/c_;
%K=j*k;
        end


end

function[EField, HField]=...
    point(Point,ConstantE,ConstantH,K,DipoleMoment,DipoleCenter)

%POINT Radiated/scattered field at a point of a dipole array 
%   or a single dipole. Gives exact near- and far-fields. Outputs
%   individual contribution of each dipole.
%
%   Observation point                   Point(1:3)         
%   Array of dipole moments             DipoleMoment(1:3,1:EdgesTotal) 
%   Array of dipole centers             DipoleCenter(1:3,1:EdgesTotal)
%   E-field at the observation point    E(1;3,1:EdgesTotal)
%   H-field at the observation point    H(1;3,1:EdgesTotal)
%
%   Copyright 2002 AEMM. Revision 2002/03/11 
%   Chapter 3

%C=4*pi;
%ConstantH=K/C;
%ConstantE=eta_/C;
    
m=DipoleMoment;
c=DipoleCenter;
%a = length(c);
%b = c(1:3,:);
r       =repmat(Point,[1 length(c)])-c(1:3,:);
PointRM =repmat(sqrt(sum(r.*r)),[3 1]);
EXP     =exp(-K*PointRM);
PointRM2=PointRM.^2;
C=1./PointRM2.*(1+1./(K*PointRM));
D=repmat(sum(r.*m),[3 1])./PointRM2;
M=D.*r;
HField=ConstantH*cross(m,r).*C.*EXP;
EField=ConstantE*((M-m).*(K./PointRM+C)+2*M.*C).*EXP;

end