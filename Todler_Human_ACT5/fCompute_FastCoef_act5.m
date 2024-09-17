function [U0, AtA, Fast_Coef]=fCompute_FastCoef_act(R0, h, L, b)

%% Compute the Forward problem, A matrix, and Regularization

%%  Input parameters
unit = 100;           % meters
square_elec = 1;

J = 2;                % number of layers
Sig0 = 1;             % sigma 0
K=L-1;                % number of current pattern
KK=K*K; 
Eh=2.54/unit;   Ew=2.54/unit;   %Electrode Measurements
Zlayer1 = 4+Eh/2;        % distance (cm) from the bottom of the tank to the center of the lowest layer
Zlayer2 = Zlayer1+b - Eh/2;         % distance (cm) from the bottom of the tank to the center of the uppermost layer
% Zlayer1 = Eh/2;        % distance (cm) from the bottom of the tank to the center of the lowest layer
% Zlayer2 = b - Eh/2;    
shift = 0;            % shift(cm) in case the electrodes were moved in the z axis
current_amp = 3.5e-4;   % max current amplitude

DegreesFreedom=496*J;     % Degrees of freedom
Voxel_N=496*J;            % the number of voxels in the mesh element
FN = 16;              % Number of terms in the Fourier coefficient

% Convert to (m)
R0 = R0/unit;
h =  h/unit;             
Eh = Eh/unit;               
Ew = Ew/unit;
Zlayer1 = Zlayer1/unit;        
Zlayer2 = Zlayer2/unit;         
shift = shift/unit;      

if square_elec == 0
    Ea = pi*(Ew/2)^2;
else
    Ea = Eh*Ew;    
end
%% Compute the theta angles and Z 
Z1 = Zlayer1+shift;       % First layer + shift
Z2 = Zlayer2+shift;       % Second layer + shift

Zlz=[repmat(Z1,L/J,1); repmat(Z2,L/J,1)];   % electrode location in z

Theta=(1:1:L/2)*2*pi/(L/J);
Theta_L=repmat(Theta',J,1);                 % electrode location about Theta (radians)  
%% Load the Joshua tree mesh 
load JoshMesh3D_496x2.mat                   % Joshua tree mesh size 496 x 6 (R-, R+, Theta-, Theta+, Z-, Z+)

mz = [Z1-Eh/2, Z1+Eh/2; Z2-Eh/2, Z2+Eh/2];  % mesh for z axis [-z1 +z1, -z2 +z2]

Joshmesh(:,1:2) = Joshmesh(:,1:2)*R0;       % Joshmesh radius goes from 0 to 1
mr = mr*R0;                                 % mesh for the inner radius from 0 to R0

%make z mesh for Joshmesh
Joshmesh(1:Voxel_N/2,5) = mz(1,1);
Joshmesh(1:Voxel_N/2,6) = mz(1,2);

Joshmesh(Voxel_N/2+1:end,5) = mz(2,1);
Joshmesh(Voxel_N/2+1:end,6) = mz(2,2);
%% Load current patterns
load CP32_16x2_M1.mat                           % load the applied current patterns

current_patterns = Cur_pat3D*current_amp;  

CP = current_patterns;                          % current in (A)
 
%% Compute FR matrices
disp('Computing FR...')
tic
[mFR11, mFR12, mFR121, mFR122, mFR141, mFR142, mFR143, mFR144, mFR14, mFR3] = Compute_Matrix_FR(FN,mr,h);
toc
disp('Done!')
%% Compute FZ matrices
disp('Computing FZ...')
tic
[mFZ11, mFZ12, mFZ13, mFZ14, mFZ3] = Compute_Matrix_Fz(FN,mz,h);
toc
disp('Done!')
%% Compute Fourier Coefficients and Itheta integrals
disp('Computing Fourier Coefficients and Itheta integrals...')
tic
[Anm_xk, Bnm_xk, an_xk, bn_xk, ITheta_1, ITheta_2, ITheta_3, ITheta_4] = Compute_FourierCoeff_Itheta(R0,h,Sig0,L,Eh,Ew,Zlz,Theta_L,FN,CP,Joshmesh,Voxel_N, Ea);
toc
disp('Done!')
%% Compute Forward Solution
disp('Computing Forward Solution...')
tic
U0 = forward_gap_model3D(L,R0,h,FN,Theta_L,Zlz,Anm_xk, Bnm_xk, an_xk, bn_xk);
% U0 = U0;
disp('Saving!')
save forward_gap3D_act5_subj15 U0
disp('Done!')
%% ************** Compute the Jacobian matrix A *************
disp('Computing A matrix...')
A_matrix = Compute_A_Matrix(mFR11, mFR12, mFR122, mFR144, mFR14, mFR3,mFZ11, mFZ12, mFZ13, mFZ14, mFZ3,Anm_xk, Bnm_xk, an_xk, bn_xk, ITheta_1, ITheta_2, ITheta_3, ITheta_4,mr,mz,FN,K,Voxel_N,Joshmesh);

disp('Saving...')
save Amatrix_act5_subj100.mat A_matrix
disp('Done!')
%% Regularization
disp('Computing Regularization...')
R=A_matrix;
Matrix_A=reshape(R,KK,DegreesFreedom);
A=Matrix_A;
AtA=A'*A;
% Choose regularization parameters.  There are two
Par_1=0.5;
Par_2=5*10^(-4);
Max_Diag=max(diag(AtA));
% Fast_Coef=inv(AtA+Par_1*diag(diag(AtA))+Par_2*Max_Diag*(eye(size(AtA))))*A' ;
Fast_Coef = (AtA+Par_1*diag(diag(AtA))+Par_2*Max_Diag*(eye(size(AtA))))\A';

disp('Saving...')
save Fast_Coef_act5_subj100.mat Fast_Coef  % Need to come up with a naming that describes the parameters of A or save as a mat file with those parameters
disp('Done!')
