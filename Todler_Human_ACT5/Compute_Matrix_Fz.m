% Program to generate the look-up table for integrals involving z which
% will be used to compute the 3D cylinder geometry A_matrix.

% =========================================================================

% Output: several matrixs - mFZ11, mFZ12, mFZ13, mFZ14 and mFZ3.
% =========================================================================

function [mFZ11, mFZ12, mFZ13, mFZ14, mFZ3] = Compute_Matrix_Fz(N,mz,h)

% mz=[ 8.25, 13.75; 13.75, 19.25;]; % 2 layers in z direction in mesh
% mz = mz + 4;
% N = 32; % Fourier number
% h = 31.4325 ; % height of the cylinder (cm)
Nz = length(mz);

%% %%%%%% compute the matrix mFZ11 %%%%%%%%%%%%%%
mFZ11 = NaN(N,N,Nz);
for iz = 1:Nz
    z1 = mz(iz,1);
    z2 = mz(iz,2);
    mFZ11(:,:,iz) = ones(N,N)*(z2-z1);
   
end
%%

%% %%%%%% compute the matrix mFZ12 and mFZ13%%%%%%%%%%%%%%
mFZ12 = NaN(N,N,N+1,Nz);
for iz = 1:Nz
    z1 = mz(iz,1);
    z2 = mz(iz,2);
        for n = 1:N
            for mp = 1:N
                for np = 0:N
                    mFZ12(n,mp,np+1,iz) = h/(mp*pi)*(sin(mp*pi*z2/h)-sin(mp*pi*z1/h)); 
                end          
            end
        end

end

mFZ13 = mFZ12; 
%%

%% %%%%%% compute the matrix mFZ14 & mFZ3 %%%%%%%%%%%%%%
mFZ14 = NaN(N,N+1,N,N+1,Nz);
mFZ3 = NaN(N,N+1,N,N+1,Nz);
for iz = 1:Nz
    z1 = mz(iz,1);
    z2 = mz(iz,2);
        for m = 1:N
            for n = 0:N
                for mp = 1:N
                    for np = 0:N
                        if m==mp
                           mFZ14(m,n+1,mp,np+1,iz) = h/(4*m*pi)*(sin(2*m*pi*z2/h)-sin(2*m*pi*z1/h))+(z2-z1)/2;
                           mFZ3(m,n+1,mp,np+1,iz) = -h/(4*m*pi)*(sin(2*m*pi*z2/h)-sin(2*m*pi*z1/h))+(z2-z1)/2;
                        else
                           mFZ14(m,n+1,mp,np+1,iz) = h/(2*(m+mp)*pi)*(sin((m+mp)*pi*z2/h)-sin((m+mp)*pi*z1/h))+ ...
                                                     h/(2*(m-mp)*pi)*(sin((m-mp)*pi*z2/h)-sin((m-mp)*pi*z1/h));
                           mFZ3(m,n+1,mp,np+1,iz) = -h/(2*(m+mp)*pi)*(sin((m+mp)*pi*z2/h)-sin((m+mp)*pi*z1/h))+ ...
                                                     h/(2*(m-mp)*pi)*(sin((m-mp)*pi*z2/h)-sin((m-mp)*pi*z1/h));
                        end
                    end
                end          
            end
        end

end


end