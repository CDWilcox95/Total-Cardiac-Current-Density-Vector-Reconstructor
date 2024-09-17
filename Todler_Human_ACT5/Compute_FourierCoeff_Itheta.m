% pre-compute F_theta lookuptable for an_xk,bn_xk,Anm_xk,A0m_xk,Bnm_xk,B0m_xk,terms;
% Pass parameters - made this a function
% This is for two rows of 16 electrodes 

function [Anm_xk, Bnm_xk, an_xk, bn_xk, ITheta_1, ITheta_2, ITheta_3, ITheta_4] = Compute_FourierCoeff_Itheta(R0,H,Sig0,L,Eh,Ew,Zlz,Theta_L,FN,CP,Joshmesh,Voxel_N, Ea)

%  Electrode parameters

XK=L-1;      % number of the current patterns
% Eh=5.3975 ;      % electrode height (cm)
% Ew=5.3975;       % electrode width (cm)
% Ea=Eh*Ew;   % area of an electrode (assumes square)

% Zlz=[repmat(15,16,1); repmat(20.5,16,1); ];   % Distance: from the bottom of the tank to the center of the electrodes on the lzth layer (***)
% Theta=[0:1:15]*2*pi/16;
% Theta_L=repmat(Theta',2,1);    % electrode location about Theta (radians)        
Delta_E=Ew/R0;            % The change in theta angle Theta Theta_l-Deta_E/2 <= Theta <= Theta_l+Deta_E/2

% FN=128;      % number of Fourier basis functions: okay to change

%% Load Current pattern
% load CP32_16x2_M8 ;
% CP=(CP32_16x2_M8/8);


%% ************** pre-compute an_xk,bn_xk,Anm_xk,A0m_xk,Bnm_xk,B0m_xk,term *************

an_xk=zeros(FN+1, XK);
bn_xk=zeros(FN+1, XK);
Anm_xk=zeros(FN+1, FN, XK);
Bnm_xk=zeros(FN+1, FN, XK);


for jxk=1:XK  % loop over the current patterns
    
    CP_xk=CP(:, jxk);
    
    for n=1:FN  % loop over the Fourier coefficients to be computed
        an=0;
        bn=0;           % compute an's and bn's
        for i=1:L
            an=an+CP_xk(i)*cos(n*Theta_L(i));
            bn=bn+CP_xk(i)*sin(n*Theta_L(i));
        end
    
        Const_term=2*Eh*sin(n*Delta_E/2)/(Ea*Sig0*pi*(n^2)*H*R0^(n-1));
        an_xk(n+1,jxk)=Const_term*an;
        bn_xk(n+1,jxk)=Const_term*bn;  
    end
    
    %% ************ compute A0m*In*cos(m*(pi)z/H, where In is the modified Bessel function of the first kind, besseli *************

    for m=1:FN
        A0m=0;
        for i=1:L
            A0m=A0m+CP_xk(i)*cos(m*pi*Zlz(i)/H);
        end
        Bessel_term=besseli(1,(m*pi*R0/H));
        Const_term=2*H*Delta_E*sin((m*pi*Eh)/(2*H))/(Ea*Sig0*pi^3*m^2*Bessel_term);
        Anm_xk(1,m,jxk)=A0m*Const_term;
    end


%% ************ compute Anm*In*cos(m*pi*z/H) & Bnm*In*cos(m*pi*z/H) term

    for m=1:FN
        for n=1:FN
            Anm=0;
            Bnm=0;
            for i=1:L
                Anm=Anm+CP_xk(i)*cos(m*pi*Zlz(i)/H)*cos(n*Theta_L(i));
                Bnm=Bnm+CP_xk(i)*cos(m*pi*Zlz(i)/H)*sin(n*Theta_L(i));
            end
            Bessel_term=besseli(n+1,(m*pi*R0/H)) + besseli(n,(m*pi*R0/H)) * (n*H)/(m*pi*R0);
            Const_term=8*H*sin(n*Delta_E/2)*sin(m*pi*Eh/H/2)/(Ea*Sig0*pi^3*m^2*n*Bessel_term);
            Anm_xk(n+1,m,jxk)=Anm*Const_term;
            Bnm_xk(n+1,m,jxk)=Bnm*Const_term;
        end
    end
end  % jxk loop


%% Loading the Joshua tree mesh 
% load 3DMesh496x2 ; % Joshua tree mesh size 496 x 6 (R-, R+, Theta-, Theta+, Z-, Z+)
% Voxel_N=496*2 ;          % the number of the mesh element

ITheta_1=zeros(Voxel_N, FN+1, FN+1);
ITheta_2=zeros(Voxel_N, FN+1, FN+1);
ITheta_3=zeros(Voxel_N, FN+1, FN+1);
ITheta_4=zeros(Voxel_N, FN+1, FN+1);

for voxel=1:Voxel_N
    th1=Joshmesh(voxel,3);
    th2=Joshmesh(voxel,4);
    
    for n=0:FN
        for np=0:FN
            if ((n==np) && (n==0))
                ITheta_1(voxel, n+1, np+1)=th2-th1;
                ITheta_2(voxel, n+1, np+1)=0;
                ITheta_3(voxel, n+1, np+1)=0;
                ITheta_4(voxel, n+1, np+1)=0;
            elseif ((n==np) && (n~=0))
                ITheta_1(voxel, n+1, np+1)=(sin(2*n*th2) - sin(2*n*th1))/(4*n) + (th2-th1)/2 ;
                ITheta_2(voxel, n+1, np+1)=(th2-th1)/2 + (sin(2*n*th1)-sin(2*n*th2))/(4*n) ;
                ITheta_3(voxel, n+1, np+1)=(cos(2*n*th1)-cos(2*n*th2))/(4*n);
                ITheta_4(voxel, n+1, np+1)=(cos(2*np*th1)-cos(2*np*th2))/(4*np);
            elseif (n~=np)
%                 ITheta_1(voxel, n+1, np+1)=(sin((n+np)*th2)-sin((n+np)*th1))/(2*(n+np)) + (sin((np-n)*th2)-sin((np-n)*th1))/(2*(np-n)) ;
                ITheta_1(voxel, n+1, np+1)=(sin((n+np)*th2)-sin((n+np)*th1))/(2*(n+np)) + (sin((n-np)*th2)-sin((n-np)*th1))/(2*(n-np)) ;
%                 ITheta_2(voxel, n+1, np+1)=(sin((n+np)*th1)-sin((n+np)*th2))/(2*(n+np)) + (sin((np-n)*th2)-sin((np-n)*th1))/(2*(np-n));
                ITheta_2(voxel, n+1, np+1)=(sin((n+np)*th1)-sin((n+np)*th2))/(2*(n+np)) + (sin((n-np)*th2)-sin((n-np)*th1))/(2*(n-np));
%                 ITheta_3(voxel, n+1, np+1)=(cos((n+np)*th1)-cos((n+np)*th2))/(2*(n+np)) + (cos((np-n)*th1)-cos((np-n)*th2))/(2*(np-n));
                ITheta_3(voxel, n+1, np+1)=(cos((n+np)*th1)-cos((n+np)*th2))/(2*(n+np)) + (cos((n-np)*th2)-cos((n-np)*th1))/(2*(n-np));
%                 ITheta_4(voxel, n+1, np+1)=(cos((n+np)*th1)-cos((n+np)*th2))/(2*(n+np)) + (cos((n-np)*th1)-cos((n-np)*th2))/(2*(n-np));
                ITheta_4(voxel, n+1, np+1)=(cos((n+np)*th1)-cos((n+np)*th2))/(2*(n+np)) + (cos((np-n)*th2)-cos((np-n)*th1))/(2*(np-n));
            else
                fprintf('Error')


            end

%             if ((n==np) && (n==0))
%                 ITheta_1(voxel, n+1, np+1)=th1-th2;
%                 ITheta_2(voxel, n+1, np+1)=0;
%                 ITheta_3(voxel, n+1, np+1)=0;
%                 ITheta_4(voxel, n+1, np+1)=0;
%             elseif ((n==np) && (n~=0))
%                 ITheta_1(voxel, n+1, np+1)=(sin(2*n*th1) - sin(2*n*th2))/(4*n) + (th1-th2)/2 ;
%                 ITheta_2(voxel, n+1, np+1)=(th1-th2)/2 - (sin(2*n*th1)-sin(2*n*th2))/(4*n) ;
%                 ITheta_3(voxel, n+1, np+1)=(cos(2*n*th1)-cos(2*n*th2))/(4*n);
%                 ITheta_4(voxel, n+1, np+1)=(cos(2*n*th2)-cos(2*n*th1))/(4*n);
%             elseif (n~=np)
%                 ITheta_1(voxel, n+1, np+1)=(sin((n+np)*th1)-sin((n+np)*th2))/(2*(n+np)) + (sin((n-np)*th1)-sin((n-np)*th2))/(2*(n-np)) ;
%                 ITheta_2(voxel, n+1, np+1)=(-1)*(sin((n+np)*th1)-sin((n+np)*th2))/(2*(n+np)) + (sin((n-np)*th1)-sin((n-np)*th2))/(2*(n-np));
%                 ITheta_3(voxel, n+1, np+1)=(cos((n+np)*th1)-cos((n+np)*th2))/(2*(n+np)) + (cos((n-np)*th1)-cos((n-np)*th2))/(2*(n-np));
%                 ITheta_4(voxel, n+1, np+1)=(cos((n+np)*th2)-cos((n+np)*th1))/(2*(n+np)) - (cos((n-np)*th1)-cos((n-np)*th2))/(2*(n-np));
%             else
%                 fprintf('Error')
%             end
        end
    end
end

end



