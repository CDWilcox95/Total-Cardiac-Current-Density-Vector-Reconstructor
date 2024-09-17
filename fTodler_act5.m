function [gamma_recon_real, gamma_recon_imag, sigma_b]=fTodler_act5(frame_voltage, cur_pattern, model_info, new_sbj_param, plt_recon)

%% Human Data - Todler reconstructions for ACT 5 

% Nilton Rosa, Jan 2022.

close all;
% clear all;
format long;
%%  Input parameters
[L,K,end_frame]=size(frame_voltage);
% end_frame=min(end_frame, 2000);
L=32;                   % number of electrodes on the surface of the body
K=L-1;                  % number of current pattern

fixed_range = 1;
video = 1;              % 1 - Make a video, 0 - don't
scale = 0.03;           % scaling factor for the range
abs_imag = 0;           % 1 - absolute image, 0 - diff image
start_frame =1;         % start frame
gamma1 = 5;             % gamma for "noser" regularization
gamma2 = 5e-4;          % gamma for Tikhonov regularization


%% Load Dataset
% path = "../ACT5_Data/22_06_21";
% file_name = "/sbj7_11kHz_3D_22_06_21_15_18_54";
% 
% path = path + "\" + file_name;
% load(path);

%% Assign Current Pattern
Trig_Pat = cur_pattern(:,1:L-1);
cd Todler_Human_ACT5;
load CP32_16x2_M1.mat                           % load the applied current patterns
current_amp=3.5*10^-4;
current_patterns = Cur_pat3D*current_amp;  

load JoshMesh3D_496x2;  Joshmesh(:,1:2)=model_info.R.*Joshmesh(:,1:2); 
Joshmesh(1:496,5:6)=[zeros(496,1),(model_info.H/2).*ones(496,1)];   JT_mesh_vox(497:992,5:6)=[(model_info.H/2).*ones(496,1),model_info.H.*ones(496,1)];

%% Load Direct Problem Solution

U0save_str="Fwd_AveGap_3DModel_sbj"+model_info.sbj_num+".mat";
Asave_str="Amatrix_sbj"+model_info.sbj_num+".mat";
if new_sbj_param
    fprintf("Running FastCoef algorithm for predicted voltages, U0, and matrix A ... \n");
    [U, A, FastCoef]=fCompute_FastCoef_act5(model_info.R, model_info.H, model_info.num_elec, model_info.vert_gap);
    save(U0save_str, 'U');
    fprintf("Computed U. Saved to: " +U0save_str+"\n");
    
    save(Asave_str, 'A');
    fprintf("Computed A.  Saved to: "+Asave_str+"\n");

    
else
    load(U0save_str);   load(Asave_str);
end




cd ../;
[X, basis_coef]=ChangeBasis(current_patterns, cur_pattern(:,1:L-1));
U=ChangeVoltageCurrent(basis_coef, U);


DegreesFreedom=496*2;
AtA=A'*A;
% Choose regularization parameters.  There are two
Par_1= 0.5;
Par_2= 5*10^-4;
Max_Diag=max(diag(AtA));
Fast_Coef=inv(AtA+Par_1*diag(diag(AtA))+Par_2*Max_Diag*(eye(size(AtA))))*A' ;

% [Usvd, Ssvd, Vsvd]=svd(AtA+Par_1*diag(diag(AtA)));


%% ********** load Joshua mapping file *******
load Todler_Human_ACT5/joshmap3D_circle;

%% Find the reference frame based on the maximum value of sigma best

sigma_b = zeros(1,length(start_frame:end_frame));

index = 1;
for i = start_frame:end_frame
    
    Data_Volts = frame_voltage(1:L-1,:,i);
    Data_Volts = Data_Volts.';  

    Data_Volts = Data_Volts - repmat(sum(Data_Volts)/L,L,1);
    
    V_ell=Data_Volts;

    num = sum(V_ell.*U);
    den = sum(U.^2);

    rhobest=sum(num)/sum(den);  
    gammabest=(1/rhobest);                   %convert to mS/M
    sigma_b(index)= gammabest;
    
    index = index + 1;
end

[max_val, Ref_frame] = min(real(sigma_b));
% Ref_frame=1695;
Ref_frame=1702;
% Ref_frame = 200;          % reference frame
%% Compute ABS image for Reference frame

Data_Volts = frame_voltage(1:L-1,:,start_frame+Ref_frame-1);

Data_Volts = Data_Volts.';  

Data_Volts = Data_Volts - repmat(sum(Data_Volts)/L,L,1);

V_ell=Data_Volts;

num = sum(V_ell.*U);
den = sum(U.^2);

rhobest=sum(num)/sum(den);  
gammabest=(1/rhobest);                   %convert to mS/M

U_ell_gammabest= rhobest*U;

n=0;
for x=1:K
    for k=1:K
      n=n+1;
%       D_todler(n,:) = Trig_Pat(:,x)'*U_ell_gammabest(:,k) - Trig_Pat(:,k)'*V_ell(:,x);
        D_todler(n,:) = U_ell_gammabest(:,x)'*Trig_Pat(:,k) - V_ell(:,k)'*Trig_Pat(:,x);

    end
end     
% [reg_corner,rho,eta,reg_param] = l_curve(Usvd,diag(Ssvd),D_todler,'tikh');

eta = (1/rhobest).*Fast_Coef * D_todler * scale;   %convert to mS/M
% [l1, l2]=Todler2ParamLCurve(A, D_todler);


gammavector_ref = (gammabest + eta);         %absolute reconstruction of ref frame
% gammavector_ref = (eta);         %absolute reconstruction of ref frame

%% initialize variables
num=zeros(1,31);
den=zeros(1,31);
[M,N]=size(Jash);
Real_image=zeros(M,N); 
Imag_image=zeros(M,N); 

real_vec = zeros(643,302,length(start_frame:end_frame));
imag_vec = zeros(643,302,length(start_frame:end_frame));
sigma_b = zeros(1,length(start_frame:end_frame));
gamma_real_vec = zeros(DegreesFreedom, length(start_frame:end_frame));
gamma_imag_vec = zeros(DegreesFreedom, length(start_frame:end_frame));

Rmin = 10e10;
Rmin_vec = zeros(1,length(start_frame:end_frame));
Rmax = -10e-10;
Rmax_vec = zeros(1,length(start_frame:end_frame));
Imin = 10e10;
Imin_vec = zeros(1,length(start_frame:end_frame));
Imax = -10e-10;
Imax_vec = zeros(1,length(start_frame:end_frame));
%% Compute todler for all frames
index = 1;
for S=start_frame:1:end_frame
    
    
    Data_Volts = frame_voltage(1:L-1,:,S);
   
    Data_Volts = Data_Volts.';     

    Data_Volts = Data_Volts - repmat(sum(Data_Volts)/L,L,1);
                     
    V_ell=Data_Volts;
                    
    num = sum(V_ell.*U);
    den = sum(U.*U);                
                    
    rhobest=sum(num)/sum(den);  % best resistivity fit
    gammabest=(1/rhobest);  % best conductivity fit
    sigma_b(index)= gammabest;
    U_ell_gammabest= rhobest*U;
    
    n=0;
    for x=1:K
        for k=1:K
          n=n+1;
%           D_todler(n,:) = Trig_Pat(:,x)'*U_ell_gammabest(:,k) - Trig_Pat(:,k)'*V_ell(:,x);
            D_todler(n,:) = (U_ell_gammabest(:,x)'*Trig_Pat(:,k) - V_ell(:,k)'*Trig_Pat(:,x));

        end
    end      


    eta = (1/rhobest).*Fast_Coef * D_todler * scale;      %convert to mS/M

    % [eta, hist_obj, hist_res] = alm_qp(gammabest.*AtA,c,zeros(size(A)),D_todler,10^-5,beta,gammabest.*ones(992,1));
    % l2_err=norm(A*eta-D_todler)              
    if abs_imag == 1
        gammavector = gammabest+eta;            %absolute image
    else
        gammavector =(gammabest + eta) - gammavector_ref;   % Diff of 2 absolute images
    end

    % Calculate the conductivity vector as the final output.
    gamma_recon_real(:,index)=real(gammavector);
    gamma_recon_imag(:,index)=imag(gammavector);
    sigma_b(index)=gammabest;

    gamma_real = real(gammavector);
    gamma_imag = imag(gammavector);

    gamma_real_vec(:,index) =  gamma_real;
    gamma_imag_vec(:,index) =  gamma_imag;    
    
    Rmin_aux = min(gamma_real);
    Rmax_aux = max(gamma_real);
    Imin_aux = min(gamma_imag);
    Imax_aux = max(gamma_imag);    
    
    if fixed_range == 1
    
        if Rmin_aux < Rmin
           Rmin = Rmin_aux; 
        end

        if Rmax_aux > Rmax
           Rmax = Rmax_aux; 
        end

        if Imin_aux < Imin
           Imin = Imin_aux; 
        end

        if Imax_aux > Imax
           Imax = Imax_aux; 
        end
    
    else
        
        Rmin_vec(index) = Rmin_aux;
        Rmax_vec(index) = Rmax_aux;  
        Imin_vec(index) = Imin_aux; 
        Imax_vec(index) = Imax_aux;  
    
    end
    
    index = index+1;
end

%% Address the pixels for plotting
index = 1;
for S=start_frame:1:end_frame

    gamma_real = gamma_real_vec(:,index);
    gamma_imag = gamma_imag_vec(:,index);
    
    if fixed_range ~= 1
       
        Rmin = Rmin_vec(index);
        Rmax = Rmax_vec(index);
        Imin = Imin_vec(index); 
        Imax = Imax_vec(index);          
        
    end    
    
    Range_R= Rmax - Rmin;
    Range_I= Imax - Imin;
                    
    RBackground=Rmin-0.005*Range_R;
    IBackground=Imin-0.005*Range_I;
                    
    Real_image(:,:) = RBackground;
    Imag_image(:,:) = IBackground;      

    Real_image(Jash > 0) = gamma_real(Jash(Jash > 0));
    Imag_image(Jash > 0) = gamma_imag(Jash(Jash > 0));    
              
    Cbar_R(1,:) = ((Range_R).*(1:N))/N+Rmin;
    Cbar_I(1,:) = ((Range_I).*(1:N))/N+Imin;                    

    Cbar_Rout=[NaN(20,N)*RBackground;repmat(Cbar_R,10,1);Cbar_R;NaN(8,N)*RBackground];
    Cbar_Iout=[NaN(20,N)*IBackground;repmat(Cbar_I,10,1);Cbar_I;NaN(8,N)*IBackground];

%     %% blurring the image, smoothing % image loww pass filter takes about 15 ms  (From Felix's program)
%     
%     PSF = fspecial('disk', 1.5);
%     Blurred = imfilter(Real_image,PSF,'replicate','conv');
%     Real_image=Blurred;
%     Blurred = imfilter(Imag_image,PSF,'symmetric','conv');
%     Imag_image=Blurred;    
       
    Real_image(Jash == 0) = NaN;
    Imag_image(Jash == 0) = NaN;

    Image_Rout = [Cbar_Rout; Real_image];
    Image_Iout = [Cbar_Iout; Imag_image];
    
    real_vec(:,:,index) = Image_Rout;
    imag_vec(:,:,index) = Image_Iout;  
    
    index = index+1;
    
    msg = "Computing Frame: " + string(S);
    disp(msg);
    
    
end
                
%% Plot Reconstructions
if plt_recon
    fig1 = figure('color','w','Position',[10 50 1400 720]);

    % max_real = max(max(max(real_vec)));
    % min_real = min(min(min(real_vec)));
    %
    % max_imag = max(max(max(imag_vec)));
    % min_imag = min(min(min(imag_vec)));
    % Rmin = min_real;
    % Rmax = max_real;
    % Imin = min_imag;
    % Imax = max_imag;
    file_name="sbj7_supine_134kHz";
    if video == 1
        v = VideoWriter(file_name + "_" + string(Ref_frame) + "_" + string(start_frame)+ "_" + string(end_frame));
        v.FrameRate = 30;
        open(v);
    end

    index = 1;

    for i=start_frame:1:end_frame
        % Conductivity

        clf;

        if fixed_range ~= 1

            Rmin = Rmin_vec(index);
            Rmax = Rmax_vec(index);
            Imin = Imin_vec(index);
            Imax = Imax_vec(index);

        end
        Tle = ['Conductivity Image, Sigma = ', num2str(real(sigma_b(index)),'%4.1f'), ' mS/M'];

        Image_Rout = real_vec(:,:,index);

        XLab=['min: ',num2str(Rmin,'%4.2f'),'     (mS/M)     ', 'max:  ',num2str(Rmax,'%4.2f')];
        subplot(2,3,[1,4]);
        surface(Image_Rout);
        setupfigure;
        title(Tle,'FontSize',12);
        text(0,-10, XLab,'FontSize',12);
        %         caxis([Rmin Rmax]);

        %Susceptance

        Tle=['Susceptance Image, Sigma = ', num2str(imag(sigma_b(index)),'%4.1f'), ' mS/M'];
        Image_Iout = imag_vec(:,:,index);
        XLab=['min: ',num2str(Imin,'%4.2f'),'         (mS/M)         ', 'max: ',num2str(Imax,'%4.2f')];
        subplot(2,3,[2,5]);
        surface(Image_Iout);
        setupfigure;
        text(0,-10, XLab,'FontSize',12);
        title(Tle,'FontSize',12);
        %         caxis([Imin Imax]);

        %Sigma best real
        subplot(2,3,3);
        plot(real(sigma_b));
        xlim([1 length(sigma_b)]);
        hold on
        plot(index,real(sigma_b(index)),'-o');
        hold off;
        title("Sigma Best - Real");

        %Sigma best imag
        subplot(2,3,6);
        plot(imag(sigma_b));
        xlim([1 length(sigma_b)]);
        hold on
        plot(index,imag(sigma_b(index)),'-o');
        hold off;
        title("Sigma Best - Imag");

        drawnow limitrate;
        msg = "Plotting Frame: " + string(i);
        disp(msg);

        if video == 1

            frame = getframe(gcf);
            writeVideo(v,frame);

        end

        index = index + 1;
    end

    if video == 1
        close(v);
    end
end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                


