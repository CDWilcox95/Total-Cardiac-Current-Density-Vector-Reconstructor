function plt_onto_JoshMesh3D(gamma_real_vec, gamma_imag_vec, ecg_wave, query_frames, sigmab, real_title, imag_title, new_fig, sbj_num)

    [M,numSlides]=size(gamma_real_vec);
    load joshmap3D_real_shape_ACT5;
   cd Todler_Human_ACT5; 
   
DegreesFreedom=496*2;
    video=1;
    start_frame=1;  end_frame=numSlides;
    num=zeros(1,31);
den=zeros(1,31);
[M,N]=size(Jash);
Real_image=zeros(M,N); 
Imag_image=zeros(M,N); 

%Circle
% real_vec = zeros(643,302,length(start_frame:end_frame));
% imag_vec = zeros(643,302,length(start_frame:end_frame));

%Conformal Map
real_vec = zeros(229,248,length(start_frame:end_frame));
imag_vec = zeros(229,248,length(start_frame:end_frame));
% sigma_b = zeros(1,length(start_frame:end_frame));


Rmin = min(min(gamma_real_vec))*0.75;
Rmin_vec = zeros(1,length(start_frame:end_frame));
Rmax = max(max(gamma_real_vec))*0.75;
Rmax_vec = zeros(1,length(start_frame:end_frame));
if isempty(gamma_imag_vec)
    Imin=-1e10; Imax=1e10;
else
    Imin = min(min(gamma_imag_vec))*0.75;
    Imax = max(max(gamma_imag_vec))*0.75;

end
Imin_vec = zeros(1,length(start_frame:end_frame));
Imax_vec = zeros(1,length(start_frame:end_frame));
    %% Address the pixels for plotting
    fixed_range = 1;

    
    index = 1;
for S=start_frame:end_frame

    if ~isempty(gamma_real_vec)
    gamma_real = gamma_real_vec(:,index);
    else
        gamma_real=zeros(992,1);
    end

    if length(gamma_imag_vec)>3
        remove2col=false;   plt_idx=1;
        gamma_imag = gamma_imag_vec(:,index);
    elseif isempty(gamma_imag_vec)
        remove2col=true;    plt_idx=1;
        gamma_imag=zeros(992,1);
    elseif length(gamma_imag_vec)==3
        remove2col=true;
        plt_idx=gamma_imag_vec(1);
        % Rmax=0.25*gamma_imag_vec(2);
        % Rmin=0.25*gamma_imag_vec(3);
        gamma_imag=zeros(992,1);
    end
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
              
    Cbar_R(1,:) = ((Range_R).*(1:N))/N+Rmin;    Cbar_R(1,:)=flipud(Cbar_R);
    Cbar_I(1,:) = ((Range_I).*(1:N))/N+Imin;                    

%     colormap( flipud(Cbar_R) );
    Cbar_Rout=[NaN(20,N)*RBackground;repmat(Cbar_R,10,1);Cbar_R;NaN(8,N)*RBackground];  Cbar_Rout=flipud(Cbar_Rout);
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


    %DICOM
%     Real_image(1:302,:)=flip(Real_image(1:302,:));  Real_image(303:end,:)=flip(Real_image(303:end,:));
%     Imag_image(1:302,:)=flip(Imag_image(1:302,:));  Imag_image(303:end,:)=flip(Imag_image(303:end,:));
    Image_Rout = [Cbar_Rout; Real_image];
    Image_Iout = [Cbar_Iout; Imag_image];
    
    real_vec(:,:,index) = Image_Rout;
    imag_vec(:,:,index) = Image_Iout;  
    
    index = index+1;
    
    msg = "Computing Frame: " + string(S);
%     disp(msg);
    
    
end



                
%% Plot Reconstructions
if new_fig
fig1 = figure('color','w','Position',[10 50 900 450]);

end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% max_real = max(max(max(real_vec)));
% min_real = min(min(min(real_vec)));
% 
% max_imag = max(max(max(imag_vec)));
% min_imag = min(min(min(imag_vec)));
% Rmin = min_real;
% Rmax = max_real;
% Imin = min_imag;
% Imax = max_imag;
 file_name="VPidx";
 Ref_frame=1;
if video == 1
    v = VideoWriter("Vid_HeartActivation_sbj"+sbj_num+"_Todler", 'MPEG-4');
    v.FrameRate = 30;
    open(v);
end

index = query_frames(1);
sigma_b = zeros(1,length(start_frame:end_frame));
        strt_idx=1; end_idx=30;

for i=query_frames
   % Conductivity

        clf;
        
        if fixed_range ~= 1

            Rmin = Rmin_vec(index);
            Rmax = Rmax_vec(index);
            Imin = Imin_vec(index); 
            Imax = Imax_vec(index);          

        end        

        Tle=real_title;
        Image_Rout = real_vec(:,:,index);

        XLab=['min: ',num2str(Rmin,'%4.2f'),'     (mS/M)     ', 'max:  ',num2str(Rmax,'%4.2f')];
        subplot(2,2,[1,3]);
%         subplot(8,3,[1,10]);
        surface(Image_Rout);
        setupfigure;                               
        title(Tle,'FontSize',12);
        text(0,-10, XLab,'FontSize',12);
%         caxis([Rmin Rmax]);

    %Susceptance             

    if ~isempty(gamma_imag_vec)
        Tle=imag_title;
        Image_Iout = imag_vec(:,:,index);
        XLab=['min: ',num2str(Imin,'%4.2f'),'         (mS/M)         ', 'max: ',num2str(Imax,'%4.2f')];
        subplot(2,3,2);
        subplot(8,3,[2,11]);
        surface(Image_Iout);
        setupfigure;
        text(0,-10, XLab,'FontSize',12);
        title(Tle,'FontSize',12);
        %         caxis([Imin Imax]);
    end
        
        if ~isempty(ecg_wave)
            if ~isempty(ecg_wave) && ~remove2col
                subplot(2,2,2);
            elseif ~isempty(ecg_wave) && remove2col
                subplot(2,2, 2);
            end
            xlabel("Frame Index");  ylabel("Voltage (mV)"); hold on;
            plot(1:length(ecg_wave),ecg_wave, 'b-'); hold on;
            hold on
            plot(plt_idx,real(ecg_wave(plt_idx)),'-o');
            hold off;
            title("Lead 2 ECG");
        end


        
    %Sigma best imag    
        if ~isempty(sigmab) && ~remove2col
            subplot(2,2,4);
%             xlabel("Frame Index"); ylabel("Volume (mL)"); hold on;
            xlabel("Frame Index"); ylabel("Voltage (mV)"); hold on;

            plot(real(sigmab), 'r-');
            xlim([1 length(sigmab)]);
            hold on
            plot(plt_idx,real(sigmab(plt_idx)),'-bo');
            hold off;
            title("Max Voltage Over Measured Voltages");
%             title("(Averaged) Volume of Blood");
        elseif ~isempty(sigmab) && remove2col
            subplot(2,2, 4);
%             xlabel("Frame Index"); ylabel("Volume (mL)"); hold on;
            xlabel("Frame Index"); ylabel("Voltage (mV)"); hold on;           
            plot(real(sigmab), 'r-');
            xlim([1 length(sigmab)]);
            hold on
            plot(plt_idx,real(sigmab(plt_idx)),'-bo');
            hold off;
            title("Max Voltage Over Measured Voltages");
%             title("(Averaged) Volume of Blood");      
        end
        drawnow limitrate; 
        msg = "Plotting Frame: " + string(i);
        disp(msg);    

        if video == 1

            frame = getframe(gcf);
            writeVideo(v,frame);

        end
    plt_idx=i;
    index = i;
end

close(v);
cd ../;
end