function [mappedWaves, first_good_slide, approx_cycle_length]=MapEKG2Color(ekg_signal, L, start_slide, end_slide)

    ekg_signal=ekg_signal(start_slide: end_slide);
    N=length(ekg_signal);   mappedWaves=string(zeros(N,1));
    max_sgnl=max(ekg_signal);    min_sgnl=min(ekg_signal);
    
    Fs=864.0553;
    fps=Fs/L;
    
%     first_good_slide=5500;  
    peak_spacing=500; %Works well for adult subjects
%     peak_spacing=350;
    first_good_slide=1;
    ekg_signal(first_good_slide:end)=ekg_signal(first_good_slide:end);
    ekg_signal(1:first_good_slide-1)=min_sgnl;
    mappedWaves(1:first_good_slide-1)="N";      %Set bad slide to N for noise
    %% Find Biggest EKG Peaks---QRS Complex
    qrs_time=0.1; %Average qrs complex length (in s)
    qrs_frames=floor(qrs_time*Fs);
%     findpeaks(ekg_signal,'MinPeakDistance', 500);
    [vals_qrs, peak_ind_qrs]=findpeaks(ekg_signal,'MinPeakDistance', peak_spacing);
    Delta_frame=zeros(length(peak_ind_qrs)-1,1);
    for i=1:length(peak_ind_qrs)-1
        Delta_frame(i)=peak_ind_qrs(i+1)-peak_ind_qrs(i);
    end
    approx_cycle_length_frame=mean(Delta_frame);    approx_cycle_length=approx_cycle_length_frame*(1/Fs);
    
    for i=1:length(peak_ind_qrs)
        mappedWaves(max(1,peak_ind_qrs(i)-qrs_frames):min(peak_ind_qrs(i)+qrs_frames,N))="qrs";
    end
    
    %Set Peaks to 0
    for i=1:length(peak_ind_qrs)
        ekg_signal(1,max(peak_ind_qrs(i)-qrs_frames,1):min(peak_ind_qrs(i)+qrs_frames,N))=min_sgnl;
    end
    
    %% Find Second Biggest EKG Peaks---T Wave
    t_time=0.17;
    t_frames=floor(t_time*Fs);
    [vals_t, peak_ind_t]=findpeaks(ekg_signal,'MinPeakDistance',peak_spacing);
    for i=1:length(peak_ind_t)
        mappedWaves(max(1,peak_ind_t(i)-t_frames):min(peak_ind_t(i)+t_frames,N))="t";
    end
    %Set Peaks to 0
    for i=1:length(peak_ind_t)
        ekg_signal(1,max(1,peak_ind_t(i)-t_frames):min(peak_ind_t(i)+t_frames,N))=min_sgnl;
    end

    
    
   %% Fix Discontinuities in Waves
   lbl0=mappedWaves(1);
   for i=1:length(mappedWaves)
       lbl=string(mappedWaves(i));
       
       if lbl0=="qrs" && lbl=="0"
           mappedWaves(i)="t";
       elseif lbl=="qrs"
           lbl0="qrs";
       elseif lbl=="t"
           lbl0="t";
       elseif lbl0=="t" && lbl=="0"
           mappedWaves(i)="p";
       end

   end
%    %% Label Rest of Signal--P Wave
%    mappedWaves(ekg_signal~=0)="p";


end