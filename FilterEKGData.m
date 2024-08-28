function FilterEKGData(frame_ECG, sbj_num)
   [L, samples, num_frames]=size(frame_ECG);
    
    Fs = 864.0553;                                      %Sampling Frequency (1/s)
   filter_iter=1;

   num_frames=min(num_frames,2000);
   frame_ECG=frame_ECG(:,:, 1:num_frames);

   total_data_pts=samples*num_frames;

%    EKG_save_file="AdaptiveFilter_ECG_subj705.mat";
%    EKG_save_file="AdaptiveFilter_ECG_subj63.mat";
   EKG_save_file="AdaptiveFilter_ECG_subj"+sbj_num+".mat";

   params=[0.1, 1000, 0.05, 1];
   %% Record ECG Wave
   while true
       ecg_wave=zeros(total_data_pts,1);

       ecgwave_save_file="AdaptiveFilter_ecgwave_sbj63.mat";
       fprintf("Calculating EKG Wave...\n");
       [ecg_wave, params]=ECG_Adaptive_Filter_Nilton(frame_ECG, 1, num_frames, 16, 25, true, params);





       [mappedWaves, first_good_slide, approx_cycle_length]=MapEKG2Color(ecg_wave', 32, 1, length(ecg_wave));
       [first_ecg_idx, end_ecg_idx]=isolateCardiacCycle(mappedWaves, 1);  num_ecg_slides=end_ecg_idx-first_ecg_idx+1;

       plot(first_ecg_idx:end_ecg_idx, ecg_wave(first_ecg_idx:end_ecg_idx));

       ecg_OK=input("Is the ECG Filtered Correctly? (y/n):  ", "s");
       
       if strcmp(ecg_OK, "y")
           break;
       elseif strcmp(ecg_OK, "end")
           break;
       else
           fprintf("Parameters:  C=%.2g,  |  filter order=%.2g,  |   step size=%.2g,  |   # Iterations=%d \n", params(1), params(2), params(3), params(4));
           C=input("New C=");
           LMS_order=input("New Order=");
           mu_LMS=input("New Step Size=");
           N_iterat=input("New Number of Iterations=");
           params=[C, LMS_order, mu_LMS, N_iterat];
       end

   end

   save(ecgwave_save_file, 'ecg_wave');
   fprintf("Filtered EKG wave saved to:  "+ecgwave_save_file);


   %% Filter Voltages Recorded on Each Electrode
   Vekg=zeros(L,total_data_pts);
   fprintf("Filtering ECG Data... \n");

   for l=1:L
       linelength=fprintf("Filtering Data on e_l=%d", l);
       [filtered_signal,~]=ECG_Adaptive_Filter_Nilton(frame_ECG, 1, num_frames, l, 0, true, params);

%        filtered_signal=RemoveFreqRange(filtered_signal, Fs, 59, 61);
       Vekg(l,:)=filtered_signal;
       fprintf(repmat('\b',1,linelength))

   end
    fprintf("\n");

   %% Normalize Data
   fprintf("Normalizing Data...")
   for s=1:total_data_pts
        Vekg(:,s)=Vekg(:,s)-sum(Vekg(:,s))/L;
   end
   fprintf("\n ECG data is Finished Filtering");
   
   
   save(EKG_save_file, 'Vekg');
   fprintf("Filtered EKG signal saved to file: "+EKG_save_file +"\n");

   load(EKG_save_file);



end