clear; close all;
load joshmap3D_circle 

%% Load in Data and set sbj_num variable
% load("D:\Cardiac_Imaging\Data\Sbj708_23_09_25\Sbj708_post_23_09_25_16_44_05.mat");   sbj_num="708";
load("D:\Cardiac_Imaging\Data\Subject63_22_06_29\sbj63_111kHz_supine_22_06_29_15_54_17.mat"); sbj_num="63"; 
% load("D:\Cardiac_Imaging\Data\Subject17_22_07_07\sbj17_111kHz_prone_22_07_07_14_32_24.mat");    sbj_num="17";
% load("D:\Cardiac_Imaging\Data\Subject7_22_07_11\Sbj7_3D_22_11_07_15_38_54.mat");    sbj_num="7";
% load("D:\Cardiac_Imaging\Data\Sbj705_post_23_09_06_15_42_02.mat");  sbj_num="705";
% load("D:\Cardiac_Imaging\Data\Subject23_22_10_24\Sbj23_supine_22_10_24_15_34_32.mat");  sbj_num="23";
% load("D:\Cardiac_Imaging\Data\Subject42_22_01_21\Sbj42_prone_22_10_21_15_13_52.mat"); sbj_num="42";
% load("D:\Cardiac_Imaging\Data\Subject333_22_10_24\Sbj333_prone_22_10_24_16_35_57.mat"); sbj_num="333";

%% Sets paramaters for program execution
% Construct new FEM mesh to compute reconstruction 
new_test_map=false;     
% Use Adaptive Filter to filter ECG voltages
filter_ECG=false;
% Compute LI, LII, LIII Electrocardiograms from data [default(false) = LII]
all_leads=false;
% Compute conductivity reconstructions using TodLeR
recon_conductivity=false;
% Interpolate Partial Body Surface Map (pBSM) and display
new_BSM_interp=false;


%% Plotting options
video_plt=false;    %Creates a video of TCV over given frames
plt_all=false;      %Plots all cardiac information from data


%% Set model paramaters (cm)

%Set vertical_gap between electrode centers (z-coordinate) [some old data files do not contain vertical_gap variables]
if ~exist('vertical_gap', 'var')
    vertical_gap=7;
end
R=2.54*circumference/(2*pi);    H=vertical_gap+2*4;
Fs=864.0553;    %Sampling Frequency




%% Filter ECG Data Using Adaptive Filter (Filter is performed over each electrode)
if filter_ECG
    FilterEKGData(frame_ECG, sbj_num); 
end
load("AdaptiveFilter_ECG_subj"+sbj_num+".mat");

%% Isolate data frames used in analysis [chosen far enough from initial data collection so that data is clean]
% s0ekg=1500*32;  sNekg=1900*32; %sbj 63 & 17
s0ekg=1500*32;  sNekg=1746*32;    %sbj 7
% s0ekg=100*32;   sNekg=300*32;    %sbj 708

%% Converts voltages from (mV) to (V)
Vekg=(10^-3).*Vekg(:,s0ekg:sNekg);


%% Compute and Display ECG lead for chosen frames [s0ekg:sNekg]
if all_leads
    %Compute and plot Lead 2 (L2) Electrocardiogram
    L3_ecg=Vekg(32,:)-Vekg(9,:);
    L2_ecg=Vekg(25,:)-Vekg(16,:);
    L1_ecg=Vekg(25,:)-Vekg(32,:);
    aVR=Vekg(25,:); aVL=Vekg(32,:); aVF=Vekg(16,:);
    check_ecg=L1_ecg+L3_ecg;
else
    L2_ecg=Vekg(25,:)-Vekg(16,:);
end

%% Map each sample to P, QRS, or T wave and truncate data to number of cardiac cycles being analyzed
num_cycles=1;  % Number of cardiac cycles analyzed
% Q0=zeros(1,3);  Pt=model_info.FEM_Mesh.nodes;   elec_nodes=model_info.FEM_Mesh.electrode_nodes; elec_pts=Pt(elec_nodes,:);
[mappedWaves, first_good_slide, approx_cycle_length]=MapEKG2Color(L2_ecg, 32, 1, length(L2_ecg));
[first_ecg_idx, end_ecg_idx, int_cycles]=isolateCardiacCycle(mappedWaves,num_cycles);  num_ecg_slides=end_ecg_idx-first_ecg_idx; cardiac_length=num_ecg_slides/Fs;  

%Sets sample range for plotting
event_frames=first_ecg_idx:end_ecg_idx;

figure;   plt_scale=10^3;
if all_leads
    tiledlayout(3,1);
    title("Bipolar Limb Leads");
    nexttile(1);    plot(plt_scale.*L1_ecg(event_frames));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated Lead I Recording"); grid on;
    nexttile(2);    plot(plt_scale.*L2_ecg(event_frames));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated Lead II Recording"); grid on;
    nexttile(3);    plot(plt_scale.*L3_ecg(event_frames));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated Lead III Recording"); grid on;
    
    figure;tiledlayout(3,1);
    title("Unipolar Limb Leads");
    nexttile(1);    plot(plt_scale.*aVR(first_ecg_idx:end_ecg_idx));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated aVR Recording"); grid on;
    nexttile(2);    plot(plt_scale.*aVL(first_ecg_idx:end_ecg_idx));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated aVL Recording"); grid on;
    nexttile(3);    plot(plt_scale.*aVF(first_ecg_idx:end_ecg_idx));    xlabel("frame");    ylabel("voltage (mV)");  title("Approximated aVF Recording"); grid on;
else
    title("Approximated Lead II Recording");
    plot(plt_scale.*L2_ecg(event_frames));    xlabel("frame");    ylabel("voltage (mV)"); grid on;
    
end

%% Compute Conductivity Reconstructions 
if recon_conductivity
    % fTodler_act5 returns conductivity (gamma_real_vec) permittivity
    % (gamma_image_vec) and best-fitted constant conductivity (sigma_b)
    [gamma_real_vec, gamma_imag_vec, sigma_b]=fTodler_act5(frame_voltage(:,:,1:2000), cur_pattern, false);
    sigmab_R=real(sigma_b);
    sigma_b=zeros(length(sigma_b),1);
    for i=1:length(sigmab_R)-1
        sigma_b((i-1)*32+1:i*32+32)=sigmab_R(i);
    end
    save("sbj"+sbj_num+"_sigma0.mat", 'sigma_b');
end

load("sbj"+sbj_num+"_sigma0.mat");

%% Isolate quadrants of Joshua Tree mesh and prepare for use in source analysis
[LL_lung_voxels, LR_lung_voxels, UL_lung_voxels, UR_lung_voxels, ...
          LA_lung_voxels, LP_lung_voxels, UA_lung_voxels, UP_lung_voxels, plot_mesh] =define_lung_lobes(Jash);
     
UAL=intersect(UL_lung_voxels, UA_lung_voxels); %Upper anterior region
LAL=intersect(LL_lung_voxels, LA_lung_voxels); %Lower anterior region

% Load in JT mesh and get midpoint of each voxel
load("mesh_ele.dat");
load("D:\Cardiac_Imaging\Todler_Human_ACT5\JoshMesh3D_496x2.mat");
mesh_idx=1:992; top_vox=1:496;  bot_vox=497:992;
%Scale Joshua Tree mesh to subjects radius (R) 
JT_mesh_vox=Joshmesh;   JT_mesh_vox(:,1:2)=R.*JT_mesh_vox(:,1:2);   JT_mesh_vox(top_vox,5:6)=[zeros(496,1),(H/2).*ones(496,1)];   JT_mesh_vox(bot_vox,5:6)=[(H/2).*ones(496,1),H.*ones(496,1)];
mid_pts=GetVoxelCenters(JT_mesh_vox);

%% Set and Define Test Region
sub_vox=DefineSubRegion(JT_mesh_vox);
test_region=mesh_idx(sub_vox);

test_pts=(10^-2) .*mid_pts(test_region,:);
plt_seg=zeros(992,1);   plt_seg(test_region)=1; 
% non_test_region=setxor(mesh_idx, Heart);  plt_seg(non_test_region,:)=nan;
% plt_onto_JoshMesh3D(plt_seg, [], L2_ecg(event_frames(1:length(event_frames))), 1, [], "Sample Region", "", true, sbj_num); 




%% Setup model_info paramater w/ subject info  (info in meters)
model_info.R=10^-2 *R;   model_info.H=10^-2 *H;   Bh(4)=model_info.R/3; model_info.num_elec=32; model_info.shape='c'; model_info.elec_planes=[0.04, 0.04+vertical_gap];
fmdl=GetReconMesh(model_info.R,model_info.H, new_test_map, "const", 'c', []);   model_info.FEM_Mesh=fmdl;   elec_pts=fmdl.nodes(fmdl.electrode_idx,:); 

%% Finds best source point in test region 
num_recon_src=1;    Bh(1)=Bh(1)+0.02; Bh(2)=Bh(2)+0.02;
[Qc, multi_Qc, M0, src_vox, frame_errs, frame_idx, mesh_errs_t]=FindBestQ0([], Bh, test_pts, test_region, elec_pts, model_info, Vekg(:, event_frames), "FitFwdSims","FEM", new_test_map, num_recon_src, "sbj"+sbj_num+"_data_");      save("BestFitQc_multi"+num2str(num_recon_src)+".mat", 'multi_Qc'); save("BestFitQ0c"+num2str(num_recon_src)+".mat", 'Qc');
save("DipoleSrc_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'Qc'); save("rTCV_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'M0');
save("MeshErrorIntensity_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'mesh_errs_t');
save("SourceVoxels_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'src_vox');

load("DipoleSrc_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'Qc'); load("rTCV_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'M0');
load("MeshErrorIntensity_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'mesh_errs_t');
load("SourceVoxels_sbj"+sbj_num+"_K"+num2str(num_recon_src)+".mat", 'src_vox');


% Sets specific frames (typically peak wave frames)
query_frames=[447, 513, 545, 578, 756]; %sbj 63


%Plots voxel(s) where source is located onto Joshua Tree mesh (red) inside
%test region (green)
plt_src=zeros(length(mesh_idx),1); plt_src(test_region)=1;  plt_src(src_vox)=2;
plt_onto_JoshMesh3D(plt_src, [], [],1, [], "Test Region w/ Dipole Sources", "", true, sbj_num); 


% %% Generate Movie
% if num_recon_src>1
%     MultiSourcePlotter(sigma_b(event_frames), Qc, M0, mappedWaves, event_frames, L2_ecg, video_plt);
% else
%     TCVPlotter(sigma_b(event_frames), Qc, M0, mappedWaves, event_frames, L2_ecg, video_plt)
% end


%% Set threshold and look at set of voxels with smallest l2 error over time
prcnt_t=0.3; 
Heart=zeros(992, length(event_frames)); sub_voxels=mesh_idx(sub_vox);
plt_seg=zeros(992,length(event_frames));   Mc=zeros(length(event_frames),3);
for t=1:length(event_frames)
    range=max(mesh_errs_t(:,t))-min(mesh_errs_t(:,t));
    thresh=min(mesh_errs_t(:,t))+prcnt_t*range;
    voxel_region=mesh_errs_t(:,t)<thresh;
    sub_vox_region=sub_voxels(voxel_region);
    plt_seg(sub_vox_region,t)=mesh_errs_t(voxel_region,t);
end

plt_onto_JoshMesh3D(plt_seg(:, 1:length(event_frames)), [], L2_ecg(event_frames(1:length(event_frames))),query_frames, [], "UAL", "", false, sbj_num); 


%% Compute pBSM by interpolating measured voltages on each electrode over a cylinder
% This takes a long time to compute, needs to be fixed (8/27/24)
if new_BSM_interp
    Vq=InterpolateVekgBSM(model_info, Vekg, true);  save("interpV_sbj"+sbj_num+".mat", 'Vq', '-v7.3');
end
% load("interpV_sbj"+sbj_num+".mat");


%% Plot Total Cardiac Vector

TCVicons=["+", "*", "o", "hexagram", "diamond"];

figure;
tiledlayout(4,2);
ax1=nexttile(1, [4 1]);
title("Path of Total Cardiac Vector During Selected Cardiac Cycle(s)");
xlabel("M_{x} (Amp meters)");    ylabel("M_{Y} (Amp meters)");    zlabel("M_{Z} (Amp meters)");    view([45, 45]);
ax1.FontSize=16;    yL=0.2.*[min(M0(1,2,:)), max(M0(1,2,:))];   xL=0.2.*[min(M0(1,1,:)), max(M0(1,1,:))];   zL=0.2.*[min(M0(1,3,:)), max(M0(1,3,:))];
%Plot Axis originating from origin
txt_shift=-.1*10^-6;
line([0 0], yL , [0,0], 'LineWidth', 2, 'Color', 'k'); line(xL, [0 0], [0,0], 'LineWidth', 2, 'Color', 'k');  line([0,0], [0,0], zL, 'LineWidth', 2, 'Color', 'k');
text(xL(2)+txt_shift, txt_shift, txt_shift, "+X"); text(xL(1)+txt_shift, txt_shift, txt_shift, "-X"); text(txt_shift, yL(2)+txt_shift, txt_shift, "+Y"); text(txt_shift, yL(1)+txt_shift, txt_shift, "-Y"); text(txt_shift, txt_shift, zL(2)-txt_shift, "+Z");  text(txt_shift, txt_shift, zL(1)+txt_shift, "-Z");  
grid on;hold on;


ax3=nexttile(4, [2 1]);
title("Lead 2 ECG");
xlabel("Time (s)"); ylabel("Voltage (V)");
ax3.FontSize=12;
grid on; hold on;



cnt=1; Mc_plt0=[0,0,0]; end_cycle=end_ecg_idx;
for t=1:length(event_frames)
    if ~isempty(int_cycles)
        
        if first_ecg_idx+t<=end_cycle
            wave_marker=TCVicons(cnt);
        else
            if cnt==length(int_cycles)
                end_cycle=end_ecg_idx;
            else
                cnt=cnt+1;
                end_cycle=int_cycles(cnt);
            end
        end
        
    else
        wave_marker=TCVicons(1);

        end_cycle=end_ecg_idx;
    end

    [wave_color, ekg_wave_color, color]=EKG2Color(mappedWaves(event_frames(t)));
    TCV_marker=color+wave_marker;

    Mc_plt=sigma_b(event_frames(t)).*OrderSrcVectors(num_recon_src, M0(1,:,t)); 
%     Mc_plt=0.09*OrderSrcVectors(num_recon_src, M0(1,:,t));

    for k=1:num_recon_src
        plot3(ax1, [Mc_plt0(k,1), Mc_plt(k,1)], [Mc_plt0(k,2),Mc_plt(k,2)], [Mc_plt0(k,3),Mc_plt(k,3)], color+'-', 'LineWidth', 4); hold on;
    end

    if t+1<=length(event_frames)
        plot(ax3, (t:t+1)./Fs, L2_ecg(event_frames(t):event_frames(t+1)), ekg_wave_color);hold on;
    end
    Mc_plt0=Mc_plt;
end

grid on; hold on;


%% Plots all sets of cardiac informtion we obtain from ACT 5 system. 
% Once BSM computation time has been fixed, this will hopefully be able to
% be displayed and computed in real time (8/27/24)
if plt_all
    PlotData(Vq, ecg_wave, Mc, Mvc, [], fmdl, mappedWaves, s0, sN, model_info, plot_idx)
end
