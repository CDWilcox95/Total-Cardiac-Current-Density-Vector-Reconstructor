%% SNR Analysis
close all;
clear all;
format long;

%%  Input parameters
L=32;                   % number of electrodes on the surface of the body
K=L-1;                  % number of current pattern

fixed_range = 1;
video = 1;              % 1 - Make a video, 0 - don't
start_frame =1;         % start frame
end_frame = 1000;         % end frame

%% Load dataset
path = "C:\Users\nilto\OneDrive - Colostate\backup\Matlab\ACT5_data\Subject007_22_02_10";
file_name = "Subj007_140kHz__22_02_10_10_55_24";

path = path + "\" + file_name;
load(path);

%% Load Direct Problem Solution
load forward_gap3D_act5_subj007
U = U0;

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
    gammabest=(1/rhobest)*1000;                   %convert to mS/M
    sigma_b(index)= gammabest;
    
    index = index + 1;
end

[max_val, Ref_frame] = max(real(sigma_b));

plot(real(sigma_b))