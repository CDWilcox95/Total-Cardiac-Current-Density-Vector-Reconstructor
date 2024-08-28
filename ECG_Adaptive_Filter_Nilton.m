function [filter_ecg_wave, params]=ECG_Adaptive_Filter_Nilton(frame_ECG, start_frame, end_frame, el_1, el_2, first_iter, params)

% Multi-frequency FIR LMS filter to remove sinusoidal noise signals
% Nilton Rosa, July 2022.

%% Get frames
n_patterns = 32;

num_frames=end_frame-start_frame+1;



if first_iter

    [L,n_patterns, num_frames]=size(frame_ECG);

    s1 = squeeze(frame_ECG(:,el_1,(start_frame:end_frame)));
    d1= reshape(s1, [n_patterns*size(s1,2),1]);

    if el_2~=0
        s2=squeeze(frame_ECG(:,el_2, (start_frame:end_frame)));
        d2= reshape(s2, [n_patterns*size(s2,2),1]);

    end

    if el_2 ~=0
        signal = d2-d1; % signal
    else
        signal=d1;
    end


else
    d1=frame_ECG;

    signal=d1;

end

%% Inputs
%;
fn = [27.0017 54.0035 60.0248  81.0052 108.007 98.6554 120.08 135.008];  % Noise frequencies
Fs = 864.0553;            % sampling frequency
t2 = 1/Fs:1/Fs:(n_patterns*num_frames)/Fs;



signal = signal - mean(signal);
N = size(t2,2);
size_fn = size(fn,2);
s = signal;                 %ref input signal
%% Create Reference signal for each input frequency

C = params(1);                      % reference signal amplitude

ref_signal = zeros(size_fn,N);
for i=1:size_fn
    ref_signal(i,:) = C*cos(2*pi*fn(i)*t2);     % create all reference signals
end

%% LMS parameters
LMS_order =params(2);    % filter order
mu_LMS = params(3);       % step size
N_iterac_LMS =params(4);    % number of iterations for the LMS
params=[C, LMS_order, mu_LMS, N_iterac_LMS];

IIR_order = 40;      % Butterworth order
Butter_CoF = 150;    % Butterworth cuttoff frequency

%% Initialize variables
essemble_avg = zeros(N,size_fn);
filtered_signal = zeros(N,1);
error_array = zeros(N,size_fn);
W = zeros(LMS_order,size_fn);
Filt_samples = zeros(LMS_order,size_fn);
y = zeros(1,N);


essemble_avg = zeros(N,size_fn);
error_array = zeros(N,size_fn);
W_LMS = zeros(LMS_order,size_fn);
u_LMS = zeros(LMS_order,size_fn);

ecg_avg = 0;
ecg_avg_samples = 1;
LMS_out = NaN(size(signal));
%% Run LMS
t_LMS  = 0;
for i=1:N
    tic
    aux_signal = signal(i);
    for k = 1: size_fn

        u_LMS(1,k) = C*cos(2*pi*fn(k)*t_LMS);
        
        xn_LMS = aux_signal;   % get the current signal sample
        e_e = 0;
        for j=1:N_iterac_LMS    % run over the iterations

            yn_LMS = (W_LMS(:,k)'*u_LMS(:,k));  % compute the filter output
            e_LMS = xn_LMS-yn_LMS;                      % compute the error
            e_e = e_e + (e_LMS)^2;
            
            % Compute the next weight
            W_LMS(:,k) = (W_LMS(:,k)+(mu_LMS*e_LMS*u_LMS(:,k)));   

        end
        u_LMS(2:end,k) = u_LMS(1:end-1,k);
        error_array(i,k) = error_array(i,k)+e_LMS;  % assign the error
        aux_signal = e_LMS;                         % filter output
        
        % Essemble averaging method - Learning curve
        essemble_avg(i,k) = essemble_avg(i,k)+e_e/N_iterac_LMS;   
    end
    
    t_LMS = t_LMS + (1*1)/Fs; % advance time
    
    LMS_out(i) = aux_signal;  % save LMS output
    
end
y=LMS_out;
%% Butterworth filter

Fn=Fs/2;        % Nyqusit freq
% implement 2nd order LPF
wo=Butter_CoF/Fn; % filter BW 40Hz/150Hz 
[b_L,a_L]=butter(IIR_order,wo); % second order butterworth LPF

% Filtfilt does not add phase shift
y=filtfilt(b_L,a_L,y); % LPF filter 

disp("Done!");








%% Low-Pass IIR filter

% Fn=Fs/2;        % Nyqusit freq
% % implement 2nd order LPF
% wo=150/Fn; % filter BW 40Hz/150Hz (passive EIT works with both, Active works only with 40 Hz)
% [b_L,a_L]=butter(20,wo); % second order butterworth LPF. May need to increase the order depending on the leads
% 
% y=filtfilt(b_L,a_L,y); % LPF filter
%% Result of adaptive noise cancellation
x0=10;
y0=80;
width=850;
height=550;
% figure(1)
% set(gcf,'position',[x0,y0,width,height])
% subplot(3,1,1)
% plot(s(end-5000:end-1000))
% xlim([1 4000]);
% title("Raw ECG Signal")
% 
% subplot(3,1,2)
% plot(y(end-5000:end-1000))
% xlim([1 4000]);
% title("Adaptive LMS Filter")
% 
% subplot(3,1,3)
% plot(s(end-5000:end-1000))
% hold on
% plot(y(end-5000:end-1000))
% hold off
% xlim([1 4000]);
% title("Raw and Filtered signal comparison")
% legend("Raw","Filtered")

%% FFT
% x1 = 5000;       	%   steady state
% F = y(x1:end);
% Y = fft(F);
% % L = ceil(size(F,2));
% L = ceil(size(F,1));
% 
% P2 = abs(Y/(L));
% P1 = P2(1:round((L)/2+1));
% P1(2:end-1) = 2*P1(2:end-1);
% % f = linspace(0,Fs/2,size(P1,2));
% f = linspace(0,Fs/2,size(P1,1));

% figure(2)
% title('FFT')
% xlabel('f (Hz)')
% plot(f,P1)
% ylabel('|P1(f)|')
%% Plot error performance

% f_idx = 4;          % error for a given frequency
% 
% figure(3)
% plot(essemble_avg(:,end))
% ylabel('ensemble-average squared error')
% xlabel('Sample')
% title('LMS - Convergence rate ')

y = (y/(2^18))*1000;
y = y - sum(y)/length(y);

filter_ecg_wave=y;