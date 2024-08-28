function TCVPlotter(sigma0, Q0, M0, mappedWaves, event_frames, ecg, video_plt)

[K, dim]=size(Q0);
[~,~,num_frames]=size(M0);
x_bnds=Q0(1)+mean(sigma0).*[min(M0(1,1,:)), max(M0(1,1,:))];    y_bnds=Q0(2)+mean(sigma0).*[min(M0(1,2,:)), max(M0(1,2,:))];   z_bnds=Q0(3)+mean(sigma0).*[min(M0(1,3,:)), max(M0(1,3,:))];
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(1,2);

Fs=864.0553;

axecg=nexttile(2);
xlabel("Time (s)"); ylabel("Voltage (mV)"); grid on; hold on;
for t=1:num_frames-1
    [wave_color, ekg_wave_color, color]=EKG2Color(mappedWaves(event_frames(t)));
    plot(axecg, [t:t+1]./Fs, (10^3).*ecg(event_frames(t:t+1)), ekg_wave_color);
end
hold on;
vidfile=VideoWriter("Multi_Src_Reconstructions.mp4", 'MPEG-4');
vidfile.FrameRate=10;

open(vidfile);


for t=1:num_frames
    [wave_color, ekg_wave_color, color]=EKG2Color(mappedWaves(event_frames(t)));
    
    vec_plt_t=[];
    for k=1:K
        mk=M0(1,:,:);
        nexttile(k);
        title("Dipole Source "+num2str(k));
        xlabel(" m_{k_{x}} (Amp x Meters)");    ylabel(" m_{k_{y}}  (Amp x Meters)");      zlabel(" m_{k_{z}}  (Amp x Meters)");   view([45, 45]);
        xlim(x_bnds);   ylim(y_bnds);   zlim(z_bnds);
%         view([60 35]);

        view([65 35]);
        hold on;
        
        grid on;
        
        plot3(Q0(1),Q0(2),Q0(3), 'g+', 'MarkerSize', 10);
        plot3(Q0(1)+sigma0(t).*mk(1,1,t), Q0(2)+sigma0(t).*mk(1,2,t), Q0(3)+sigma0(t).*mk(1,3,t), wave_color, 'MarkerSize', 5, 'MarkerFaceColor', color);
        qk=quiver3(Q0(1),Q0(2),Q0(3), sigma0(t)*mk(1,1,t), sigma0(t)*mk(1,2,t), sigma0(t)*mk(1,3,t), 'c', 'LineWidth', 2);
        vec_plt_t(k)=qk;
        hold on;    

        
    end
    
    marker_plt_ecg=plot(axecg, t/Fs, (10^3).*ecg(event_frames(t)), 'bo', 'MarkerSize', 6);
    
    
    frame = getframe(gcf);
    writeVideo(vidfile,frame);
    
    delete(marker_plt_ecg);
    delete(vec_plt_t);
    
end




close(vidfile);
end