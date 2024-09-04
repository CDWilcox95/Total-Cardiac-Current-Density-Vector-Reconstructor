function Vq=InterpolateVekgBSM(model_info, Vekg)

[L, num_frames]=size(Vekg);
if compute_new
    Z=model_info.elec_planes;
    top_elecs=L/2+1:L;  bot_elecs=1:L/2;
    Nz=51;                              % number of faces along z-axis = Nz-1
    Nt=1000; %ceil((2*pi*R)/((Z(2)-Z(1))*Nz)); % number of faces along perimeter
    z=linspace(Z(1),Z(2),51);
    t=linspace(0,2*pi,Nt+1); t(end)=0;
    theta_pts=[(2*pi/(L/2)).*[0:L/2], (2*pi/(L/2)).*[0:L/2]];
    z_pts=[Z(1).*ones(L/2+1,1); Z(2).*ones(L/2+1,1)];

    [Tq,Zq]=meshgrid(t,z);
    % for s=1:end_ecg_idx-first_ecg_idx+1
    for s=1:num_frames
        V_idx=Vekg(:,s);  V_bot=V_idx(bot_elecs); V_top=V_idx(top_elecs);
        V_bot(end+1)=V_bot(1);      V_top(end+1)=V_top(1);
        V_idx=[V_bot; V_top];

        Vq(:,:,s)=griddata(theta_pts, z_pts, V_idx, Tq, Zq, "cubic");
    end

end

end