function [Q0_new, multi_Q0, frame_err, frame_idx, mesh_errs_t]=FindBestQ0Sphere(sigmab, Bh, heart_pts, elec_pts, model_info, Vd, method, fwd_map, construct_map, num_src, file_ext)

R=model_info.R; H=model_info.H; L=model_info.num_elec;  fmdl=model_info.FEM_Mesh;   
electrode_nodes=fmdl.electrode_idx;
% load("SourcePts_K"+num2str(num_src)+".mat");

el_area=pi*0.0254^2;
f0=768;
switch method

    case "FitFwdSims"

        e1=[1,0,0]'; e2=[0,1,0]'; e3=[0,0,1]';  E=[e1, e2, e3];         %Elementary Basis

        

        T=fmdl.elems;                           % Mesh Elements       (T)
        Pt=fmdl.nodes;                          % Mesh Points         (Pt)
        fmdl.nodes=Pt;  fmdl.electrode_idx=electrode_nodes;

        [C, vol_tetr]=CalculateMeshCoef(T, Pt);
        N=length(Pt);
        mesh_data.elem_labels=fmdl.mat_idx;
        mesh_data.Pts=Pt;   mesh_data.radius=R; mesh_data.height=H; mesh_data.total_elems=T;


        fmdl=GetReconMesh(R,H, false, "const", 'c', []);  Pt=fmdl.nodes;

        %Initial Guess
        theta_l13=(2*pi/16)*13; x13=-(model_info.R/2)*cos(theta_l13);   y13=-(model_info.R/2)*sin(theta_l13); z13=model_info.H/2;
        x0=[x13, y13, z13];

        if isempty(heart_pts)
            [heart_pts, ~]=GenerateSourceOnSphere(Bh(4), Bh(1:3), 250);
        end
            
        numHpts=length(heart_pts);



        temp_pts=heart_pts;
        for i=1:numHpts
            if norm(heart_pts(i,1:2),2)>=model_info.R
                temp_pts(i,:)=-1;
            end
        end
        heart_pts=temp_pts;
        numHpts=length(heart_pts);
        [~,num_data_volt]=size(Vd);

        switch fwd_map
            case "FEM"

                if construct_map

                    fmdl=GetReconMesh(R,H, false, "const", 'c', []);    


                    Pt=fmdl.nodes; T=fmdl.elems;    C=fmdl.coef;    vol_tetr=fmdl.volumes;  N=length(Pt);
                    electrode_nodes=fmdl.electrode_idx; elec_pts=Pt(electrode_nodes,:); 

                    
                    fprintf("Generating Forward Maps Using FEM Construction... \n");
                    
                    G=zeros(L,3, numHpts);

                    sigma_dist=1;
                    if length(sigma_dist)==1

                        rhs_setting="const";
                        sigma_Omega=ones(length(T),1);
%                         sigma_Omega=1;
                        sigma0=1;
                    elseif length(sigma_dist)==4
                        rhs_setting="var";
                        sigma_Omega=SetCondBodySims(R,H, fmdl, sigma_dist(1), sigma_dist(2), sigma_dist(3), sigma_dist(4));

                        % true_sigma=sigmaOmega_t;
                        elem_ind=1:length(T);
                        fmdl.heart_elems=elem_ind(sigma_Omega==sigma_dist(2));
                        sigma0=sigma_dist(2);
                    else
                        rhs_setting="var";
                        sigma_Omega=sigma_dist;
                        Bh=model_info.seg_info.H;
                        Helems=GetHeartElemsFromSpecPt(Bh(1:3), Bh(4), fmdl);

                        fmdl.heart_elems=Helems;
                        sigma0=mean(sigma_Omega(Helems));
                        sigma_Omega(Helems)=sigma0;
                    end


                    mesh_data.elem_labels=fmdl.mat_idx;
                    mesh_data.Pts=Pt;   mesh_data.radius=R; mesh_data.height=R; mesh_data.total_elems=T;
                    switch model_info.shape
                        case 'c'
%                             S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("cyl_FWD_MAPFEMStiffMatrix.mat", 'S');
                            load("cyl_FWD_MAPFEMStiffMatrix.mat");
            
                        case 'e'
                            S=Create3DStiffMatrix(C, vol_tetr, sigma_Omega, T, Pt, "FWD_MAP");  save("ell_FWD_MAPFEMStiffMatrix.mat", 'S');
                            load("ell_FWD_MAPFEMStiffMatrix.mat");
            
                    end
                    g=CreateUniqueVector("pt_elec", electrode_nodes, fmdl, Pt, T);
                    % cond_setting="Constant";    rhs_setting="const";file_ext="BestQ0";
                    % S0=Create3DStiffMatrix(C, vol_tetr, model_info.ref_cond, T, Pt,"BestQ0");    save("FEMStiffMatrix_1cond_O2.mat", 'S0');
                    % load("FEMStiffMatrix_1cond_O2.mat");          
                    S0=S;
                    S0(1:N, N+1)=g;    S0(N+1,1:N)=g';


                    elec_nodes=fmdl.electrode_idx;
                    line_length=fprintf("Number of Forward Maps Computed:  %d/%d", 0, numHpts);

                    for k=1:numHpts
                        Qk=heart_pts(k,:);

                        G0=FWD_EKG_infmap(elec_pts, Qk);

                        for i=1:3
                            b=CalculateRHSvec(rhs_setting, T, Pt, Qk, sigma_dist, C, vol_tetr, fmdl, E(:,i)');  b(N+1)=0;

                            U_ki=S0\b;      V_ki=U_ki(electrode_nodes);

                            u0_mesh=Computeu0OnMesh(fmdl, 1, Qk, E(:,i)');


                            U0=ComputeSolnVoltages(el_area, fmdl, U_ki, u0_mesh);

                            Gi=U0;

                            % Gi=(G0(:,i)+V_ki);
                            G(:,i,k)=Gi-sum(Gi)/L;
                        end

                        fprintf(repmat('\b',1,line_length));
                        line_length=fprintf("Number of Forward Maps Computed:  %d/%d", k, numHpts);
                    end


                    fprintf("\n Forward Map Generation Completed! \n");
                    test_pt_gen="";
                    save(file_ext+test_pt_gen+"BestFit_procG.mat", 'G');
                else
                    test_pt_gen="";
                    load(file_ext+test_pt_gen+"BestFit_procG.mat");
%                     load(file_ext+"BestFit_procGvarcond.mat");
                end
            case "inf"

                
                if construct_map
                    elec_nodes=fmdl.electrode_idx;
                    cnt=1;


                    for k=1:numHpts
                        Q=heart_pts(k,:);
                        G(:,:,cnt)=FWD_EKG_infmap(elec_pts, Q);
                        cnt=cnt+1;
                    end

                    save("BestFit_procG0.mat", 'G');
                else
                    load("BestFit_procG0.mat");
                end
        end

        fprintf("Simulating Source Points for "+num2str(num_src)+" Source Problem \n");
        line_length=fprintf("Percentage of Points Simulated:  %.4g %%",100*(0/numHpts));
        if num_src>1

            cnt=1;  src_idx=zeros(num_src,1);
            for k=1:num_src

                for t=1:num_data_volt
                    for j=1:numHpts
                        if k>1 && j==src_idx(k-1,t)
                            break;
                        end
                        Gj=G(:,:,j);
                        Gtest=zeros(L,3*k);

                        for i=1:k-1
                            Gtest(:,3*i-2:3*i)=G(:,:,src_idx(i));
                        end
                        Gtest(:,3*k-2:3*k)=Gj;
                        xt=inv(Gtest'*Gtest)*Gtest'*Vd(:,t);


                        V_sim(:,t)=Gtest*xt;
                        l2_err(j,t)=norm(V_sim(:,t)-Vd(:,t),2);


                        fprintf(repmat('\b',1,line_length));
                        line_length=fprintf("Percentage of Frames Simulated:  %.4g %%",100*(j/numHpts));
                    end
                    [val, idx]=min(sum(l2_err, 2));
                    src_idx(k,t)=idx;
                end
                
    
            end
            err_val=val;
            Q0=heart_pts(src_idx(:,1),:);

        else
            cycle_err=zeros(numHpts,1);
            l2_err=zeros(numHpts, num_data_volt);
            f0=768-334;
%             [Uref_sim, Vref]=SimulateRefFrameECG(Vd, G, f0, sigmab)
            for j=1:numHpts
                Gtest=G(:,:,j);

                V_sim=zeros(L, num_data_volt);
                for t=1:num_data_volt
                    

                    [xt, ~]=RegSystem(Gtest, Vd(:,t), 'none', []);
                    V_sim(:,t)=Gtest*xt;
                    
                    l2_err(j,t)=norm(V_sim(:,t)-Vd(:,t),2);
                end
                cycle_err(j)=CC_RelativeError(V_sim, Vd);
                

                fprintf(repmat('\b',1,line_length));
                line_length=fprintf("Percentage of Points Simulated:  %.4g %%",100*(j/numHpts));

            end

            [err_val, minCC_idx]=min(cycle_err);


            Q0=heart_pts(minCC_idx,:);
        end




        frame_err=zeros(num_data_volt,1); frame_idx=zeros(num_data_volt,1);
        for t=1:num_data_volt
            [err_t, mint_idx]=min(l2_err(:,t));
            frame_err(t)=err_t; frame_idx(t)=mint_idx;
            mesh_errs_t(:,t)=l2_err(:,t);
        end
        save("source_errs_t.mat", 'mesh_errs_t');
        
        multi_Q0=heart_pts(frame_idx,:);    [max_err, max_idx]=max(frame_err);

        Q0_new=Q0;
        save("BestQ0w_"+fwd_map+"_K"+num2str(num_src)+".mat", 'Q0');

        fprintf("\n Best Fit Q0 Found! \n");
        fprintf("Cycle Best Fit Specs: \n ___________________________________________\n");
        fprintf("Distance from Bh center to Q0: ||Bh-Q0||=%.3g \n", norm(Bh(1:3)-Q0,2));
        fprintf("l2 Error of Q0: %.4g \n\n", err_val);

        fprintf("Moving Dipole Best Fit Specs: \n ___________________________________________\n");
        fprintf("Maximum Error: %.4g at frame number %d \n\n", max_err,max_idx);
end


