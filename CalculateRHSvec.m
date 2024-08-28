function b=CalculateRHSvec(cond_setting, T, Pt, Q, sigma_dist, C, vol_tetr, fmdl, M)
if length(sigma_dist)==length(T)
    sigma0=mean(sigma_dist(fmdl.heart_elems));
    sigma_dist(fmdl.heart_elems)=sigma0;
elseif length(sigma_dist)>1
    sigma0=sigma_dist(2);
else
    sigma0=sigma_dist;
end
[num_elems, elem_nodes]=size(T);
if elem_nodes==4
    order=1;
else
    order=2;
%     load("bndry_nodes_to_tetrahedra_elem.mat");
    map=fmdl.face2elem;
end

N=length(Pt);
[K,~]=size(Q);
Heart_sigma=sigma0;

boundary_nodes=fmdl.boundary;
[T_lbl, B_lbl, S0_lbl]=DetermineBoundaryLabels(fmdl, Pt, 0.15, max(Pt(:,3)));

if order==2
    map_T=map(T_lbl,:);
    map_S0=map(S0_lbl,:);
    map_B=map(B_lbl,:);
end


%% Calculate S0 Vector (source with corrective term conductivity)
S0=zeros(N,K);
if strcmp(cond_setting, "var")
    if length(sigma_dist)==4
        R=0.15; H=0.15;
        sigma_Omega=SetCondBodySims(R,H, fmdl, sigma_dist(1), sigma_dist(2), sigma_dist(3), sigma_dist(4), true);
    else
        sigma_Omega=sigma_dist;
    end
        for k=1:K
            S0(:,k)=CalculateS0Vector(Pt, T, Q(k,:), fmdl, C, vol_tetr, sigma_Omega, M(k,:), sigma0);
        end
end

%         load("InfMap_StiffMat.mat");

%% Calculate Neumann Boundary Condition Vector
du0_dx=@(p, Q0, m)  (1/(4*pi*sigma0))*(m(1)*norm(p-Q0,2)^3-dot(p-Q0,m)*(3*(p(1)-Q0(1))*norm(p-Q0,2)))/norm(p-Q0,2)^6;
du0_dy=@(p, Q0,m)  (1/(4*pi*sigma0))*(m(2)*norm(p-Q0,2)^3-dot(p-Q0,m)*(3*(p(2)-Q0(2))*norm(p-Q0,2)))/norm(p-Q0,2)^6;
du0_dz=@(p, Q0,m)  (1/(4*pi*sigma0))*(m(3)*norm(p-Q0,2)^3-dot(p-Q0,m)*(3*(p(3)-Q0(3))*norm(p-Q0,2)))/norm(p-Q0,2)^6;
grad_u0=@(p, Q0,m) [du0_dx(p,Q0,m), du0_dy(p,Q0,m), du0_dz(p,Q0,m)];





n_circle=@(p) [p(1),p(2),0]./norm([p(1),p(2)],2);

b=zeros(N,1);
for v=1:K
    Qk=Q(v,:);  mk=M(v,:);  bk=zeros(N,1);
    b1=zeros(N,1);   b2=zeros(N,1);  b3=zeros(N,1);
    T_nodes=boundary_nodes(T_lbl,:);   num_T_elems=length(T_nodes);


    for j=1:num_T_elems
        surf_tri=T_nodes(j,:);
        tri_pts=Pt(surf_tri,:);
        brcntr=[sum(tri_pts(1:3,1)), sum(tri_pts(1:3,2)), sum(tri_pts(1:3,3))]./3;

        v1=tri_pts(2,:)-tri_pts(1,:);
        v2=tri_pts(3,:)-tri_pts(1,:);
        tri_normal=cross(v1,v2)/norm(cross(v1,v2),2);

        A=norm(cross(v1,v2), 2)/2;
        h=0;
        switch order
            case 1
                for i=1:3
                    pt_i=tri_pts(i,:);
                    b1(surf_tri(i))=b1(surf_tri(i))+(Heart_sigma*A/3)*du0_dz(pt_i, Qk,mk);

                    % b1(surf_tri(i))=b1(surf_tri(i))+(Heart_sigma*A/3)*dot(grad_u0(pt_i, Qk,mk), tri_normal);

                end

            case 2
                tetr_elem=map_T(j,1);   tetr_idx=map_T(j,2:end);
                coef_t=C(tetr_elem,:,:);    coef_t=reshape(coef_t, [elem_nodes,elem_nodes]);
                

                for k=1:6
                    c_k=coef_t(:,tetr_idx(k));
                    int_sum=0;

                    for i=1:6
                        pt_i=tri_pts(i,:);
                        if i<=3
                            int_sum=int_sum+(A/20)*du0_dz(pt_i, Qk, mk)*phi2(c_k, pt_i);
                        else
                            int_sum=int_sum+(2*A/15)*du0_dz(pt_i, Qk, mk)*phi2(c_k, pt_i);
                        end

                    end
                    int_sum=int_sum+(9*A/20)*du0_dz(brcntr, Qk, mk)*phi2(c_k, brcntr);
                    b1(surf_tri(k))=b1(surf_tri(k))+Heart_sigma*int_sum;

                end


        end
    end


    % Calculate Side Surface Vector
    %         S0_nodes=boundary_nodes(boundary_labels==S0_lbl,:);      num_S0_elems=length(S0_nodes);
    S0_nodes=boundary_nodes(S0_lbl,:);      num_S0_elems=length(S0_nodes);

    for j=1:num_S0_elems
        surf_tri=S0_nodes(j,:);
        tri_pts=Pt(surf_tri,:);
        brcntr=[sum(tri_pts(1:3,1)), sum(tri_pts(1:3,2)), sum(tri_pts(1:3,3))]./3;
        
        if norm(brcntr(1:2),2)<0.1495
            continue;
        end

        v1=tri_pts(2,:)-tri_pts(1,:);
        v2=tri_pts(3,:)-tri_pts(1,:);

        tri_normal=cross(v1,v2)/norm(cross(v1,v2),2);
        A=norm(cross(v1,v2), 2)/2;
        h=0;
        switch order
            case 1
                for i=1:3
                    pt_i=tri_pts(i,:);
                    b2(surf_tri(i))=b2(surf_tri(i))+(Heart_sigma*A/3)*dot(grad_u0(pt_i, Qk,mk), tri_normal);
                end

            case 2
                tetr_elem=map_S0(j,1);   tetr_idx=map_S0(j,2:end);
                coef_t=C(tetr_elem,:,:);    coef_t=reshape(coef_t, [elem_nodes,elem_nodes]);
                
                for k=1:6
                    int_sum=0;
                    c_k= coef_t(:,tetr_idx(k));

                    for i=1:6
                        pt_i=tri_pts(i,:);
                        if i<=3
                            int_sum=int_sum+(A/20)*dot(grad_u0(pt_i, Qk,mk), tri_normal)*phi2(c_k, pt_i);
                        else
                            int_sum=int_sum+(2*A/15)*dot(grad_u0(pt_i, Qk,mk), tri_normal)*phi2(c_k, pt_i);
                        end

                    end
                    int_sum=int_sum+(9*A/20)*dot(grad_u0(brcntr, Qk,mk), tri_normal)*phi2(c_k, brcntr);

                    b2(surf_tri(k))=b2(surf_tri(k))+Heart_sigma*int_sum;
                end

        end

    end


    % Calculate Bottom Surface Vector
    %         B_nodes=boundary_nodes(boundary_labels==B_lbl,:);   num_B_elems=length(B_nodes);
    B_nodes=boundary_nodes(B_lbl,:);   num_B_elems=length(B_nodes);
    for j=1:num_B_elems
        surf_tri=B_nodes(j,:);
        tri_pts=Pt(surf_tri,:);
        brcntr=[sum(tri_pts(1:3,1)), sum(tri_pts(1:3,2)), sum(tri_pts(1:3,3))]./3;

        v1=tri_pts(2,:)-tri_pts(1,:);
        v2=tri_pts(3,:)-tri_pts(1,:);

        tri_normal=cross(v1,v2)/norm(cross(v1,v2),2);
        A=norm(cross(v1,v2), 2)/2;
        
        h=0;
        switch order
            case 1
                for i=1:3
                    pt_i=tri_pts(i,:);
                    % ans1=-dot(grad_u0(pt_i, Qk,mk), tri_normal);
                    % ans2=du0_dz(pt_i, Qk,mk);

                    b3(surf_tri(i))=b3(surf_tri(i))+(Heart_sigma*A/3)*dot(grad_u0(pt_i, Qk,mk), tri_normal);

                    % b3(surf_tri(i))=b3(surf_tri(i))+(Heart_sigma*A/3)*du0_dz(pt_i, Qk,mk);
                end

            case 2
                tetr_elem=map_B(j,1);   tetr_idx=map_B(j,2:end);
                coef_t=C(tetr_elem,:,:);    coef_t=reshape(coef_t, [elem_nodes,elem_nodes]);
                
                for k=1:6
                    int_sum=0;
                    c_k=coef_t(:,tetr_idx(k));

                    for i=1:6
                        pt_i=tri_pts(i,:);
                        if i<=3
                            int_sum=int_sum+(A/20)*du0_dz(pt_i, Qk, mk)*phi2(c_k, pt_i);
                        else
                            int_sum=int_sum+(2*A/15)*du0_dz(pt_i, Qk, mk)*phi2(c_k, pt_i);
                        end

                    end
                    int_sum=int_sum+(9*A/20)*du0_dz(brcntr, Qk, mk)*phi2(c_k, brcntr);

                    b3(surf_tri(k))=b3(surf_tri(k))-Heart_sigma*int_sum;
                end

        end

    end

    bk=-(b1+b2+b3);

    b=b+bk-S0(:,v);
end



end