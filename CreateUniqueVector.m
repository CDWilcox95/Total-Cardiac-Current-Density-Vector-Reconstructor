function g=CreateUniqueVector(unique_cond, elec_nodes, fmdl, Pt, T)
    N=length(Pt); 
    M=length(T);

    g=zeros(N,1);
    switch unique_cond

        case "pt_elec"
            g(elec_nodes)=1;

        case "AveGap"



        case "integral"
            boundary_nodes=fmdl.boundary;
            boundary_labels=fmdl.boundary_numbers;

            T_nodes=boundary_nodes(boundary_labels==3,:);   num_T_elems=length(T_nodes);
            for i=1:num_T_elems
                surf_tri=T_nodes(i,:);
                tri_pts=Pt(surf_tri,:);

                v1=tri_pts(2,:)'-tri_pts(1,:)';
                v2=tri_pts(3,:)'-tri_pts(1,:)';

                A=norm(cross(v1,v2), 2)/2;
                g(surf_tri)=g(surf_tri)+A/3;
            end

            S0_nodes=boundary_nodes(boundary_labels==1,:);   num_S0_elems=length(S0_nodes);
            for i=1:num_S0_elems
                surf_tri=S0_nodes(i,:);
                tri_pts=Pt(surf_tri,:);

                v1=tri_pts(2,:)'-tri_pts(1,:)';
                v2=tri_pts(3,:)'-tri_pts(1,:)';

                A=norm(cross(v1,v2), 2)/2;
                g(surf_tri)=g(surf_tri)+A/3;
            end


            B_nodes=boundary_nodes(boundary_labels==2,:);   num_B_elems=length(B_nodes);
            for i=1:num_B_elems
                surf_tri=B_nodes(i,:);
                tri_pts=Pt(surf_tri,:);

                v1=tri_pts(2,:)'-tri_pts(1,:)';
                v2=tri_pts(3,:)'-tri_pts(1,:)';

                A=norm(cross(v1,v2), 2)/2;
                g(surf_tri)=g(surf_tri)+A/3;
            end

    end

    save("g_vec.mat", 'g');

end