function [T_lbl, B_lbl, S0_lbl, intsurf_lbl]=DetermineBoundaryLabels(fmdl, Pt, R, H)
    boundary_nodes=fmdl.boundary;

    num_boundary_elems=length(boundary_nodes);
    T_lbl=zeros(num_boundary_elems,1);  B_lbl=T_lbl;    S0_lbl=T_lbl;   intsurf_lbl=T_lbl;
    for t=1:num_boundary_elems
        surf_tri=boundary_nodes(t,:);
        tri_pts=Pt(surf_tri,:);

        if tri_pts(:,3)==0
            B_lbl(t)=t;
        elseif tri_pts(:,3)==H
            T_lbl(t)=t;
        else 
%             cnt=0;
            S0_lbl(t)=t;
%             for i=1:3
%                 if norm(tri_pts(i,1:2))>=R-0.1
%                     cnt=cnt+1;
%                 end
%             end
%             if cnt==3
%                 S0_lbl(t)=t;
%             else
%                 intsurf_lbl(t)=t;
%             end
        end

    end

    B_lbl=nonzeros(B_lbl);  T_lbl=nonzeros(T_lbl);  S0_lbl=nonzeros(S0_lbl);    intsurf_lbl=nonzeros(intsurf_lbl);


end