function bndry_elec_map=GetElectrodeElems(fmdl)
    
    T=fmdl.elems;   Pt=fmdl.nodes;
    M=length(T);    N=length(Pt);
    bndry_idx=fmdl.boundary;    M0=length(bndry_idx);

    bndry2elem=zeros(M0,2);
    L=32;

    bndry_elec_map=zeros(M0,1);
    for j=1:M0
        b_idx=bndry_idx(j,:);
        tri_pts=Pt(b_idx,:);

        v1=tri_pts(2,:)-tri_pts(1,:);
        v2=tri_pts(3,:)-tri_pts(1,:);

        A=norm(cross(v1,v2), 2)/2;

        
        for l=1:L
            el=fmdl.electrode(l).nodes;
            if length(intersect(b_idx, el))==3
                bndry_elec_map(j,1)=l;
                bndry_elec_map(j,2)=A;
                break;
            end
        end

    end


    
end