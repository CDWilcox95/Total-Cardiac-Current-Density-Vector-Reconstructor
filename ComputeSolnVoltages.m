function V=ComputeSolnVoltages(el_area, fmdl, Uh, U0)


    T=fmdl.elems;   Pt=fmdl.nodes;
    M=length(T);    N=length(Pt);

    bndry_idx=fmdl.boundary;
    L=length(fmdl.electrode);

    % ucoef=ComputeSolnCoef(fmdl, Uh);


    bndry_elec_map=fmdl.bndry2elecmap;




    V=zeros(L,1);
    for l=1:L
        el_bndryidx=bndry_idx(bndry_elec_map(:,1)==l,:);

        el_elemA=bndry_elec_map(bndry_elec_map(:,1)==l, 2);

        uc_int=0;   el_area=0;  
        for t=1:length(el_bndryidx)
            if ~isempty(Uh)
                uc_int=uc_int+(el_elemA(t)/3)*( sum(Uh(el_bndryidx(t,:))+U0(el_bndryidx(t,:))) ) ;
                el_area=el_area+el_elemA(t);
            else
                uc_int=uc_int+(el_elemA(t)/3)*(sum(U0(el_bndryidx(t,:))));
                el_area=el_area+el_elemA(t);
            end
        end

        V(l)=(1/el_area)*uc_int;
    end

    


end