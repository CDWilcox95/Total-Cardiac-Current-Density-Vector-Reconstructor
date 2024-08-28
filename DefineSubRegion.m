function sub_vox=DefineSubRegion(JTmesh)

    [N,~]=size(JTmesh);
    sub_vox=zeros(N,1);
    
    
    Ta=0;       Tb=3*pi/4;
    
    Ta=5*pi/4;   Tb=2*pi;   %DICOM
    tau=0.8;
    Ra=0;   Rb=tau*max(JTmesh(:,2));
    for n=1:N
        vox=JTmesh(n,:);
        ra=vox(1);  rb=vox(2);  ta=vox(3);  tb=vox(4);
        
        if ta >= Ta && tb <= Tb && rb <= Rb && ra >= Ra
            sub_vox(n)=1;
        end
        
        
    end

    sub_vox=logical(sub_vox);

end