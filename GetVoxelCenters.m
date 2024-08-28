function mid_pts=GetVoxelCenters(Jt_mesh)
    
    N=length(Jt_mesh);
    mid_pts=zeros(N,3);
    for n=1:N
       vox=Jt_mesh(n,:);
       ra=vox(1);   rb=vox(2);  ta=vox(3);  tb=vox(4);  za=vox(5);   zb=vox(6);
       r0=0.5*(ra+rb);   t0=0.5*(ta+tb);    z0=0.5*(za+zb); 
       
       if t0>2*pi
           t0=mod(t0,2*pi);
       end
        
       x0=r0*cos(t0);    y0=r0*sin(t0); 
       mid_pts(n,:)=[x0,y0,z0];
    end


end