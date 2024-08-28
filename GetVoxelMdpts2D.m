function mid_pts=GetVoxelMdpts2D(Jt_mesh)

    N=length(Jt_mesh);
    mid_pts=zeros(N,2);
    for n=1:N
       vox=Jt_mesh(n,:);
       ra=vox(1);   rb=vox(2);  ta=vox(3);  tb=vox(4);
       r0=0.5*(ra+rb);   t0=0.5*(ta+tb);  
       
       if t0>2*pi
           t0=mod(t0,2*pi);
       end
        
       x0=r0*cos(t0);    y0=r0*sin(t0); 
       mid_pts(n,:)=[x0,y0];
    end
end