function U0=Computeu0OnMesh(fmdl, sigma0, Q, m)

Pt=fmdl.nodes;N=length(Pt);
[K,~]=size(Q);

U0=zeros(N,1);


u0=@(p, Qk) (1/(4*pi*sigma0))*dot(p-Qk, m)/norm(p-Qk,2)^3;


for k=1:K
    Qk=Q(k,:);
    for n=1:N
        pt=Pt(n,:);

        U0(n)=U0(n)+u0(pt, Qk);
    end

end

end