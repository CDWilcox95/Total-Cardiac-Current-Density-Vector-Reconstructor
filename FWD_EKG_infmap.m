function A=FWD_EKG_infmap(Xl, Q)

[numSrc,~]=size(Q);
[L,dim]=size(Xl);
A=zeros(L,dim*numSrc);


theta_l=[(2*pi/(L/2)).*[1:L/2], (2*pi/(L/2)).*[1:L/2]];
f=0.5;  dist_gap=theta_l(2)-theta_l(1);
h=0.0254;   R=max(vecnorm(Xl(:,1:2),2,2));  

area_el=pi*h^2;

u0p=@(Qp, M, r, theta, z) (1/(4*pi))*(M(1).*(r.*cos(theta)-Qp(1).*cos(Qp(2)))+M(2).*(r.*sin(theta)-Qp(1).*sin(Qp(2)))+ M(3).*(z-Qp(3)))./(r.^2 + Qp(1).^2-2.*r.*Qp(1).*cos(theta-Qp(2))+(z-Qp(3)).^2).^(3/2);


for k=1:numSrc
    Qk=Q(k,:);
    Qp=[norm(Qk(1:2),2), atan2(Qk(2), Qk(1)), Qk(3)];

    for l=1:L

        for n=1:dim
            src_vec=zeros(3,1); src_vec(n)=1;
            Vkl=(1/area_el)*integral2(@(theta,z) u0p(Qp, src_vec,R, theta,z), theta_l(l)-f*dist_gap, theta_l(l)+f*dist_gap, Xl(l,3)-h/2, Xl(l,3)+h/2);
            A(l,(k-1)*dim + n)=Vkl;


            % A(l,(k-1)*dim + n)=(1/(4*pi))*((Xl(l,n)-Q(k,n))/norm(Xl(l,:)-Q(k,:),2)^3);
        end
    end

end

for k=1:numSrc*dim
    A(:,k)=A(:,k)-sum(A(:,k))/L;
    A(:,k)=A(:,k)-sum(A(:,k))/L;
    A(:,k)=A(:,k)-sum(A(:,k))/L;
end


end