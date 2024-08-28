function e=CC_RelativeError(U, V)
[L, num_frames]=size(U);

e=0;
num=0;  den=0;
c=zeros(num_frames,1);
for s=1:num_frames
%     c(s)=dot(U(:,s), V(:,s))/norm(V(:,s),2)^2;

    for l=1:L
        num=num+abs(U(l,s)-V(l,s))^2;
        den=den+abs(V(l,s))^2;
    end
    
end

e=sqrt(num/den);

end