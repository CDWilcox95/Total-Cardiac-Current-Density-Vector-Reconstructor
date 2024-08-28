function [C, vol_tetr]=CalculateMeshCoef(T, Pt)

[M, elem_nodes]=size(T);
if elem_nodes==4
    order=1;
else
    order=2;
end
vol_tetr=zeros(M,1);
C=zeros(M, elem_nodes,elem_nodes);  %rows=# elements, columns=#coefficients, slices=#coefficients
switch order
    case 1
        for t=1:M
            t_ind=T(t,:);
            pts=Pt(t_ind,:);
            B=[1, pts(1,:);1, pts(2,:);1, pts(3,:);1, pts(4,:)];
            vol_tetr(t)=abs((1/6)*det(B));

            coef_t=inv(B)';
            C(t,:,1)=coef_t(1,:);
            C(t,:,2)=coef_t(2,:);
            C(t,:,3)=coef_t(3,:);
            C(t,:,4)=coef_t(4,:);
        end

    case 2
        for t=1:M
            t_ind=T(t,:);
            pts=Pt(t_ind,:);
            B=zeros(10,10);
            for i=1:10
                B(i,:)=[1, pts(i,1)^2, pts(i,2)^2, pts(i,3)^2, pts(i,1), pts(i,2), pts(i,3), pts(i,1)*pts(i,2), pts(i,1)*pts(i,3), pts(i,2)*pts(i,3)];
            end

            G=[1, pts(1,:);1, pts(2,:);1, pts(3,:);1, pts(4,:)];

            vol_tetr(t)=abs((1/6)*det(G));

            coef_t=inv(B)';
            C(t,:,1)=coef_t(1,:);
            C(t,:,2)=coef_t(2,:);
            C(t,:,3)=coef_t(3,:);
            C(t,:,4)=coef_t(4,:);
            C(t,:,5)=coef_t(5,:);
            C(t,:,6)=coef_t(6,:);
            C(t,:,7)=coef_t(7,:);
            C(t,:,8)=coef_t(8,:);
            C(t,:,9)=coef_t(9,:);
            C(t,:,10)=coef_t(10,:);

        end



end


end