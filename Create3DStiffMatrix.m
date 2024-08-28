function [S, Sk]=Create3DStiffMatrix(C, vols_tetr, cond_vals, T, Pt, file_ext)

[M, elem_nodes]=size(T);
if elem_nodes==4
    order=1;
else
    order=2;
    grad_phi2=@(x,y,z, coef)[2*coef(2)*x+coef(5)+coef(8)*y+coef(9)*z, 2*coef(3)*y+coef(6)+coef(8)*x+coef(10)*z, 2*coef(4)*z+coef(7)+coef(9)*x+coef(10)*y]; 
    brcntrs=ComputeBarycenters(Pt, T);
end

N=length(Pt);

S=sparse(N,N);  Sk=zeros(elem_nodes, elem_nodes, M);
fprintf("Computing order %d Stiffness Matrix... \n", order);

line_length=fprintf("Number of elements integrated:  %d/%d", 0, M);

for t=1:M
    t_ind=T(t,:);
    

    coef_t=C(t,:,:);    coef_t=reshape(coef_t, [elem_nodes,elem_nodes]);  
    sigma_t=cond_vals(t);   Vol_t=vols_tetr(t);
    for i=1:elem_nodes
        
        c_i=coef_t(:,i);
        for j=1:elem_nodes

            c_j=coef_t(:,j);
            Sk(i,j,t)=dot(c_i(2:4), c_j(2:4))*Vol_t;
            if order==1
                S(t_ind(i), t_ind(j))=S(t_ind(i), t_ind(j))+cond_vals(t)*dot(c_i(2:4), c_j(2:4))*Vol_t;
            else
                S(t_ind(i), t_ind(j))=S(t_ind(i), t_ind(j))+cond_vals(t)*(6*Vol_t)*SIntegral2Order(c_i, c_j, Pt(t_ind(1:4),:), brcntrs(t,:));
            end
        end
    end
    fprintf(repmat('\b',1,line_length));
    line_length=fprintf("Number of Elements Integrated Over:  %d/%d", t, M);
end


fprintf("\nStiffness matrix computed! \n\n");
save(file_ext+"FEMElemStiffMatrix.mat", 'Sk');
save(file_ext+"FEMStiffMatrix.mat", 'S');

end