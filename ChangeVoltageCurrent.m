function V_new=ChangeVoltageCurrent(basis_coef, V0)

    
    [L, numCurrents]=size(V0);
    V_new=zeros(L,numCurrents);
    
    for alpha=1:numCurrents
        x_alpha=basis_coef(:,alpha)';
        for l=1:L
            Vl_alpha=dot(x_alpha,V0(l,:));
            V_new(l,alpha)=Vl_alpha;
        end
    end


end