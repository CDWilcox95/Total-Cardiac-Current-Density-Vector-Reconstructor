function [x, reg_param]=RegSystem(A, b, method, param)
[num_r, num_c]=size(A);
[U,S,V]=svd(A); sing_vals=diag(S);
switch method
    case 'LM'
        if isempty(param)
            param=0.05;
        end
        reg_param_1=5;
        reg_param_2=param;
        x=inv((A'*A)+reg_param_1.*eye(size(A'*A))+reg_param_2.*diag(A'*A).*eye(size(A'*A)))*((A')*b);

        reg_param=param;
    case 'tSVD'
        LCurve_param=l_curve(U,sing_vals,b,'tsvd');     reg_param=LCurve_param;
        param=LCurve_param;
        UT=U(:,1:param);  ST=S(1:param, 1:param);  VT=V(:,1:param);
        K=VT*inv(ST)*UT';
        x=K*b;     %% Solve for Lagrange Mult.
    case 'Lasso'
        [x, ~] = CD_Lasso(A,b,zeros(num_c, 1),param,500);
        reg_param=0.1;
    case 'l2'
        lambda=l_curve(U,sing_vals,b,'tikh');   reg_param=lambda;
        x=zeros(num_c,1);
        for i=1:length(sing_vals)
            fi=(sing_vals(i)^2)/(sing_vals(i)^2 +lambda^2);
            x=x+fi*(dot(U(:,i),b)/sing_vals(i))*V(:,i);
        end
    case 'manual_l2'
        x=zeros(num_c,1);
        for i=1:length(sing_vals)
            fi=(sing_vals(i)^2)/(sing_vals(i)^2 +param^2);
            x=x+fi*(dot(U(:,i),b)/sing_vals(i))*V(:,i);
        end
        reg_param=param;
    case 'manual_tSVD'
        UT=U(:,1:param);  ST=S(1:param, 1:param);  VT=V(:,1:param);
        K=VT*inv(ST)*UT';
        x=K*b; 
        reg_param=param;
    case 'none'
        x=A\b;  reg_param=[];
end


end