function U0 = forward_gap_model3D(L,R0,h,FN,Theta_L,Zlz,Anm_xk, Bnm_xk, an_xk, bn_xk)
% Computed from Eq. 4 of Felix's Paper
format long

r = R0;

term1 = zeros(L,L-1);
term2 = zeros(L,L-1);
for l=1:L
   
    for k =1:L-1
        sum_term = 0;
        for n=1:FN
            
            sum_term = sum_term + (r^n)*(an_xk(n+1,k)*cos(n*Theta_L(l))+bn_xk(n+1,k)*sin(n*Theta_L(l)));
            
        end
       
        term1(l,k) = sum_term;
        
    end
       
end

for l=1:L
   
    for k =1:L-1
        sum_term_2 = 0;
        for m=1:FN
            
            sum_term = 0;
            for n=0:FN

                sum_term = sum_term + besseli(n,m*pi*r/h)*cos(m*pi*Zlz(l)/h)*(Anm_xk(n+1,m,k)*cos(n*Theta_L(l))+Bnm_xk(n+1,m,k)*sin(n*Theta_L(l)));
            end
            sum_term_2 = sum_term_2 + sum_term;
            
        end
        
        term2(l,k) = sum_term_2;
        
    end
       
end

U0 = term1+term2;

end