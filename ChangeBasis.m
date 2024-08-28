function [X, basis_coef]=ChangeBasis(I0, I)
[L, numCurrents]=size(I0);

%% Populate X Matrix
B=zeros(size(I0));
for i=1:numCurrents
   B(:,i)=I0(:,i); 
    
end

%% Get unknowns for current pattern k.
basis_coef=zeros(numCurrents,numCurrents);
newCurrents=zeros(L,numCurrents);
for k=1:numCurrents
   x=B\I(:,k);
   f=0;
   for i=1:numCurrents
       f=f+x(i)*I0(:,i);
   end
   newCurrents(:,k)=f;
   basis_coef(:,k)=x;
end

X=zeros(numCurrents^2, numCurrents^2);

for i=1:numCurrents
    x_i=basis_coef(:,i);
    for j=1:numCurrents
        X_min=zeros(numCurrents,numCurrents);
        
        X_min=basis_coef';
        
        X_min=x_i(j).*X_min;
        
        X(1+(i-1)*numCurrents:(i)*numCurrents, 1+(j-1)*numCurrents:(j)*numCurrents)=X_min;
    end
    
end





X( :, ~any(X,1) ) = [];  %columns
X(~any(X,2),:)=[];

end