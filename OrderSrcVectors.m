function M=OrderSrcVectors(num_src, M)
[~, dim, num_frames]=size(M);
M0=zeros(num_src,dim);

for k=1:num_src

    M0(k,:)=M(dim*k-(dim-1):dim*k);
end

M=M0;

end