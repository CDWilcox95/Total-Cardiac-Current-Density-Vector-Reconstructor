function A_matrix = Compute_A_Matrix(mFR11, mFR12, mFR122, mFR144, mFR14, mFR3,mFZ11, mFZ12, mFZ13, mFZ14, mFZ3,Anm_xk, Bnm_xk, an_xk, bn_xk, ITheta_1, ITheta_2, ITheta_3, ITheta_4,mr,mz,FN,K,Voxel_N,Joshmesh)

    A_matrix=zeros(K,K,Voxel_N);

    FR11=zeros(FN,FN);FZ11=zeros(FN,FN);FTheta11=zeros(FN,FN);FTheta21=zeros(FN,FN);
    FR12=zeros(FN,FN,FN+1);FR122=zeros(FN,FN,FN+1);FZ12=zeros(FN,FN,FN+1);FTheta12=zeros(FN,FN,FN+1);FTheta22=zeros(FN,FN,FN+1);
    FR13=zeros(FN,FN,FN+1);FR132=zeros(FN,FN,FN+1);FZ13=zeros(FN,FN,FN+1);FTheta13=zeros(FN,FN,FN+1);FTheta23=zeros(FN,FN,FN+1);
    FR14=zeros(FN,FN+1,FN,FN+1);FR144=zeros(FN,FN+1,FN,FN+1);FR3=zeros(FN,FN+1,FN,FN+1);FZ14=zeros(FN,FN+1,FN,FN+1);FZ3=zeros(FN,FN+1,FN,FN+1);FTheta14=zeros(FN,FN+1,FN,FN+1);FTheta24=zeros(FN,FN+1,FN,FN+1);

    for voxel=1:Voxel_N
        tic;
        r1=Joshmesh(voxel,1);
        r2=Joshmesh(voxel,2);
        th1=Joshmesh(voxel,3);
        th2=Joshmesh(voxel,4);
        Z1=Joshmesh(voxel,5);
        Z2=Joshmesh(voxel,6);
        fprintf( '\n %g  ', voxel);

        [minr,index_mr] = min(abs(r1-mr(:,1))); % find the index of r for this voxel in mr
        [minz,index_mz] = min(abs(Z1-mz(:,1))); % find the index of z for this voxel in mz
        FR11(1:FN,1:FN) = mFR11(1:FN,1:FN,index_mr); %retrieve FR11 from mFR11
        FZ11(1:FN,1:FN) = mFZ11(1:FN,1:FN,index_mz); %retrieve FZ11 from mFZ11
         FR12(1:FN,1:FN,1:FN+1) = mFR12(1:FN,1:FN,1:FN+1,index_mr); %retrieve FR12 from mFR12
        FR122(1:FN,1:FN,1:FN+1) =mFR122(1:FN,1:FN,1:FN+1,index_mr); %retrieve FR122 from mFR122
         FZ12(1:FN,1:FN,1:FN+1) = mFZ12(1:FN,1:FN,1:FN+1,index_mz); %retrieve FZ12 from mFZ12
         FR13(1:FN,1:FN,1:FN+1) =  mFR12(1:FN,1:FN,1:FN+1,index_mr); %retrieve FR13 from mFR13 (note: mFR13=mFR12)
        FR132(1:FN,1:FN,1:FN+1) = mFR122(1:FN,1:FN,1:FN+1,index_mr); %retrieve FR132 from mFR132 (note: mFR132=mFR122)
         FZ13(1:FN,1:FN,1:FN+1) =  mFZ13(1:FN,1:FN,1:FN+1,index_mz); %retrieve FZ13 from mFZ13 
         FR14(1:FN,1:FN+1,1:FN,1:FN+1) = mFR14(1:FN,1:FN+1,1:FN,1:FN+1,index_mr); %retrieve FR14 from mFR14
        FR144(1:FN,1:FN+1,1:FN,1:FN+1)= mFR144(1:FN,1:FN+1,1:FN,1:FN+1,index_mr); %retrieve FR144 from mFR144
          FR3(1:FN,1:FN+1,1:FN,1:FN+1) =  mFR3(1:FN,1:FN+1,1:FN,1:FN+1,index_mr); %retrieve FR3 from mFR3
         FZ14(1:FN,1:FN+1,1:FN,1:FN+1) = mFZ14(1:FN,1:FN+1,1:FN,1:FN+1,index_mz); %retrieve FZ14 from mFZ14
          FZ3(1:FN,1:FN+1,1:FN,1:FN+1) =  mFZ3(1:FN,1:FN+1,1:FN,1:FN+1,index_mz); %retrieve FZ3 from mFZ3

        for jx=1:K            

    %        matlabpool open 4;
    %        parfor jk=1:K
            for jk=1:jx

                % *********** 1st summation term ********

                for n=1:FN
                    for np=1:FN
                       FTheta11(n,np)= an_xk(n+1,jx)*an_xk(np+1,jk)*ITheta_1(voxel,n+1,np+1) + ... 
                                       bn_xk(n+1,jx)*bn_xk(np+1,jk)*ITheta_2(voxel,n+1,np+1) + ...
                                       an_xk(n+1,jx)*bn_xk(np+1,jk)*ITheta_3(voxel,n+1,np+1) + ...
                                       bn_xk(n+1,jx)*an_xk(np+1,jk)*ITheta_4(voxel,n+1,np+1);
                       FTheta21(n,np)= bn_xk(n+1,jx)*bn_xk(np+1,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                       an_xk(n+1,jx)*an_xk(np+1,jk)*ITheta_2(voxel,n+1,np+1) - ...
                                       bn_xk(n+1,jx)*an_xk(np+1,jk)*ITheta_3(voxel,n+1,np+1) - ...
                                       an_xk(n+1,jx)*bn_xk(np+1,jk)*ITheta_4(voxel,n+1,np+1); % change last minus sign (Felix) 12062011
                    end                
                end
                Term_1=sum(sum(FR11.*FZ11.*(FTheta11+FTheta21)));
                % ******* end of 1st term *******


                % *********** 2nd summation term ********
                for n=1:FN
                    for mp=1:FN
                        for np=0:FN      
                            FTheta12(n,mp,np+1)=an_xk(n+1,jx)*Anm_xk(np+1,mp,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                                bn_xk(n+1,jx)*Bnm_xk(np+1,mp,jk)*ITheta_2(voxel,n+1,np+1) + ...
                                                an_xk(n+1,jx)*Bnm_xk(np+1,mp,jk)*ITheta_3(voxel,n+1,np+1) + ...
                                                bn_xk(n+1,jx)*Anm_xk(np+1,mp,jk)*ITheta_4(voxel,n+1,np+1);
                            FTheta22(n,mp,np+1)=bn_xk(n+1,jx)*Bnm_xk(np+1,mp,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                                an_xk(n+1,jx)*Anm_xk(np+1,mp,jk)*ITheta_2(voxel,n+1,np+1) - ...
                                                bn_xk(n+1,jx)*Anm_xk(np+1,mp,jk)*ITheta_3(voxel,n+1,np+1) - ...
                                                an_xk(n+1,jx)*Bnm_xk(np+1,mp,jk)*ITheta_4(voxel,n+1,np+1);
                        end
                    end                
                end
                Term_2=sum(sum(sum(FR12.*FZ12.*FTheta12 + FR122.*FZ12.*FTheta22)));
                % ******* end of 2nd term *******

                % *********** 3rd summation term ********

                for np=1:FN
                    for m=1:FN
                        for n=0:FN                        
                            FTheta13(np,m,n+1)=Anm_xk(n+1,m,jx)*an_xk(np+1,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                               Bnm_xk(n+1,m,jx)*bn_xk(np+1,jk)*ITheta_2(voxel,n+1,np+1) + ...
                                               Anm_xk(n+1,m,jx)*bn_xk(np+1,jk)*ITheta_3(voxel,n+1,np+1) + ...
                                               Bnm_xk(n+1,m,jx)*an_xk(np+1,jk)*ITheta_4(voxel,n+1,np+1);
                            FTheta23(np,m,n+1)=Bnm_xk(n+1,m,jx)*bn_xk(np+1,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                               Anm_xk(n+1,m,jx)*an_xk(np+1,jk)*ITheta_2(voxel,n+1,np+1) - ...
                                               Bnm_xk(n+1,m,jx)*an_xk(np+1,jk)*ITheta_3(voxel,n+1,np+1) - ...
                                               Anm_xk(n+1,m,jx)*bn_xk(np+1,jk)*ITheta_4(voxel,n+1,np+1); 
                        end
                    end
                end
                Term_3=sum(sum(sum(FR13.*FZ13.*FTheta13 + FR132.*FZ13.*FTheta23)));
                % ******* end of 3rd term *******


                % *********** 4th summation term ********

                for m=1:FN
                    for n=0:FN
                        for mp=1:FN
                            for np=0:FN                        
                                FTheta14(m,n+1,mp,np+1)=Anm_xk(n+1,m,jx)*Anm_xk(np+1,mp,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                                        Bnm_xk(n+1,m,jx)*Bnm_xk(np+1,mp,jk)*ITheta_2(voxel,n+1,np+1) + ...
                                                        Anm_xk(n+1,m,jx)*Bnm_xk(np+1,mp,jk)*ITheta_3(voxel,n+1,np+1) + ...
                                                        Bnm_xk(n+1,m,jx)*Anm_xk(np+1,mp,jk)*ITheta_4(voxel,n+1,np+1);
                                FTheta24(m,n+1,mp,np+1)=Bnm_xk(n+1,m,jx)*Bnm_xk(np+1,mp,jk)*ITheta_1(voxel,n+1,np+1) + ...
                                                        Anm_xk(n+1,m,jx)*Anm_xk(np+1,mp,jk)*ITheta_2(voxel,n+1,np+1) - ...
                                                        Bnm_xk(n+1,m,jx)*Anm_xk(np+1,mp,jk)*ITheta_3(voxel,n+1,np+1) - ...
                                                        Anm_xk(n+1,m,jx)*Bnm_xk(np+1,mp,jk)*ITheta_4(voxel,n+1,np+1);
                            end
                        end
                    end                
                end
                Term_4=sum(sum(sum(sum(FR14.*FZ14.*FTheta14 + FR144.*FZ14.*FTheta24 + FR3.*FZ3.*FTheta14)))) ;
                % ******* end of 4th term *******

                A_matrix(jx,jk,voxel)=Term_1 + Term_2 + Term_3 + Term_4 ;
                A_matrix(jk,jx,voxel)=A_matrix(jx,jk,voxel);

            end % jk loop for current pattern
    %         fprintf('%g ,', jx);

        end % jx loop for current pattern    
      toc;
    %   xx =A_matrix(:,:,1);
    %   xx =xx;
    end % voxel loop

end