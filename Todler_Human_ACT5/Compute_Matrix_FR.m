% Computing the integrals FR to form the matrices that are part of the
% matrix A for the Todler algorithm
% Jennifer Mueller 7/1/2021
% =========================================================================
% We will use Gaussian quadrature to compute the integrals with 10 terms. 

% =========================================================================
% Computes the matrices mxFR11, mxFR12, mxFR121, mxFR122, mxFR141, mxFR142, mxFR143, mxFR144, mxFR14, mxFR3

% =========================================================================

function [mFR11, mFR12, mFR121, mFR122, mFR141, mFR142, mFR143, mFR144, mFR14, mFR3] = Compute_Matrix_FR(N,mr,h)

% N = 32; % Number of Fourier coefficients
% h = 31.4325 ; % height of the cylinder (cm)
% load mr  %  This is the vector of radii in the Joshua tree mesh hard-coded for a radius of 15 cm.  This has 22 rings and therefore 496 voxels
% % The one below is for 11 rings.
% mr = [0,2.694079530;2.694079530,3.810003810;3.810003810,5.388159061;5.388159061,6.599120176;6.599120176,8.082238591;8.082238591,9.332565253;9.332565253,10.43412515;10.43412515,11.74322042;11.74322042,12.92035154;12.92035154,13.99884788;13.99884788,15]; %11 rings
Nr = length(mr);

%%  compute mxFR11 %%
mFR11 = NaN(N,N,Nr);
for ir = 1:Nr
    r1 = mr(ir,1);
    r2 = mr(ir,2);
    for n = 1:N
        for np = 1:N
            mFR11(n,np,ir) = n*np/(n+np)*(r2^(n+np)-r1^(n+np));
        end
    end
%     disp('mxFR11')
%     ir
end
%% %%%%%% set up the nodes and weights to compute the integrals involving r with Gaussian quadrature%%%%%%%%%%%%%%

%% the positive nodes and weights for Gaussian quadrature with n=10
xoc = [0.148874338981631,0.433395394129247,0.679409568299024,0.865063366688985,0.973906528517172]; % the nodes in Gaussian quadrature
woc = [0.295524224714753,0.269266719309996,0.219086362515982,0.149451349150581,0.066671344308688]; % the weights in Gaussian quadrature
% xoc = [0.064056892862605626085,0.191118867473616309159,0.315042679696163374387,0.433793507626045138487,0.545421471388839535658,0.648093651936975569252,0.740124191578554364244,0.820001985973902921954,0.886415527004401034213,0.938274552002732758524,0.974728555971309498198,0.995187219997021360180]; % the nodes in Gaussian quadrature
% woc = [0.127938195346752156974,0.125837456346828296121,0.121670472927803391204,0.115505668053725601353,0.107444270115965634783,0.097618652104113888270,0.086190161531953275917,0.073346481411080305734,0.059298584915436780746,0.044277438817419806169,0.028531388628933663181,0.012341229799987199547]; % the weights in Gaussian quadrature


%% to get the full list of nodes and weights
xoc = [xoc,-xoc]; % the negative and positive nodes are symmetric
woc = [woc,woc]; % the weights for negative node -xoc is same as that for the positive nodes xoc


%% %%%%%% compute matrices mFR122, mFR12 %%%%%%%%%%%%%%
mFR121 = NaN(N,N,N+1,Nr);
mFR122 = NaN(N,N,N+1,Nr);
% j = 1:length(xoc);
for ir = 1:Nr
    r1 = mr(ir,1);
    r2 = mr(ir,2); 
    xoc_new = (r2-r1)/2*xoc+(r2+r1)/2;
    const = (r2-r1)/2;
    for n = 1:N
        for mp = 1:N
            npBessel=besseli(0,mp*pi*xoc_new/h);
            for np = 0:N
                np1Bessel=besseli(np+1,mp*pi*xoc_new/h);
                sum121 = const*sum(woc.*(xoc_new).^n .* np1Bessel);
                sum122 = const*sum(woc.*(xoc_new).^(n-1) .* npBessel);%besseli(np,mp*pi*xoc_new/h));
                npBessel=np1Bessel;
                mFR121(n,mp,np+1,ir) = n*mp*pi/h*sum121;
                mFR122(n,mp,np+1,ir) = n*np*sum122;                
             
            end
        end
    end
%     disp('mxFR12')
%     ir
end

mFR12 = mFR121 + mFR122;

%%
% 
% 
%% %%%%%% compute matrices mFR144, mFR14 & mFR3 %%%%%%%%%%%%%%
mFR141 = NaN(N,N+1,N,N+1,Nr);
mFR142 = NaN(N,N+1,N,N+1,Nr);
mFR143 = NaN(N,N+1,N,N+1,Nr);
mFR144 = NaN(N,N+1,N,N+1,Nr);
mFR3 = NaN(N,N+1,N,N+1,Nr);

for ir = 1:Nr
    r1 = mr(ir,1);
    r2 = mr(ir,2);
    xoc_new = (r2-r1)/2*xoc+(r1+r2)/2;
    const = (r2-r1)/2;
    for m = 1:N
        for n = 0:N
            A_bessel = besseli(n+1,m*pi*xoc_new/h);
            B_bessel = besseli(n,m*pi*xoc_new/h);
            for mp = 1:N
                C_bessel = besseli(0,mp*pi*xoc_new/h);
                for np = 0:N
                    D_bessel = besseli(np+1,mp*pi*xoc_new/h);
                    sum141 = const*sum(woc.*xoc_new.*A_bessel.*D_bessel);
                    sum142 = const*sum(woc.*A_bessel.*C_bessel);
                    sum143 = const*sum(woc.*D_bessel.*B_bessel);
                    sum144 = const*sum((woc./xoc_new).*B_bessel.*C_bessel);
                    sum3   = const*sum(woc.*xoc_new.*B_bessel.*C_bessel);
                    C_bessel = D_bessel;

                    mFR141(m,n+1,mp,np+1,ir) = m*mp*pi^2/h/h*sum141;
                    mFR142(m,n+1,mp,np+1,ir) = m*np*pi/h*sum142;
                    mFR143(m,n+1,mp,np+1,ir) = mp*n*pi/h*sum143;
                    mFR144(m,n+1,mp,np+1,ir) = n*np*sum144;
                    mFR3(m,n+1,mp,np+1,ir) = m*mp*pi^2/h/h*sum3; 
                end
            end
        end
    end
%     disp('mxFR14')
%     ir
end
mFR14 = mFR141 + mFR142 + mFR143 + mFR144;

end