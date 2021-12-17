%%%%%%%%%%%%%%%%%%%%% NL-TV (NL-BiHarmonic)+ Denoising %%%%%%%%%%%%%%%%%%%%
% u(t)= NL_TV(NL_FourthOrder)(u)+lambda*(u-f)
% ****** BL biharmonic ******
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
addpath('../evaluate')

% mex nlmeans_weight_sym_topmbiggest.cpp

sigma=20;
load BarbarbaPart-20.mat
[M1 N1]=size(g1);

m=15;   %% half size of search window
bb=5;   %% half size of patch

fv = reshape(f1, [M1*N1, 1]);
gv = reshape(g1, [M1*N1, 1]);

%% Compute Nonlocal Weights
h=sigma/sqrt(2)*bb;%(2*bb+1);
a=1; c=1;

time = cputime;
w = nlmeans_weight_sym_topmbiggest(f1,h,bb,m);
w2 = w.';
w = max(w, w2);
time=cputime - time;
disp(['weight version 2, time of calculating weight: ' num2str(time)])

I=ones(size(w,1),1);
NL_grad=I;
u=f1;

%% Prepare for the NLFourthOrder
sig = 1;
siz = 2*bb+1;%sig*6;
f2 = imgaussian(f1, sig, siz);
fv2 = reshape(f2, [M1*N1, 1]);

K = 30;
K = K^2;
f_NL_grad=abs(w*(fv2.*fv2)-2*(w*fv2).*fv2+(w*I).*(fv2.*fv2));
cof = 1./(1 + f_NL_grad/K);
COF=reshape(cof, M1,N1);
figure,imshow(COF)

%% Main for NLTV/NLF
disp('-------- BiHarmonic --------');
MaxiterU=200;
dtU= 1e-5; % semi-implicit
lambda=20;

Energy = [];

uv=fv;
for IterU=1:MaxiterU
    %     %% Explicit
    %     NL_lap = w*uv - w*I.*uv;
    %     tmp1 = w*I.*cof.*NL_lap;
    %     tmp2 = w*(NL_lap.*cof);
    %     NL_4order = tmp1 - tmp2;
    %     uv=uv+dtU*(NL_4order+lambda*(fv-uv));
    
    %% Semi-Implicit
    tmpf = lambda.*(f1-u);
    tmpfv = reshape(tmpf, M1*N1, 1);
    NL_lap = w*uv - w*I.*uv;
    tmp2 = w*(cof.*NL_lap);
    tmp1_1 = cof.*(w*I).*(w*uv);
    tmp1_2 = cof.*(w*I).*(w*I);
    uv = (uv + dtU*(tmp1_1-tmp2 + tmpfv))./(1+dtU*tmp1_2);
    
    
    u=reshape(uv,[M1, N1]);
    PSNR(IterU) = psnr(u, g1);
    SNR(IterU) = snr(u, g1);
    MAE(IterU) = mae(u, g1);
    RMSE(IterU)=sqrt(sum((uv-gv).^2)/(M1*N1));
    
    %% Energy
    NL_lap = w*uv - w*I.*uv;
    Reg = cof .* (NL_lap.^2) / 2;
    reg = sum(Reg);
    Fid = (uv - fv).^2 .* lambda / 2;
    fid = sum(Fid);
    energy = reg + fid;
    if IterU~= 1
        dif_eng = energy - Energy(end);
    end
    
    Energy = [Energy energy];
    
    if IterU~= 1
        disp(['iter: ' num2str(IterU) '   energy difference: ' num2str(dif_eng)])
    end
end
[pv, pi] = max(PSNR);
[sv, si] = max(SNR);
[mv, mi] = min(MAE);
[rv, ri] = min(RMSE);
disp(['maxPSNR: ' num2str(pv) ' - iter: ' num2str(pi)]);
disp(['maxSNR: ' num2str(sv) ' - iter: ' num2str(si)]);
disp(['minMAE: ' num2str(mv) ' - iter: ' num2str(mi)]);
disp(['minRMSE: ' num2str(rv) ' - iter: ' num2str(ri)]);

figure,imshow(uint8([g1, f1, u]));title('Orig, Noisy, Denoised');
figure;plot(PSNR);title('PSNR between denoised and original');

