% 2021-12-15, Ying Wen, for JMIV experimental results
clear; close all

% mex nlmeans_weight_sym_topmbiggest.cpp
addpath('../evaluate/')

ImgName = 'lena256';
sigma = 7;
sigmag1 = 9;
sigmag2 = 2;% or 2
load(['res-JMIV/' ImgName '-' num2str(sigma) ...
      '-' num2str(sigmag1) '-' num2str(sigmag2)]);
h = K;
g = padarray(g,size(K),'symmetric');
f = padarray(f,size(K),'symmetric');

%% start to deblur
sigma1 = 10;%sigma;
% compute the weight
nwin = 4; %3*3 patch
nbloc = 12; %21*21 searching windows
tm = cputime;
w = nlmeans_weight_sym_topmbiggest(f,sigma1,nwin,nbloc);
w2 = w.';
w = max(w, w2);
tm = cputime - tm;
disp(['weight matrix time: ', num2str(tm), 's']);

preImage = f;

%% NLF debluring
lambda = 2000;
niter = 5000;
K = 40;
[uNLF, energy, SNR, PSNR, MAE, SSIM, RMSE] = NLABH_decon...
    (g, f, preImage, h, w, K, lambda, niter);
uNLF_ = uNLF(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2));
f_ = f(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2));
g_ = g(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2));

figure,imshow([g_, f_, uNLF_],[]); title('Orig,    Degenerated,   Reconstructed');

disp(['lambda:' num2str(lambda) ' K: ' num2str(K)]);
[v,ind]=max(PSNR);
disp(['Max PSNR: ' num2str(v), '  ind: ', num2str(ind)]);
disp(['PSNR: ' num2str(PSNR(end)) '  SNR: ' num2str(SNR(end)) ...
      '  MAE: ' num2str(MAE(end)) '  SSIM: ' num2str(SSIM(end))]);
disp('---------------- end ----------------')
disp('   ')

figure,
subplot(1,2,1),plot(PSNR);title('PSNR')
subplot(1,2,2),plot(energy);title('energy')

function [u, energy, SNR, PSNR, MAE, SSIM, RMSE] = NLABH_decon...
    (Io, f, u0, h, w, K, lambda, iters)

% Copyright(C)Ying Wen, 2020
% E(u) = \int_{\Omega}\alpha_(NL)(u_0)|\Laplace_{NL}(u)| + \lambda\int_{\Omega}(Ku-f)^2

% Output
% u: decon result image
% energy: E for each iteration
% SNR, PSNR, MAE, SSIM: measurement

% Input
% Io: Latent clear image
% f: Original blurry image
% u0: PreImage
% h: Parameter of blurring
% w: Nonlocal weight
% K: Parameter of ABH
% lambda: \lambda
% iters: Max number of iteration


dt = 0.0001; % Barbara, Barbara2, lena256
% dt = 0.001 % PDE_seismic
eps = 1e-7;
eps2 = 1e-6;
[M,N] = size(f);
MN = M*N;
I = ones(size(w,1),1);
energy = [];
SNR = [];
PSNR = [];
MAE = [];
SSIM = [];
RMSE = [];
k=1;

u = u0; % Initional

%% for alpha % only once
u0v = reshape(u0, MN, 1);
K = K^2;
f_NL_grad=abs(w*(u0v.*u0v)-2*(w*u0v).*u0v+(w*I).*(u0v.*u0v));
cof = 1./(1 + f_NL_grad/K);
COF=reshape(cof, M,N);
figure,imshow(COF);title('COF')

%% Initional energy
NL_lap2 = w*u0v - w*I.*u0v;
E1 = cof.*(NL_lap2.^2);
e1 = sum(E1(:));
difu = f - imfilter(u0,h,'circular');
E2 = difu.^2;
e2 = lambda*sum(E2(:))/2;
energy0 = (e1 + e2)/MN;
energy = [energy, energy0];

v = reshape(u, MN, 1);
u_old = u;
v_old = v;

egy = energy0 + 1;

while ((k<iters) & (dt>eps2) & (abs(egy-energy0)>eps))
    %% for fidelity term
    uf = f - imfilter(u,h,'circular');
    uf = imfilter(uf,h,'circular');
    tmpf = lambda * reshape(uf, MN, 1); % E-L fidelity
    
    %% for ABH term % Semi-Implicit
    NL_lap = w*v - w*I.*v;
    tmp2 = w*(cof.*NL_lap);
    tmp1_1 = cof.*(w*I).*(w*v);
    tmp1_2 = cof.*(w*I).*(w*I); %-
    v = (v + dt*(tmp1_1-tmp2 + tmpf))./(1+dt*tmp1_2);
    
    %% for energy
    u = reshape(v, M, N);
    NL_lap2 = w*v - w*I.*v;
    E1 = cof.*(NL_lap2.^2);
    e1 = sum(E1(:));
    difu = f - imfilter(u,h,'circular');
    E2 = difu.^2;
    e2 = lambda*sum(E2(:))/2;
    egy = (e1 + e2)/MN;

%     if egy < energy0
        energy0 = energy(end);
        energy = [energy, egy];
        k = k+1;
        u_old = u;
        v_old = v;
%     else
%         u = u_old;%old
%         v = v_old;%old
%         dt = .8*dt;
%     end
%     disp(['dt: ' num2str(dt) ' dif_energy:' num2str(abs(egy-energy0)) ...
%         ' k: ' num2str(k)]);
    

    %% for measurement
    Psnr = psnr(u(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)),...
        Io(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)));
    Snr = snr(u(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)),...
        Io(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)));
    Mae = mae(u(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)),...
        Io(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)));
    Ssim = ssim(u(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)),...
        Io(1+size(h,1):end-size(h,1),1+size(h,2):end-size(h,2)));
    Rmse=sqrt((sum(sum((u-Io).^2)))/(MN));
    PSNR = [PSNR, Psnr];
    SNR = [SNR, Snr];
    MAE = [MAE, Mae];
    SSIM = [SSIM, Ssim];
    RMSE = [RMSE, Rmse];
end
disp(['dt: ' num2str(dt)])
end