% 2020-10-15 Normalized Weight Matrix,
% nlmeans_weight_sym_inpaint_Average1(f,h,D0,nwin,nbloc), 
% f: damaged image
% h: smoothing parameter for calculating Gaussian weights
% D0: damaged area
% nwin: patch radius
% nbloc: searching window radius


clear; close all

mex nlmeans_weight_sym_inpaint_Average1.c

%% Barbara for a rectangle damage
g0 = imread('Barbara.bmp');
g0 = double(rgb2gray(g0));
I = g0(144:370,200:382); % Original Image
D0 = ones(size(I));
D0(170:185, 32:120) = 0; % Damaged Area, rectangle1
f = I.*D0 + 255*(1-D0); % Damaged Image
forig = f;

figure, imshow(uint8([f, I, D0*255]),[]);title("Damaged and Original image");
[Nx, Ny] = size(I);
E = ones(Nx*Ny,1);

%% NL weights
sigma1 = 8; % sigma;
nwin = 3; % patch
nbloc = 5; % searching windows
tm = cputime;
w = nlmeans_weight_sym_inpaint_Average1(f,sigma1,D0,nwin,nbloc);
% wsum = w*E;
w = w.';
% w = NormWRow(w, wsum);

tm = cputime - tm;
disp(['weight matrix time: ', num2str(tm), 's']);

%% Alpha_NL
% cof option 1
V = ones(Nx*Ny, 1);
cof = V;
% -------------------------------------------------------------------------
    
COF = reshape(cof, Nx,Ny);
figure,imshow(COF);title('COF')

%% Parameter of NLF
dt = 1;%1e-4;
lambda = 0.75;

%% Initialization
uv = reshape(f, Nx*Ny, 1);
u = f;
D = D0;

iterN = 4;
iteru = [30, 100, 15, 20, 300, 400];%



wE  = w*E;
wEwE = wE.*wE;
we = reshape(wE, Nx, Ny);
% figure,imshow(we,[]); title('w*E')

for i=1:iterN
    %%
    u_yellow(:,:,1) = u.*D;
    u_yellow(:,:,2) = u.*D;
    u_yellow(:,:,3) = u.*D;
    d = 1-D;
    Dy(:,:,1) = d;
    Dy(:,:,2) = d;
    Dy(:,:,3) = 0*d;
    Dy = Dy*255;
    u_yellow = u_yellow+Dy;
    imwrite(uint8(u_yellow), ['yellow-Barbara-rectangle1-NNLF-lambda' num2str(lambda)...
                        '-' num2str(i-1) '.png'])
                    
    %%
    for iu=1:iteru(i)
        %% Fidelity term
        uf = (f - u).*D;
        tmpf = lambda * reshape(uf, Nx*Ny, 1); % E-L fidelity
        
        %% Semi-Implicit
        uv = reshape(u, Nx*Ny, 1);
        wuv = w*uv;
        
        wEuv = wE.*uv;
        NL_lap = wuv - wEuv;%w*uv - w*E.*uv;
        tmp2 = w*(cof.*NL_lap);
        tmp1_1 = cof.*wE.*wuv;
        tmp1_2 = cof.*wEwE; %-
        tmp1 = tmp1_1 - uv.*tmp1_2;
        uv = (uv + dt*(tmp1_1-tmp2 + tmpf))./(1+dt*tmp1_2);
%         uv = uv + dt*(tmp1-tmp2 + tmpf);
%         uv = (uv + dt*(tmp1_1-tmp2))./(1+dt*tmp1_2); % without fidelity
        
        u = reshape(uv, Nx, Ny);
        u = f.*D0 + u.*(1-D0);
        
        PSNR = psnr(u, I);
        dif = u-f;
        disp(['iter: ' num2str(iu) ...
            '---PSNR: ' num2str(PSNR)...
            '---tmp1:' num2str(max(tmp1(:))) ...
            '---tmp2:' num2str(max(tmp2(:)))]);
        figure(113), imshow(uint8([u, f, I]),[]);title("inpainted result and Original image");
        figure(112),imshow(dif,[]);title("difference");
    end
    
    imwrite(uint8(u), ['Barbara-rectangle1-NNLF-lambda' num2str(lambda)...
                        '-' num2str(i) '-' num2str(iu) '.png'])
                  
    %% update weight
    U(:,:,i) = u;
    D_(:,:,i) = D;
    
    f = u;
    D = 1 - double(u==255 & D0~=1);

    w = nlmeans_weight_sym_inpaint_Average1(u,sigma1,D,nwin,nbloc);
    w = w.';
    
    wE  = w*E;
    wEwE = wE.*wE;
    we = reshape(wE, Nx, Ny);
end

figure(113), imshow(uint8([u, forig, I]),[]);title("inpainted result and Original image");

save('NNLF-Barbara-rectangle1-new', 'I', 'forig', 'D0', 'sigma1',...
    'nwin','nbloc','dt','lambda','iterN','iteru', 'u','PSNR');


