% CNV2D_FFBTvsFFT_Test.m 
% This code compares convolution of real functions f,g using iFFBT vs. iFFT. 
%%% Inputs:
% (1) given real functions f,g in the form of a function bundle OR 
% finite values on sampling/visualization grid. 
% Hint: If f,g are given as sampled values then the function handle should 
% be replaced in the related steps, accordingly. 
% (2) the radius rf which f is well supported in B_rf
% (3) the radius rg which g is well supported in B_rg
% (3) the Fourier-Bessel approximation order (M,N)
% (4) grid size for visualization Lv*Lv
%%% Outputs: 
% (1) Approximated values of f*g using iFFBT of order K_MN on 
%     the visualization grid of size Lv*Lv
% (2) Approximated values of f*g using iFFT2 of order K_MN on 
%     the visualization grid of size Lv*Lv
% (3) Surface plots comparison of f*g approximated using iFFBT vs. iFFT2.
%%=========================================================================
% Author: A-GF, 2023. Version 1
%%=========================================================================
% Hint: This test code uses the following zero solver in Lines 50 and 123;
% Jason Nicholson, Bessel Zero Solver, MATLAB Central File Exchange.  
% (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
%% 
close all
clear all
%--------------------------------------------------------------------------
% Hint: The functions should be well-supported in B
%--------------------------------------------------------------------------
%%%----------------Initiating functions f,g--------------------------------
% Example 6.1
rf = 1; % radius covering support of f
f = @(x,y) rectpuls(sqrt(x.^2 + y.^2)./(2.*rf));
rg = 1; % radius covering support of g
g = @(x,y) rectpuls(sqrt(x.^2 + y.^2)./(2.*rg));
% Example 6.2 
%rf = 1; % radius covering support of f
%f = @(x,y) rectpuls(sqrt(x.^2 + y.^2)./(2.*rf));
%rg = 2; % radius covering support of g
%g = @(x,y) rectpuls(sqrt(x.^2 + y.^2)./(2.*rg));
%% Computing the disk covering both supports of f,g, and f*g
a =  ceil(rf + rg);% smallest radius covering support of f*g
b = 1./a;
%% Initiating iFFBT approximation orders (M,N)
M = 10;
N = 10;
%% Computing minimum Fourier sampling order Ks of iFFBT based on (M,N) 
Z = zeros(M+1,N);
for t = 1:M+1
     Z(t,:) = besselzero(t-1,N,1)./pi; % See Line 21.
end
CZ = ceil(Z);
K_MN = max(max(CZ));
Ks = K_MN; % Ks should not be less than K_MN in SMNK 
Ls = 2.*Ks+1; % finite Fourier sampling size
%% Generating Hc-matrix (Hc=H*)
Hc = zeros(Ls,Ls,2*M+1,N);
for t = M+1:2*M+1
    for n = 1:N
        for l = 1:Ks+1
            for k = 1:Ks+1
            lk = [l-1,k-1];
            Hc(l,k,t,n) = CoE(lk,t-1-M,n)./a;
            end
        end
        for l = Ks+2:Ls
            for k = 1:Ks+1
            lk = [l-Ls-1,k-1];
            Hc(l,k,t,n) = CoE(lk,t-1-M,n)./a;
            end
        end
        for l = 1:Ks+1
            for k = Ks+2:Ls
            lk = [l-1,k-Ls-1];
            Hc(l,k,t,n) = CoE(lk,t-1-M,n)./a;
            end
        end
        for l = Ks+2:Ls
            for k = Ks+2:Ls
            lk = [l-Ls-1,k-Ls-1];
            Hc(l,k,t,n) = CoE(lk,t-1-M,n)./a;
            end
        end
    end
end
% negative m
for t = 1:M
    for n = 1:N
        for l = 1:Ks+1
            for k = 1:Ks+1
            lk = [l-1,k-1];
            Hc(l,k,t,n) = (-1)^(t-1-M).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = Ks+2:Ls
            for k = 1:Ks+1
            lk = [l-Ls-1,k-1];
            Hc(l,k,t,n) = (-1)^(t-1-M).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = 1:Ks+1
            for k = Ks+2:Ls
            lk = [l-1,k-Ls-1];
            Hc(l,k,t,n) = (-1)^(t-1-M).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = Ks+2:Ls
            for k = Ks+2:Ls
            lk = [l-Ls-1,k-Ls-1];
            Hc(l,k,t,n) = (-1)^(t-1-M).*CoE(lk,t-1-M,n)./a;
            end
        end
    end
end
%% Generating visualization grid of size Lv*Lv
Lv = Ls + 21; % 
dLv = 2.*a./Lv;
VLv = 1:Lv;
xLv = -a + (VLv-1).*dLv;
%% Generating Fourier-Bessel basis values on the visualization grid Lv*Lv
P = zeros(2*M+1,N,Lv,Lv);
for t = 1:2*M+1
    ZtN = besselzero(t-1-M,N,1); % See Line 21.
    for n = 1:N
        ztn = ZtN(n);
        Dmna = a.*sqrt(pi).*abs(besselj(t-M,ztn));
        Nmna = 1./Dmna;
        for i = 1:Lv
            for j = 1:Lv
                [tij,rij] = cart2pol(xLv(i),xLv(j));
                Jtnarij = besselj(t-1-M,b.*ztn.*rij);
                Etij = exp(1i.*(t-1-M).*tij);
                P(t,n,i,j) = Nmna.*Etij.*Jtnarij.*rectpuls(b.*rij./2);
            end
        end
    end
end
%% Generating K*-matrix on the visualization grid of size Ls*Ls
KcM = zeros(Lv,Lv,Ls,Ls);
for i = 1:Lv
    for j = 1:Lv
        for l = 1:Ls
            for k = 1:Ls
                Hlk = squeeze(Hc(l,k,:,:));
                Pij = squeeze(P(:,:,i,j));
                PijT = transpose(Pij);
                Tijlk = trace(Hlk*PijT);
                KcM(i,j,k,l) = Tijlk;
            end
        end
    end
end
%% Generating F,G values on Fourier sampling grid of size Ls*Ls  
dLs = 2.*a./Ls;
VLs = 1:Ls;
xLs = -a + (VLs-1).*dLs;
[XLs,YLs] = meshgrid(xLs,xLs);
% Generating function values F on Fourier sampling grid
Fs = f(XLs,YLs);
% Generating function values G on Fourier sampling grid
Gs = g(XLs,YLs);
% Generating DFT matrices of Fs,Gs
HFs = fft2(Fs);
HGs = fft2(Gs);
% Generating Hadamard product of DFT matrics
HFG = HFs.*HGs;
%% Generating SMNK[F,G] on the visualization grid of size Lv*Lv
SMNK_FG = zeros(Lv,Lv);
for i = 1:Lv
    for j = 1:Lv
        Kij = squeeze(KcM(i,j,:,:));
        SMNK_FG(i,j) = dLs.^4.*trace(Kij*HFG); 
    end
end
%% ------------------------------------------------------------------------
%========== Computing Convolutional iFFT2 of order K_[M,N] ================
% Generating iFFT2 of order Ks of F*G
Kx = Ks;
Ky = Ks;
Lx = Ls;
Ly = Ls;
Ny = Lv;
Nx = Lv;
QHFG = zeros(Ny,Nx);
% 1
for t = 1:Ky+1
     for s = 1:Kx+1
         QHFG(t,s) = (-1)^(t+s).*HFG(t,s);
     end
 end
% 2 
for t = 1:Ks+1
    for s = Nx-Ks+1:Nx
        QHFG(t,s) = (-1)^(t+s-Nx).*HFG(t,Lx+s-Nx);
    end
end
% 3
for t = Ny-Ky+1:Ny
     for s = 1:Kx+1
         QHFG(t,s) = (-1)^(t+s-Ny).*HFG(Ly+t-Ny,s);
    end
end
% 4
for t = Ny-Ky+1:Ny
    for s = Nx-Kx+1:Nx
        QHFG(t,s) = (-1)^(t+s-Nx-Ny).*HFG(Ly+t-Ny,Lx+s-Nx);
     end
end
C = (Lv.^2)./(Ls.^2);
VK_FG = C.*(dLs.^2).*ifft2(QHFG); % S_K[f,g]
%% Visualization of results
% max
MxS = max(abs(SMNK_FG),[],'all');
MxV = max(abs(VK_FG),[],'all');
Mx = max([MxS,MxV,1]);
% min 
mxS = min(real(SMNK_FG),[],'all');
mxV = min(real(VK_FG),[],'all');
mx = min([mxS,mxV,0]);
%% Generating function values F,G on the visualization grid of size Lv*Lv 
[XLv,YLv] = meshgrid(xLv,xLv);
Fv = f(XLv,YLv);
Gv = g(XLv,YLv);
%--------------------------Surf plots--------------------------------------
figure 
tiledlayout(2,2);
%grid on;
% Visualization of f
h(1) = nexttile;
surf(XLv,YLv,Fv);
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$f(x,y)$','Interpreter','latex')
title('$f$','Interpreter','latex')
xlim([-a a]) 
ylim([-a a])
zlim([mx Mx])
daspect([1 1 1])
% Visualization of g
h(2) = nexttile;
surf(XLv,YLv,Gv);
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$g(x,y)$','Interpreter','latex')
title('$g$','Interpreter','latex')
xlim([-a a]) 
ylim([-a a])
zlim([mx Mx])
daspect([1 1 1])
% Visualization of SMNK[F,G]
h(3) = nexttile;
surf(XLv,YLv,real(SMNK_FG));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(S_{' num2str(M) ',' num2str(N) '}' ...
    '^{' num2str(a) ',' num2str(Ks) '}[f,g](x,y))$'], ...
    'Interpreter','latex')
title(['$f\ast g\approx \Re(S_{' num2str(M) ',' num2str(N) '}' ...
    '^{' num2str(a) ',' num2str(Ks) '}[f,g])$ ' ...
    'computed using iFFBT'],'Interpreter','latex')
xlim([-a a]) 
ylim([-a a])
zlim([mx Mx])
daspect([1 1 1])
% Visualization of VK[F,G]
h(4) = nexttile;
surf(XLv,YLv,real(VK_FG));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(V_{' num2str(Ks) '}^{' num2str(a) '}' ...
    '[f,g](x,y))$'],'Interpreter','latex')
title(['$f\ast g\approx \Re(V_{' num2str(Ks) '}' ...
    '^{' num2str(a) '}[f,g])$ computed using iFFT'],'Interpreter','latex')
xlim([-a a]) 
ylim([-a a])
zlim([mx Mx])
daspect([1 1 1])
% 
sgtitle(['The functions $f,g$, and approximations ' ...
    'of $f\ast g$ using iFFBT vs. iFFT'],'Interpreter','latex');
set(h,'Colormap',jet,'CLim',[mx Mx]);
cb = colorbar;
cb.Layout.Tile ='south';
