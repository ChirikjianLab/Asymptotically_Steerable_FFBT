% ACOM_fun2D_FFBTvsFFT_Test.m
% This code compares iFFBT and iFFT2 of a given complex function f vs. f
%%% Inputs:
% (1) given complex function f in the form of a function bundle OR 
% finite values on sampling/visualization grid. 
% Hint: If f is given as sampled values then the function handle should be 
% replace in the related steps, accordingly. 
% (2) the radius a which f is well supported in B_a
% (3) the Fourier-Bessel approximation order (M,N)
% (4) grid size for visualization Lv*Lv
%%% Outputs: 
% (1) Values of iFFBT of order K_MN on the visualization grid of size Lv*Lv
% (2) Values of iFFT2 of order K_MN on the visualization grid of size Lv*Lv
% (3) Surface plots comparison of f vs. iFFBT of f vs. iFFT2 of f.
%%=========================================================================
% Author: A-GF, 2023. Version 1
%%=========================================================================
% Hint: This test code uses the following zero solver in Line 32 and 46;
% Jason Nicholson, Bessel Zero Solver, MATLAB Central File Exchange.  
% (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
%%
close all
clear all
%% Initiating parameters and approximation order 
a = 1; 
b = 1./a;
M = 15;
N = 5;
%% Computing minimum approximation order of DFBT-FFT based on (M,N) 
Z = zeros(M+1,N);
for t = 1:M+1
     Z(t,:) = besselzero(t-1,N,1)./pi; % See Line 18.
end
CZ = ceil(Z);
K_MN = max(max(CZ));
K = K_MN; % K >= K_MN in S_MNK, according to Theorem 3.5 of the paper
Ls = 2.*K+1; % Finite Fourier sampling size
%% Generating visualization grid of size Lv*Lv
Lv = Ls+34; % Lv can be larger than Ls
dLv = 2.*a./Lv;
VLv = 1:Lv;
xLv = -a + (VLv-1).*dLv;
%% Generating Fourier-Bessel basis values on the visualization grid Lv*Lv
P = zeros(2*M+1,N,Lv,Lv);
for t = 1:2*M+1
    ZtN = besselzero(t-1-M,N,1); % See Line 18.
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
%% Generating H-matrix
H = zeros(Ls,Ls,2.*M+1,N);
for t = M+1:2*M+1
    for n = 1:N
        for l = 1:K+1
            for k = 1:K+1
            lk = [l-1,k-1];
            H(l,k,t,n) = (-1)^(l+k).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = K+2:Ls
            for k = 1:K+1
            lk = [l-Ls-1,k-1];
            H(l,k,t,n) = (-1)^(l+k-1).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = 1:K+1
            for k = K+2:Ls
            lk = [l-1,k-Ls-1];
            H(l,k,t,n) = (-1)^(l+k-1).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = K+2:Ls
            for k = K+2:Ls
            lk = [l-Ls-1,k-Ls-1];
            H(l,k,t,n) = (-1)^(l+k).*CoE(lk,t-1-M,n)./a;
            end
        end
    end
end
% negative m
for t = 1:M
    for n = 1:N
        for l = 1:K+1
            for k = 1:K+1
            lk = [l-1,k-1];
            H(l,k,t,n) = (-1)^(t-1-M).*(-1)^(l+k).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = K+2:Ls
            for k = 1:K+1
            lk = [l-Ls-1,k-1];
            H(l,k,t,n) = (-1)^(t-1-M).*(-1)^(l+k-1).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = 1:K+1
            for k = K+2:Ls
            lk = [l-1,k-Ls-1];
            H(l,k,t,n) = (-1)^(t-1-M).*(-1)^(l+k-1).*CoE(lk,t-1-M,n)./a;
            end
        end
        for l = K+2:Ls
            for k = K+2:Ls
            lk = [l-Ls-1,k-Ls-1];
            H(l,k,t,n) = (-1)^(t-1-M).*(-1)^(l+k).*CoE(lk,t-1-M,n)./a;
            end
        end
    end
end
%% Generating K-matrix
Kmat = zeros(Lv,Lv,Ls,Ls);
for i = 1:Lv
    for j = 1:Lv
        for l = 1:Ls
            for k = 1:Ls
                Hlk = squeeze(H(l,k,:,:));
                Pij = squeeze(P(:,:,i,j));
                PijT = transpose(Pij);
                Tijlk = trace(Hlk*PijT);
                Kmat(i,j,k,l) = Tijlk;
            end
        end
    end
end
%% ----------------- Generating Function Sample Values --------------------
% Hint: The best convergence of the iFFBT numerical algorithm is guaranteed 
% for smooth functions on O_a well-supported in B_a and vanishing on B_a.
%-- Example 4.1 (Fourier-Bessel basis functions)
U = @(x,y) Psi(1,2,a,x,y)+Psi(2,1,a,x,y);
%-- Example 4.2 (normalized approximately well-supported Gaussian)
% A1 = [0.1,0;0,0.05];
% A = inv(A1);
% B1 = [0.05,0;0,0.1];
% B = inv(B1);
% U = @(x,y) exp(-[x,y]*A*[x;y])+1i.*exp(-[x,y]*B*[x;y]);
%% Generating samples of U on the sampling grid of size Ls*Ls 
dLs = 2.*a./Ls;
VLs = 1:Ls;
xLs = -a + (VLs-1).*dLs;
ULs = zeros(Ls,Ls);
for i = 1:Ls
    for j = 1:Ls
        ULs(i,j) = U(xLs(i),xLs(j)); % if U is not a function handle then 
        % this line should be evaluated directly by the function rule
    end
end
FULs = fft2(ULs);
%% Generating SMNK on the visualization grid of size Lv*Lv
SMNK = zeros(Lv,Lv);
for i = 1:Lv
    for j = 1:Lv
        Kij = squeeze(Kmat(i,j,:,:));
        SMNK(i,j) = dLs.^2.*trace(Kij*FULs); 
    end
end
%% Generating Fourier series of order K_MN on Lv*Lv
EU = zeros(Lv,Lv);
for t = 1:K+1
    for s = 1:K+1
        EU(t,s) = FULs(t,s);
    end
end
%
for t = 1:K+1
    for s = Lv-K+1:Lv
        EU(t,s) = FULs(t,Ls+s-Lv);
    end
end
%
for t = Lv-K+1:Lv
    for s= 1:K+1
        EU(t,s) = FULs(Ls+t-Lv,s);
    end
end
%
for t= Lv-K+1:Lv
    for s = Lv-K+1:Lv
        EU(t,s) = FULs(Ls+t-Lv,Ls+s-Lv);
    end
end
% Generating S_K[U] on the visualization grid
VULv =  Lv.*Lv.*ifft2(EU)./(Ls.*Ls);
%% Generating U on the visualization grid of size Lv*Lv 
[XLv,YLv] = meshgrid(xLv,xLv);
ULv = zeros(Lv,Lv);
for i = 1:Lv
     for j = 1:Lv
         ULv(i,j) = U(XLv(i,j),YLv(i,j));% if U is not function handle then 
         % this line should be evaluated directly by the function rule
     end
 end
%% ================ Visualization of S_MNK vs. V_MNK vs. ULv ==============
%------------------------ Surface plots of real parts ---------------------
% Hint: For better visualization, color map configurations should be 
% adopted according to parameters of the given function.   
% Computing bounds for visualization of data 
MxReS = max(real(SMNK),[],'all');
MxReV = max(real(VULv),[],'all');
MxReU = max(real(ULv),[],'all');
MxRe = max([MxReS,MxReV,MxReU]);
% min
MnReS = min(real(SMNK),[],'all');
MnReV = min(real(VULv),[],'all');
MnReU = min(real(ULv),[],'all');
MnRe = min([MnReS,MnReV,MnReU]);
%----------
figure
tiledlayout(1,3);
figRe(1) = nexttile;
surf(XLv,YLv,real(ULv));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$\Re(f(x,y))$','Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
%
figRe(2) = nexttile;
surf(XLv,YLv,real(SMNK).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(S_{' num2str(M) ',' num2str(N) '}^' ...
    '{' num2str(a) ',' num2str(K) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
%
figRe(3) = nexttile;
surf(XLv,YLv,real(VULv).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(V_{' num2str(K) '}^{' num2str(a) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
view(3);
grid on
sgtitle(['The functions $\Re(f)$ and approximations using iFFBT vs. iFFT ' ...
    'on the uniform grid of size $ ' num2str(Lv) ' \times ' num2str(Lv) ' $'],'Interpreter','latex');
%set(h,'DataAspectRatio',[1 1 1])
set(figRe,'Colormap',parula,'CLim',[MnRe MxRe]);
cb = colorbar;
cb.Layout.Tile ='south';
%%%%------------------ Surface plots of imaginary parts ------------------- 
% Computing bounds for visualization of data 
% max Im
MxImS = max(imag(SMNK),[],'all');
MxImV = max(imag(VULv),[],'all');
MxImU = max(imag(ULv),[],'all');
MxIm = max([MxImS,MxImV,MxImU]);
% min Im
MnImS = min(imag(SMNK),[],'all');
MnImV = min(imag(VULv),[],'all');
MnImU = min(imag(ULv),[],'all');
MnIm = min([MnImS,MnImV,MnImU]);
%
figure
tiledlayout(1,3);
figIm(1) = nexttile;
surf(XLv,YLv,imag(ULv));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$\Im(f(x,y))$','Interpreter','latex')
daspect([1 1 1])
%
figIm(2) = nexttile;
surf(XLv,YLv,imag(SMNK).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Im(S_{' num2str(M) ',' num2str(N) '}^' ...
    '{' num2str(a) ',' num2str(K) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnIm MxIm])
daspect([1 1 1])
%
figIm(3) = nexttile;
surf(XLv,YLv,imag(VULv).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Im(V_{' num2str(K) '}^{' num2str(a) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnIm MxIm])
daspect([1 1 1])
view(3);
grid on
sgtitle(['The functions $\Im(f)$ and approximpations using iFFBT vs. iFFT ' ...
    'on the uniform grid of size $ ' num2str(Lv) ' \times ' num2str(Lv) ' $'],'Interpreter','latex');
%set(h,'DataAspectRatio',[1 1 1])
set(figIm,'Colormap',parula,'CLim',[MnIm MxIm]);
cb = colorbar;
cb.Layout.Tile ='south';

%------------------------ Surface plots in one figure --------------------- 
Mx = max(MxIm,MxRe);
Mn = min(MnIm,MnRe);

figure
tiledlayout(2,3);
fig(1) = nexttile;
surf(XLv,YLv,real(ULv));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$\Re(f(x,y))$','Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
%
fig(2) = nexttile;
surf(XLv,YLv,real(SMNK).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(S_{' num2str(M) ',' num2str(N) '}^' ...
    '{' num2str(a) ',' num2str(K) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
%
fig(3) = nexttile;
surf(XLv,YLv,real(VULv).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Re(V_{' num2str(K) '}^{' num2str(a) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnRe MxRe])
daspect([1 1 1])
%
fig(4) = nexttile;
surf(XLv,YLv,imag(ULv));
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$\Im(f(x,y))$','Interpreter','latex')
zlim([MnIm MxIm])
daspect([1 1 1])
%
fig(5) = nexttile;
surf(XLv,YLv,imag(SMNK).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Im(S_{' num2str(M) ',' num2str(N) '}^' ...
    '{' num2str(a) ',' num2str(K) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnIm MxIm])
daspect([1 1 1])
%
fig(6) = nexttile;
surf(XLv,YLv,imag(VULv).');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel(['$\Im(V_{' num2str(K) '}^{' num2str(a) '}(f)(x,y))$'],'Interpreter','latex')
zlim([MnIm MxIm])
daspect([1 1 1])
view(3);
grid on
sgtitle('The functions $f$ and approximations of $f$ using iFFBT vs. iFFT','Interpreter','latex');
%set(h,'DataAspectRatio',[1 1 1])
set(fig,'Colormap',parula,'CLim',[Mn Mx]);
cb = colorbar;
cb.Layout.Tile ='south';