% This code creates the function handle to compute the value of the 
% B_a normalized (m,n)-Fourier-Bessel basis function $\Psi_{m,n}^a(x,y)$ 
% according to at given Cartesian coordinates (x,y). 
%%=========================================================================
% Author: A-GF, 2023. Version 1.
%%=========================================================================
% Hint: This test code uses the following zero solver in Line 13;
% Jason Nicholson, Bessel Zero Solver, MATLAB Central File Exchange.  
% (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
%--------------------------------------------------------------------------
function [P] = Psi(m,n,a,x,y)
[t,r] = cart2pol(x,y);
Zmn = besselzero(m,n,1); % See Line 7.
zmn = Zmn(n);
Dmna = a.*sqrt(pi).*abs(besselj(m+1,zmn));
Nmna = 1./Dmna;
b = 1./a;
Jmnar = besselj(m,b.*zmn.*r);
P = Nmna.*exp(1i.*m.*t).*Jmnar.*rectpuls(b.*r./2);
end