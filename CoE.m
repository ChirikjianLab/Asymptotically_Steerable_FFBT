% This function computes the value c(k;m,n) according to (2.14) for given 
% integer vector k=(k1,k2), non-negative integer m, and positive integer n.
%%=========================================================================
% Author: A-GF, 2023. Version 1.
%%=========================================================================
% Hint: This test code uses the following zero solver in Line 11;
% Jason Nicholson, Bessel Zero Solver, MATLAB Central File Exchange.  
% (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
%--------------------------------------------------------------------------
function [c] = CoE(k,m,n)
ZmnV = besselzero(m,n,1); % See Line 6.
zmn = ZmnV(n);
[Pk,Rk] = cart2pol(k(1),k(2));
bmnk = (pi.*Rk).^2-zmn.^2;
Bmnk = 2.*bmnk;
amnk = besselj(m,pi.*Rk,1).*exp(-1i.*m.*Pk);
c = sqrt(pi).*(-1).^(n).*1i.^(m).*zmn.*amnk./Bmnk;
end