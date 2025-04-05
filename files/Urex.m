function Uz = Urex(alpha,zk,Rk,n,r,z)
%
% Aug 1, 2024
% [U0 , Uz ] = Urex(alpha,zk,Rk,n,r,z)
%
%%
%
x = (r^2+abs(zk)^2-Rk^2)/(2*abs(zk));
y =  sqrt(r^2-x^2);
d =  r+Rk-abs(zk);
%
%
if( nargin == 5 ) 
    Uz   = atan(sqrt((r-x)./(r+x)))./atan(((2*r-d)./d).*sqrt((r-x)./(r+x)));
elseif( nargin == 6 ) 
    z=z(:).';
    xi2 = exp(i*carg(zk))*(x+i*y);
    xi1 = exp(i*carg(zk))*(x-i*y);
    alp = 2*atan(sqrt((r-x)./(r+x)));
    nu  = 2*atan(((2*r-d)./d).*sqrt((r-x)./(r+x)));
    Uz  = (1./nu).*carg(exp(-i*alp).*(z-xi2)./(z-xi1));
end
%%
end