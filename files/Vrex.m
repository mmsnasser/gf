function Vz = Vrex(ak,r,z)
%
% Aug 14, 2024
% Vz = Vrex(alpha,zk,Rk,n,r,z)
%
%%
%
%
%
if( nargin == 2 ) 
    Vz   = (2/pi)*atan((r-ak)./(2*sqrt(r*ak)));
elseif( nargin == 3 ) 
    z=z(:).';
    wk  =  (r-ak)./(r+ak);
    wz  =  (r-z)./(r+z);
    w   =  i*sqrt(wz.^2-wk.^2);
    Vz  = (1./pi).*carg((w-wk)./(w+wk));
end
%%
end
