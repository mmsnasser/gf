function [cap,a,uz] = capm(et,etp,alphav,deltv,m,alpha,zv)
% Compute the capacity of the condensers (G,F)
% 
% Input:
% 1,2) et, etp: parametrization of the boundary and its first derivative
% 3) alphav=[alphav(1),...,alphav(mp)]: alphav(j) is an auxiliary point
% interior to \Gamma_j
% 4) m: the number of the closed sets E_k
% 8) alpha: for bounded G, alpha is an auxiliary point in G
% 
% Output:
% cap (the capacity of the condensers (G,F)).
%
% Computing the constants \h_{j,k} for j=1,2,...,m
%
n=length(et)/(m+1); 
A=et-alpha;
for k=1:m
    gamk{k}=log(abs(et-alphav(k)));
    [mu{k},h{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,100);
    for j=1:m+1
        hjk(j,k)=mean(h{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat=hjk; mat(1:m+1,m+1)=1; 
rhs  = deltv(:); delt = deltv(2:end,1);
x=mat\rhs; a=x(1:m,1); c = x(m+1);
% Computing the capacity
cap  = (2*pi)*sum(delt.*a);
%
uz   = [];
%
if( nargin == 7 ) 
    Afet = zeros(size(et));
    for k=1:m
        Afet = Afet + a(k)*(gamk{k}+h{k}+i*mu{k});
    end
    fet = Afet./A;
    fzv = fcau(et,etp,fet,zv);
    uz  = c+real((zv-alpha).*fzv);
    for k=1:m
        uz = uz-a(k)*log(abs(zv-alphav(k)));
    end
end
%
end