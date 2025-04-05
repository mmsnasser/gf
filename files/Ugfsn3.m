function V0 = Ugfs3(av,bv,r,n,ind)
%
% June 5, 2024
%
%%
a    =  av(ind);
wk   =  (r-a)./(r+a);
psi  =  @(z) i*sqrt(((r-z)./(r+z)).^2-(wk).^2);
psiv =  @(w) r*((1-sqrt(wk.^2-w.^2))./(1+sqrt(wk.^2-w.^2)));
%
%%
for k=1:ind-1
    avn(k,1) = psi(bv(k));
    bvn(k,1) = psi(av(k));
end
Lc = (avn+bvn)/2; Lk = abs(bvn-avn); thetk = zeros(size(Lc))+pi/2;
%
map = PreImageHalfPlaneRec (Lc,Lk,thetk,1,n,1e-12,100,i*sqrt(1-wk^2));
%
et      =   map.et;
etp     =   map.etp;
alphav  =   map.cent;
%    

zet     =  map.zet;
zeth    =  zet(n+1:end);
zetd    =  psiv(zeth);
%%
Vzetd   =  Vrex(a,r,zetd(:).').';
%%
gam0  = [zeros(n,1);1-Vzetd];
alpha = 0;
A     =  et-alpha;
m     =  length(alphav)
%
[~,h0]=fbie(et,etp,A,gam0,n,5,[],1e-13,100);
for j=1:m+1
    H0(j,1)=mean(h0(1+(j-1)*n:j*n));
end
%
for k=1:m
    gamk{k} = log(abs(et-alphav(k)));
    [~,hk{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-13,100);
    for j=1:m+1
        Hm(j,k)=mean(hk{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat=[Hm ones(m+1,1)]; 
rhs  = -H0;
x=mat\rhs; a=x(1:m,1); cm = x(m+1);
% 
V0 = cm-sum(a.*log(abs(alphav)));
%%
end