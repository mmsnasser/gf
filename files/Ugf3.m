function V0 = Ugf3(alpha,cv,Rv,n,r,ind)
%
% May 23, 2024
%
%%
c = cv(ind); R = Rv(ind);
%
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%      
cr = real(exp(-i*angle(c))*c);
x = (r^2+cr^2-R^2)/(2*cr);
y =  sqrt(r^2-x^2);
%
vv = exp(i*angle(c))*[x-i*y   x+i*y];
cc = [c           0];
dd = [-1          1];
[eto,etop]=plgsegcirarcp(vv,cc,dd,n/2);
%
eti = []; etip = [];
for k=1:ind-1
    eti  = [eti;cv(k)+Rv(k)*exp(-i*t)];
    etip = [etip;  -i*Rv(k)*exp(-i*t)];
end
%
% Ueti = Ugf1(alpha,cv(ind),Rv(ind),n,r,eti(:).').';
Ueti = Urex(alpha,cv(ind),Rv(ind),n,r,eti(:).').';
%
et  = [eto;eti];
etp = [etop;etip];
gam0  = [zeros(n,1);1-Ueti];
A    =  et-alpha;
alphak = [cv(1:ind-1)]; alphak=alphak(:);
m=length(alphak);
%
[~,h0]=fbie(et,etp,A,gam0,n,5,[],1e-14,100);
for j=1:m+1
    H0(j,1)=mean(h0(1+(j-1)*n:j*n));
end
%
for k=1:m
    gamk{k} = log(abs(et-alphak(k)));
    [~,hk{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,100);
    for j=1:m+1
        Hm(j,k)=mean(hk{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat=[Hm ones(m+1,1)]; 
rhs  = -H0;
x=mat\rhs; a=x(1:m,1); cm = x(m+1);
% 
V0 = cm-sum(a.*log(abs(alphak-alpha)));
%
end