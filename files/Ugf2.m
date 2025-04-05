function U0 = Ugf2(alpha,cv,Rv,n,r,z)
%
% May 23, 2024
%
%%
z=z(:).';
%
t       =  (0:2*pi/n:2*pi-2*pi/n).';
% 
m     =  length(cv);
%
eto     =   r.*exp(i.*t);
etop    =   r*i.* exp(i.*t);
%
deltv = [0]; et=eto; etp=etop;
for j=1:m
    deltv =[deltv;1];
    et  =[et;  cv(j)+Rv(j).*exp(-i.*t)];
    etp =[etp;        -i*Rv(j).*exp(-i.*t)];
end
%
[~,~,U0] = capm(et,etp,cv,deltv,m,alpha,z);
%%
end