function Uz = Ugfs2(av,bv,r,n,alphao,zo)
%
% June 4, 2024
%
%%
zo=zo(:).';
%
t = (0:2*pi/n:2*pi-2*pi/n).';
[s,sp] =  deltw(t,2,2);
%
etc = r.*exp(i.*t);
Lc = (av+bv)/2; Lk = bv-av; thetk = zeros(size(Lc));
%
map = PreImageStrSlit(Lc,Lk,thetk,1,n,1e-12,100);
%
zeti    =   map.zet;
zetip   =   map.zetp;
eti     =   map.et;
etip    =   map.etp;
alphav  =   map.cent;
%    
eto     =  (etc.'+fcau(zeti,zetip,eti-zeti,etc.',n,0)).';
etop    =   derfft(real(eto))+i*derfft(imag(eto));
% 
et      =  [eto ;eti ];
etp     =  [etop;etip];
%
m     =  length(av);
deltv =  [0   ; ones(m,1)];
%
alpha =  (alphao+fcau(zeti,zetip,eti-zeti,alphao,n,0));
z     =  (zo+fcau(zeti,zetip,eti-zeti,zo,n,0));
%
[~,~,Uz] = capm(et,etp,alphav,deltv,m,alpha,z);
%%
end