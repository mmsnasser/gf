function map = PreImageStrSlit(Lc,Lk,thetk,r,n,tol,Maxiter)
%
%
%
%%
t         =   (0:2*pi/n:2*pi-2*pi/n).'; 
m         =   length(thetk);
for k=1:m
    thet(1+(k-1)*n:k*n,1)  =  thetk(k);
end
cent     =  Lc;
radx     = (1-0.5*r).*Lk;
%%

%%
err = inf;
itr = 0;
while (err>tol)
    itr  =itr+1;  
%%
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*radx(k).*( cos(t)-i*r.*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*radx(k).*(-sin(t)-i*r.*cos(t));    
end
%%


%%
A       =  exp(i.*(pi/2-thet));
gam     =  imag(exp(-i.*thet).*et);
%%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-13,200);
%%
fnet            = (gam+h+i.*mun)./A;
zet             =  et+fnet;
rotzet          =  exp(-i.*thet).*zet;
for k=1:m
    wnL         =  rotzet((k-1)*n+1:k*n,1);
    centk(k,1)  =  exp(i.*thetk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
    radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
end
cent  =  cent -1.0.*(centk-Lc);
radx  =  radx -(1-0.5*r).*(radk-Lk) ;
err   = (norm(centk-Lc,1)+norm(radk-Lk,1))/m;
% 'iteration no.   error'
[itr err]
error (itr,1) = err;
itrk  (itr,1) = itr;
%%
if itr>=Maxiter
    'No convergence after Maximunm number of iterations'
    break;
end
end
%%
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*radx(k).*( cos(t)-i*r.*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*radx(k).*(-sin(t)-i*r.*cos(t));    
end
A        =  exp(i.*(pi/2-thet));
gam      =  imag(exp(-i.*thet).*et);
[mun,h]  =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
fnet     = (gam+h+i.*mun)./A;
zet      =  et+fnet;
%
for k=1:m
    Jk=(k-1)*n+1:k*n;
    zetp(Jk,1) = derfft(real(zet(Jk)))+i*derfft(imag(zet(Jk)));
end
%%
map.zet   =  zet;
map.zetp  =  zetp;
map.et    =  et;
map.etp   =  etp;
map.cent  =  cent;
map.radx  =  radx;
map.fnet  =  fnet;
map.r     =  r;
%%
end