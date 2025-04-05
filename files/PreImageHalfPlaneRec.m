function map = PreImageHalfPlaneRec (Lc,Lk,thetk,ratio,n,tloerance,Maxiter,w0)
% PreImageHalfPlaneRec.m
% Nasser, Aug 14, 2024
% Let Omega be the upper half-plane slitted along m segments slits where
% m=length(thetk)-1 (i.e., Omega is an unbounded multiply connected domain 
% of connectivity m+1)
% Lc: is a vector of the center points of these slits
% Lc: is a vector of the length of these slits
% thetk: is a vector of the angles between these slits and the positive
% real line.
% 
% This function computes: et , etp , cent , wet; i.e., a preimage 
% bounded multiply connected domain G border by the unit disk and
% m pseudellipses interior to the unit disk:
% et: the parametruzation of the boundary of G. 
% etp: the derivative of et.
% cent: the center of the unit disk and the m-1 pseudellipses
% wet: the boundary values of the conformal mapping from G onto Omega,
% 
% ratio, 0<ratio<=1 is the ration of lengthg of minor axix to the major 
% axis of the ellipses. It is better to choose the ratio close to 1. 
% However, when the two slits are close to each other, we need to chose 
% the ratio small but not very small. May be >=0.01. But for small value 
% we need to increase the value of n.
% 
% n: the number of discretization points in each boundary component of G
% tolerance: is the tolerance of the iterative method
% Maxiter: is the maximum number of iterations for the iterative method
%
%
% This code of computing the preimage domain G and the conformal mapping
% is based on the iterative method presented in Section 4 in the paper:
% M. Nasser and C. Green, A fast numerical method for ideal fluid flow 
% in domains with multiple stirrers, Nonlinearity 31 (2018) 815-837.
%
%
% OOOOOLD
Lc = [0;Lc];
Lk = [0;Lk];
thetk = [0;thetk];
%
% image of alpha=0
alpha   =  0;
w0i     =  w0/i;
%
t             =   (0:2*pi/n:2*pi-2*pi/n).'; 
%
Psi       =  @(z)(1i.*(2i./(i-z)-1)); %i(i+z)/(i-z)
Psip      =  @(z)((-2./((i-z).^2)));
Psiv      =  @(z)(i+2./(i+z));
Psivp     =  @(z)(-2./((i+z).^2));
%
m             =   length(thetk)-1
thet(1:n,1)   =   0;
for k=2:m+1
    thet(1+(k-1)*n:k*n,1)  =  thetk(k);
end
cent     =  Lc;
radx     = (1-0.5*ratio).*Lk;
rady     =  ratio.*radx;
%
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
err = inf;
itr = 0;
while (err>tloerance)
    itr  =itr+1;  
    %
    for k=2:m+1
        xt(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
        xtp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
    end
    %
    for k=2:m+1
        et(1+(k-1)*n:k*n,1)    =  Psiv(xt(1+(k-1)*n:k*n,1));
        etp(1+(k-1)*n:k*n,1)   =  Psivp(xt(1+(k-1)*n:k*n,1)).*xtp(1+(k-1)*n:k*n,1);
    end
    %    
    A       =  exp(i.*(pi/2-thet)).*(et-alpha);
    gam(1:n,1)   =   0;
    for k=2:m+1
        gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*Psi(et((k-1)*n+1:k*n,1))); 
    end
    %
    [mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
    h0            =  sum(h(1:n,1))/n;
    %
    fnet          = (gam+h+i.*mun)./A;
    wet           =  w0i*(Psi(et)+(et-alpha).*fnet+i*h0)/(1+h0);
    rotwn         =  exp(-i.*thet).*wet;
    for k=2:m+1
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*thetk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    cent   =  cent-1.0.*(centk-Lc);
    radx   =  radx-(1-0.5*ratio).*(radk-Lk) ;
    rady   =  ratio.*radx;
    err    = (norm(centk-Lc,1)+norm(radk-Lk,1))/(m);
    [itr err]
    error (itr,1) = err;
    itrk  (itr,1) = itr;
    %
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%
for k=2:m+1
    xt(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    xtp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  Psiv(xt(1+(k-1)*n:k*n,1));
    etp(1+(k-1)*n:k*n,1)   =  Psivp(xt(1+(k-1)*n:k*n,1)).*xtp(1+(k-1)*n:k*n,1);
end
A       =  exp(i.*(pi/2-thet)).*(et-alpha);
gam(1:n,1)   =   0;
for k=2:m+1
    gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*Psi(et((k-1)*n+1:k*n,1))); 
end
[mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
h0            =  sum(h(1:n,1))/n;
fnet          = (gam+h+i.*mun)./A;
zet           =  w0i*(Psi(et)+(et-alpha).*fnet+i*h0)/(1+h0);
%
%%
map.et   =  et;
map.etp  =  etp;
map.cent =  Psiv(cent(2:end));
map.zet  =  zet;
%%
end