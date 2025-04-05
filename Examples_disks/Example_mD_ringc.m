clear
clc
%
addpath ../bie; addpath ../fmm; addpath ../files;
%%
%
n       =  3*5*2^8;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
% 
m  =  10;
cv =  5.*exp(2i*pi*[0:m-1]'./m);
R  =  1;
Rv = R+zeros(size(cv));
%
rv =   [0,1,2,3,3.9999,4.02,4.05:0.1:11.95,12].';
alpha   =  0;
z = alpha;
%%
for kk=1:length(rv)
    r = rv(kk);
    r
    %
    if r<cv(1)-R
        gr(kk,1) = 0;
    elseif r>cv(1)-R & r<cv(1)+R
        for k=1:m
            cc = real(exp(-i*angle(cv(k)))*cv(k));
            x = (r^2+cc^2-R^2)/(2*cc);
            y =  sqrt(r^2-x^2);
            vv = exp(i*angle(cv(k)))*[x-i*y   x+i*y];
            a(k) = vv(1); b(k) = vv(2);
        end
        %
        for k=1:m
            vert(2*k-1)=  a(k);
            vert(2*k)  =  b(k);
            cent(2*k-1) =  cv(k);
            dr(2*k-1)   = -1;
            cent(2*k)   =  0;
            dr(2*k)     =  1;
        end
        [et,etp]=plgsegcirarcp(vert,cent,dr,n/(2*m));
        %        
        if r==5
            figure
            plot(real(et),imag(et))
            axis equal
        end
        %
        A     =  et-alpha;  
        gam   = -log(abs(et-alpha));
        [mu,h]=  fbie(et,etp,A,gam,n,5,[],0.5e-14,100);
        fet   = (gam+h+i*mu)./A;   cm = exp(-mean(h));
        %%
        f_z   =  fcau(et,etp,fet,z(:).');
        Phiofz = cm.*(z-alpha).*exp((z-alpha).*f_z); 
        % Compute the arc \hat L
        zet   = cm.*(et-alpha).*exp((et-alpha).*fet);
        if r==5
            figure
            plot(real(zet),imag(zet))
            axis equal
        end
        %
        ai = zet(1); a0 = zet(n/2+1); a1 = zet(3*n/4+1); cs = (a1-ai)/(a1-a0);
        %
        Mobf = @(z)cs*(1+(ai-a0)./(z-ai));
        %
        for j=1:2*m-1
            zr(j,1) = zet(j*n/(2*m)+1);
        end
        %
        xr = Mobf(zr);
        MPhiofz = Mobf(Phiofz); 
        %
        for j=1:m
            kr(2*j-1)= 1;
            kr(2*j)  = 0;
        end
        U0=0;
        for j=1:2*m-1
            U0=U0+(((-1)^(j+1))/pi)*angle(MPhiofz-xr(j));
        end
        gr(kk,1) = U0;
        %
    elseif r>cv(1)+R
        gr(kk,1) = Ugf2(alpha,cv,Rv,n,r,z);
    end
end
%%
[rv gr]
%
figure;
hold on; box on
plot(rv,gr,'-k','LineWidth',1.5);
axis square
grid on
axis([0 12 0 1])
%%
m = length(cv);
eto     =   max(rv).*exp(i.*t);
etc=eto; 
for j=1:m
    etc  =[etc;  cv(j)+Rv(j).*exp(-i.*t)];
end
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=2:m+1
    crv = etc(1+(k-1)*n:k*n); crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5)
end
plot(real(alpha),imag(alpha),'ok','MarkerFaceColor','k','MarkerSize',4)
set(gca,'FontSize',14)
axis([-7 7 -7 7])
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
%%