clear
clc
%
addpath ../bie; addpath ../fmm; addpath ../files;
%%
n =  2^12;
%
t       =  (0:2*pi/n:2*pi-2*pi/n).';
% 
cv =  [1;2;3;4;5].*exp(2i*pi*[0:4].'/5);
Rv =  [0.4;0.4;0.4;0.4;0.4];
%
rv=[];
rv=[0;abs(cv(1))-Rv(1)-0.0001];
for k=1:length(cv)-1
    rv=[rv;linspace(abs(cv(k))-Rv(k)+0.02,abs(cv(k))+Rv(k)-0.02,10)'];
    rv=[rv;linspace(abs(cv(k))+Rv(k)+0.02,abs(cv(k+1))-Rv(k+1)-0.02,10)'];
end
rv=[rv;linspace(abs(cv(end))-Rv(end)+0.01,abs(cv(end))+Rv(end)-0.01,10)'];
rv=[rv;linspace(abs(cv(end))+Rv(end)+0.01,8,10)'];
% rv =[-0.025:0.05:8.025].';
alpha   =  0;
z = alpha;
%%
for kk=1:length(rv)
    r = rv(kk);
    r
    %
    if r<abs(cv(1))-Rv(1)
        r
        gr(kk,1) = 0;
    elseif r>abs(cv(end))+Rv(end)
        gr(kk,1) = Ugf2(alpha,[cv(1:end)],[Rv(1:end)],n,r,z);
    elseif r>abs(cv(1))-Rv(1) & r<abs(cv(end))+Rv(end)
        for jj=1:length(cv)
            if r>abs(cv(jj))-Rv(jj) & r<abs(cv(jj))+Rv(jj)
                U0   = Urex(alpha,cv(jj),Rv(jj),n,r);
                if jj>1
                    V0   = Ugf3(alpha,cv,Rv,n,r,jj);
                else 
                    V0=0;
                end
                gr(kk,1) = U0+V0;
            end
        end
        for jj=1:length(cv)-1
            if r>abs(cv(jj))+Rv(jj) & r<abs(cv(jj+1))-Rv(jj+1)
                gr(kk,1) = Ugf2(alpha,[cv(1:jj)],[Rv(1:jj)],n,r,z);
            end
        end
    end
end
%%
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(rv,gr,'-k','LineWidth',1.5);
xlabel('$r$','interpreter','latex')
ylabel('$g(r)$','interpreter','latex')
set(gca,'FontSize',14)
axis([0 8 0 1])
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc Fig5Dsp
%%
m = length(cv);
eto     =   max(rv).*exp(i.*t);
et=eto; 
for j=1:m
    et  =[et;  cv(j)+Rv(j).*exp(-i.*t)];
end
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(real(alpha),imag(alpha),'ok','MarkerFaceColor','k','MarkerSize',4)
k=1;crv = et(1+(k-1)*n:k*n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),':k','LineWidth',1.5)
for k=2:m+1
    crv = et(1+(k-1)*n:k*n); crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5)
end
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal
%%