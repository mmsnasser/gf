clear
clc
%
addpath ../bie; addpath ../fmm; addpath ../files;
%%
n   =  2^12;
%
a  =  1; b = 8; 
%
rv =[0.0001,0.9999,1.01,1.05:0.1:7.95,7.99,8.01,8.05:0.1:9.95,10].';
%%
t = (0:2*pi/n:2*pi-2*pi/n).';
[s,sp] =  deltw(t,2,2);
%
%%
for kk=1:length(rv)
r = rv(kk)
%%
if r<a
    gr(kk,1)=0;
elseif r>a & r<b
    %
    gr(kk,1) = Vrex(a,r);
elseif r>b
    %
    alphao  =  0;  zo = 0;
    Uz = Ugfsn2(a,b,r,n,alphao,zo);
    gr(kk,1) = Uz;
end
end
%%
gre = (2/pi)*atan((rv-a)./(2*sqrt(rv*a)));
gre(rv<a)=0;
gre(rv>b)=NaN;
%%
[rv gr gre]
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(rv,gr,'-k','LineWidth',1.5);
plot(rv,gre,':r','LineWidth',1.5);
legend({'Approximate values','Exact values'},'Location','northwest')
xlabel('$r$','interpreter','latex')
ylabel('$g(r)$','interpreter','latex')
set(gca,'FontSize',14)
axis([0 12 0 1])
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc Fig8DRing
%%
