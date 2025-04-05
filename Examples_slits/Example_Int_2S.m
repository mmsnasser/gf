clear
clc
%
addpath ../bie; addpath ../fmm; addpath ../files;
%%
n   =  2^13;
t = (0:2*pi/n:2*pi-2*pi/n).';
%
av = [1;3];
bv = [2;4];
%
rv=[];
rv=[0;av(1)-0.0001];
for k=1:length(av)-1
    rv=[rv;linspace(av(k)+0.01,bv(k)-0.001,40)'];
    rv=[rv;linspace(bv(k)+0.01,av(k+1)-0.001,40)'];
end
rv=[rv;linspace(av(end)+0.01,bv(end)-0.001,40)'];
rv=[rv;linspace(bv(end)+0.01,10,40)'];
%%
for kk=1:length(rv)
    r = rv(kk);
    r
    %
    if r<av(1)
        r
        gr(kk,1) = 0;
    elseif r>bv(end)
        alphao  =  0;  zo = 0;
        gr(kk,1) = Ugfsn2([av(1:end)],[bv(1:end)],r,n,alphao,zo);
    elseif r>av(1) & r<bv(end)
        for jj=1:length(av)
            if r>av(jj) & r<bv(jj)
                U0      =  Vrex(av(jj),r);
                if jj>1
                    V0   = Ugfsn3(av,bv,r,n,jj);
                else 
                    V0=0;
                end
                gr(kk,1) = U0+V0;
            end
        end
        for jj=1:length(av)-1
            if r>bv(jj) & r<av(jj+1)
                alphao  =  0;  zo = 0;
                gr(kk,1) = Ugfsn2([av(1:jj)],[bv(1:jj)],r,n,alphao,zo);
            end
        end
    end
end
%%
[rv  gr]
gre = (2/pi)*atan((rv-av(1))./(2*sqrt(rv*av(1))));
gre(rv<av(1))=0;
gre(rv>bv(1))=NaN;
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(rv,gr,'-k','LineWidth',1.5);
plot(rv,gre,':r','LineWidth',1.5);
% legend({'Approximate values','Exact values'},'Location','northwest')
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
