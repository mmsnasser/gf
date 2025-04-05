clear
clc
%
addpath ../bie; addpath ../fmm; addpath ../files;
%%
n =  2^12;
%
t       =  (0:2*pi/n:2*pi-2*pi/n).';
% 
zk =  [5];
Rk =  [1];
%
rv =   [0,1,2,3,3.9999,4.02,4.05:0.1:11.95,12].';
alpha   =  0;
z = alpha;
%%
for kk=1:length(rv)
    r = rv(kk);
    r
    %
    if r<abs(zk(1))-Rk(1)
        r
        gr(kk,1) = 0;
    elseif r>abs(zk(end))+Rk(end)
        gr(kk,1) = Ugf2(alpha,[zk(1:end)],[Rk(1:end)],n,r,z);
    elseif r>abs(zk(1))-Rk(1) & r<abs(zk(end))+Rk(end)
        for jj=1:length(zk)
            if r>abs(zk(jj))-Rk(jj) & r<abs(zk(jj))+Rk(jj)
                U0   = Urex(alpha,zk(jj),Rk(jj),n,r);              
                if jj>1
                    V0   = Ugf3(alpha,zk,Rk,n,r,jj);
                else 
                    V0=0;
                end
                gr(kk,1) = U0+V0;
            end
        end
        for jj=1:length(zk)-1
            if r>abs(zk(jj))+Rk(jj) & r<abs(zk(jj+1))-Rk(jj+1)
                gr(kk,1) = Ugf2(alpha,[zk(1:jj)],[Rk(1:jj)],n,r,z);
            end
        end
    end
end
%%
for kk=1:length(rv)
    r = rv(kk);
    %
    if r<=abs(zk(1))-Rk(1)
        gre(kk,1) = 0;
    elseif r>abs(zk(end))+Rk(end)
        gre(kk,1) = NaN;
    elseif r>abs(zk(1))-Rk(1) & r<abs(zk(end))+Rk(end)
        x = (r^2+abs(zk)^2-Rk^2)/(2*abs(zk));
        y =  sqrt(r^2-x^2);
        %
        xi1 = x-i*y;
        xi2 = x+i*y;
        %
        alpha = atan((r-x)/y);
        nu    = atan((y^2-(r-x)*(x+Rk-abs(zk)))/(y*(r+Rk-abs(zk))));
        %
        gre(kk,1) = atan((y^2-x*(r-x))/(r*y))/atan((y^2-(r-x)*(x+Rk-abs(zk)))/(y*(r+Rk-abs(zk))));
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
plot(rv,gre,'or','LineWidth',1.5);
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
%%





