clear all
close all
clc

Result = [];
rmax=0.5;
a = 0.01;
b = 1;
w=3;
for r = rmax/10:0.001:rmax/2;
    dC = DerivativeC(r,rmax,a,w,b);
    TC = getTotalCost(r,rmax,a,w,b);
    Result =[Result;r,dC,TC];
end

plot(Result(:,1),Result(:,3),'-k','LineWidth',2)
ylabel('\itC','FontSize',30);
xlabel('\itr','FontSize',30);
set(gca,'YTick',[1.06:0.04:1.2])
set(gca,'YTickLabel',[{'1.06'},{'1.10'},{'1.14'},{'1.18'}])
set(gca,'FontSize',25);
set(gcf,'paperposition',[0 0 8 8]);
saveas(1,'FigS2','png');
