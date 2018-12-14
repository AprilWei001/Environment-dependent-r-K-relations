clear all
close all
clc
A=importdata('R_K_strain_indexed.txt');
Data=A.data;
All_env_corr=[];
Q=[];
max_Col=[];
Proportion=[];
figure
for i=1:2:18;
    data=Data(:,i:i+1);
    ind=find(data(:,1)==-1);
    data(ind,:)=[];
    subplot(3,3,i/2+1/2);
    [data2 ind]=sort(data(:,1));
    mi=min(data(:,1));
    ma=max(data(:,1));
    Collect=[];
    for j=500:500:size(data,1);
        Collect=[Collect;mean(data(ind(j-199:j),:))];
    end

%     for j=1:100;
%         ind1=find(data(:,1)<j*0.004);
%         if length(ind1)>0;
%         ind2=find(data(ind1,1)>(j-1)*0.004);
%         if length(ind2)>200;
%             Collect=[Collect;mean(data(ind1(ind2),1)),mean(data(ind1(ind2),2))];
%         end
%         end
%     end
    plot(Collect(:,1),Collect(:,2),'ok','MarkerFaceColor','k','MarkerSize',15);
    [m indm]=max(Collect(:,2));
    max_Col=[max_Col,Collect(indm,1)];
    [rho p]=corr(data(:,1),data(:,2),'type','Spearman');
    hold on
    plot([0.1076 0.1076]',[0,max(Collect(:,2))+10^6]','.-k','LineWidth',2);
    ylim([min(Collect(:,2)),max(Collect(:,2))]);
    xlim([0.04 0.18])
    set(gca,'FontSize',25);
    xlabel('{\itr} (bin)');
    ylabel('{\itK} (bin)');
    All_env_corr=[All_env_corr;rho,p];
    Q=[Q;mean(data)];
    %text(Q(i/2+1/2,1),Q(i/2+1/2,2),strcat('rho=',num2str(rho)),'FontSize',25);
    indp=find(data(:,1)<0.1074);
    Proportion=[Proportion;length(indp)/length(data(:,1))];
end

Env=[ {'hydroxyurea'},{'NaCl'},{'allantoin'},{'caffeine'},{'galactose'},{'glycine'},{'isoleucine'},{'phleomycin'},{'rapamycin'}];
mean(max_Col)
mean(max_Col)-sqrt(var(max_Col))
mean(max_Col)+sqrt(var(max_Col))
set(gcf,'paperposition',[0 0 20 20]);
subplot(3,3,1)
ylim([1.2*10^6 2.4*10^6])
set(gca,'YTick',[1.4:0.4:2.2]*10^6)
text(0.002,2.4*10^6,'a','FontSize',50)
title(Env(1));
subplot(3,3,2)
ylim([1.4*10^6 2.0*10^6])
title(Env(2));
set(gca,'YTick',[1.4:0.3:1.9]*10^6)
text(0.002,2.0*10^6,'b','FontSize',50)
subplot(3,3,3)
ylim([1.0*10^6 2.0*10^6])
set(gca,'YTick',[1.2:0.2:1.8]*10^6)
text(0.002,2.0*10^6,'c','FontSize',50)
title(Env(3));
subplot(3,3,4)
text(0.002,1.2*10^6,'d','FontSize',50)
ylim([0.7*10^6 1.2*10^6])
title(Env(4));
subplot(3,3,5)
ylim([2*10^6 5*10^6])
set(gca,'YTick',[2:1:5]*10^6)
title(Env(5));
text(0.002,5*10^6,'e','FontSize',50)
subplot(3,3,6)
ylim([0.8*10^6 1.4*10^6])
title(Env(6));
set(gca,'YTick',[0.9:0.2:1.4]*10^6)
text(0.002,1.45*10^6,'f','FontSize',50)
ylim([.85*10^6 1.45*10^6])
subplot(3,3,7)
text(0.002,1.88*10^6,'g','FontSize',50)
ylim([1.*10^6 1.88*10^6])
title(Env(7));
subplot(3,3,8)
ylim([2.8*10^6 4.8*10^6])
title(Env(8));
set(gca,'YTick',[2.8:1:4.8]*10^6)
text(0.002,4.8*10^6,'h','FontSize',50)
subplot(3,3,9)
ylim([1*10^6 1.65*10^6])
title(Env(9));
set(gca,'YTick',[1.1:0.2:1.6]*10^6)
text(0.002,1.65*10^6,'i','FontSize',50)

saveas(1,'Fig5','png');
