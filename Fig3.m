clear all
close all
clc

[QTL36 txt]=xlsread('QTL36_fitting_result.xlsx');
QTL36=QTL36(:,1);

[simQTL36 txt]=xlsread('simulatedQTL_fitting_result_36.xlsx');

simQTL36=simQTL36(:,1);
simQTL36=reshape(simQTL36,[18,100]);
Sim_mean=[];
Sim_SD=[];

for i=1:18;
    Sim_mean=[Sim_mean;mean(simQTL36(i,:))];
    Sim_SD=[Sim_SD;std(simQTL36(i,:))];
end

subplot('position',[0.1 0.25 0.22 0.63])
hold on

plot([1:9],QTL36(1:2:18),'.r','MarkerSize',30)%rQTL explain r
plot([1:9],QTL36(20:2:36),'sb','MarkerFaceColor','b','MarkerSize',10)%KQTL explain r
plot([1:9],Sim_mean(1:2:18),'vk','MarkerFaceColor','k','MarkerSize',10)
h=legend('{\itr}QTL','{\itK}QTL','random SNPs','location','SouthWest')
set(h,'FontSize',11.5)
xlim([0 10])
ylim([0 0.9])
set(gca,'FontSize',20);
box on
text(-3,0.99,'a','FontSize',50)
set(gca,'YTick',[0:0.3:0.9]);
set(gca,'XTick',[1:9]);
set(gca,'XTickLabel',[ {'hydroxyurea'},{'NaCl'},{'allantoin'},{'caffeine'},{'galactose'},{'glycine'},{'isoleucine'},{'phleomycin'},{'rapamycin'}]);
set(gca,'XTickLabelRotation',45)
set(gca,'YTickLabel',[{'0.0'},{'0.3'},{'0.6'},{'0.9'}]);

ylabel([{'Fraction of {\itr} variance'}, {'explained by 36 SNPs'}])
subplot('position',[0.41 0.25 0.22 0.63])
hold on
plot([1:9],QTL36(19:2:36),'.r','MarkerSize',30)%rQTL explain K
plot([1:9],QTL36(2:2:18),'sb','MarkerFaceColor','b','MarkerSize',10)%KQTL explain K
plot([1:9],Sim_mean(2:2:18),'vk','MarkerFaceColor','k','MarkerSize',10)
%legend('rQTL','KQTL','location','SouthEast')

xlim([0 10])
ylim([0 0.9])
set(gca,'YTick',[0:0.3:0.9]);
set(gca,'XTick',[1:9]);
set(gca,'XTickLabel',[ {'hydroxyurea'},{'NaCl'},{'allantoin'},{'caffeine'},{'galactose'},{'glycine'},{'isoleucine'},{'phleomycin'},{'rapamycin'}]);

set(gca,'XTickLabelRotation',45)

ylabel([{'Fraction of {\itK} variance'}, {'explained by 36 SNPs'}])
set(gca,'YTickLabel',[{'0.0'},{'0.3'},{'0.6'},{'0.9'}]);
set(gca,'FontSize',20);
box on
text(-3,0.99,'b','FontSize',50)
A=importdata('R_K_strain_indexed.txt');
Data=A.data;
All_env_corr=[];
Q=[];
for i=1:2:18;
    data=Data(:,i:i+1);
    ind=find(data(:,1)>0);
    data=data(ind,:);
    [rho p]=corr(data(:,1),data(:,2),'type','Spearman');
    All_env_corr=[All_env_corr;rho,p];
    Q=[Q;mean(data)];
end


[Data txt]=xlsread('QTL36_fitting_result.xlsx');
R2=Data(:,1);
Data=Data(:,3:end);

Predict1=[];
Sign1=[];
for i=1:2:18;
    data1=Data(i,:);
    data2=Data(i+18,:);
    ind=find(isnan(data1)==0);
    Sign_test=data1.*data2;
    ind2=find(Sign_test<0);
    Sign1=[Sign1;length(ind2)];
    [rho p]=corr(data1(ind)',data2(ind)','type','Spearman');
    Predict1=[Predict1;rho,p];
end
Predict2=[];
Sign2=[];
for i=2:2:18;
    data1=Data(i,:);
    data2=Data(i+18,:);
    ind=find(isnan(data1)==0);
    Sign_test=data1.*data2;
    ind2=find(Sign_test<0);
    Sign2=[Sign2;length(ind2)];
    [rho p]=corr(data1(ind)',data2(ind)','type','Spearman');
    Predict2=[Predict2;rho,p];
end

subplot('position',[0.72 0.25 0.22 0.63])
hold on
[rho p]=corr(Q(:,1),Sign1,'type','Spearman')
find(Sign1/36<0.5)
find(Sign1/36==0.5)
find(Sign2/36<0.5)
[rho p]=corr(Q(:,1),Sign2,'type','Spearman')

plot(Q(:,1),Sign2/36,'sb','MarkerFaceColor','b','MarkerSize',10);
plot(Q(:,1),Sign1/36,'.r','MarkerSize',30);
set(gca,'FontSize',20);
xlim([0.06 0.16])
xlabel('Environment quality \itQ')
ylabel([{'Fraction of'},{'antagonistic QTLs'}])
set(gca,'XTick',[0.06:0.03:0.16])
set(gca,'YTick',[0.2:0.2:1])
set(gca,'YTickLabel',[{'0.2'},{'0.4'},{'0.6'},{'0.8'},{'1.0'}])
box on
text(0.035,1.08,'c','FontSize',50);
text(0.095, 0.48, '{\it\rho} = 0.94, {\itP} = 4.9\times10^{-4}', 'color','r','FontSize',18)
text(0.095, 0.38, '{\it\rho} = 0.74, {\itP} = 0.018', 'color','b','FontSize',18)

set(gca,'FontSize',20);
set(gcf,'paperposition',[0 0 18 6]);
saveas(1,'Fig3','png');
