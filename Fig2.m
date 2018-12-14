clear all
close all
clc
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
ind_pos = find(All_env_corr(:,1)>0);
A.textdata(1,ind_pos*2)
subplot(2,2,1)

[Max ind ] = max(Q);
A.textdata(1,ind*2:ind*2+1)
data=Data(:,ind*2-1:ind*2);
All_env_corr(ind,:)
ind=find(data(:,1)>0);
data=data(ind,:);
plot(data(:,1),data(:,2),'.k')
xlim([0.1 0.2])
set(gca,'FontSize',25,'LineWidth',1);
xlabel('\itr');
ylabel('\itK');
text(0.07,4*10^6,'a','FontSize',50);
text(0.14, 4.48*10^6,'Allantoin','FontSize',30)
box off
text(0.12,4*10^6,'{\it\rho}_{rK} = -0.52, {\itP} < 10^{-250}','FontSize',23);
subplot(2,2,2)
A.textdata(1,ind*2:ind*2+1)
[Min ind ] = min(Q);
All_env_corr(ind,:)
data=Data(:,ind*2-1:ind*2);
ind=find(data(:,1)>0);
data=data(ind,:);
plot(data(:,1),data(:,2),'.k')
xlim([0.03 0.1])
set(gca,'FontSize',25,'LineWidth',1);
xlabel('\itr');
ylabel('\itK');
text(0.05,3*10^6,'{\it\rho}_{rK} = 0.32, {\itP} < 10^{-166}','FontSize',23);
text(0.05,3.4*10^6,'Caffeine','FontSize',30)
text(0.007,3*10^6,'b','FontSize',50);
box  off
set(gca,'YTick',[0:1:3]*10^6)


subplot(2,2,3)
[rho p]=corr(All_env_corr(:,1),Q,'type','Spearman')
plot(Q(:,1),All_env_corr(:,1),'.k','MarkerSize',40);
set(gca,'FontSize',25,'LineWidth',1);
xlim([0.06 0.16])
ylim([-0.6 0.6])
xlabel('Environment quality ({\itQ})');
ylabel('{\it\rho}_{rK}');
text(0.025,0.6,'c','FontSize',50);
set(gca,'YTick',[-0.6:0.3:0.6])
set(gca,'YTickLabel',[{'-0.6'},{'-0.3'},{'0.0'},{'0.3'},{'0.6'}]);
set(gca,'XTick',[0.06:0.03:0.16])
box off
text(0.08,0.61,'{\it\rho} = -0.88, {\itP} = 0.0031','FontSize',23);

subplot(2,2,4)
plot(Q(:,2),All_env_corr(:,1),'.k','MarkerSize',40);
set(gca,'FontSize',25,'LineWidth',1);
xlabel('Mean {\itK} of all genotypes');
ylabel('{\it\rho}_{rK}');
xlim([0.5*10^6 3.5*10^6])
set(gca,'XTick',[1:3]*10^6);
set(gca,'XTickLabel',[{'1'},{'2'},{'3'}])
text(3.5*10^6,-0.7,'\times 10^{6}','FontSize',25)
ylim([-0.6 0.6])
text(-0.6*10^6,0.6,'d','FontSize',50);
text(10^6,0.61,'{\it\rho} = 0.18, {\itP} = 0.64','FontSize',23);
set(gcf,'paperposition',[0 0 14 14]);
set(gca,'YTick',[-0.6:0.3:0.6])
set(gca,'YTickLabel',[{'-0.6'},{'-0.3'},{'0.0'},{'0.3'},{'0.6'}]);
box off
saveas(1,'Fig2','png');
% Strain_E_corr=[];
% All_P=[];
% for i=1:row;
%     data1=Data(i,1:2:18);
%     data2=Data(i,2:2:18);
%     ind=find(data1~=-1);
%     [rho11 p1]=corr(data1(ind)',Q(ind,1),'type','Spearman');
%     [rho12 p2]=corr(data1(ind)',Q(ind,2),'type','Spearman');
%     [rho21 p3]=corr(data2(ind)',Q(ind,1),'type','Spearman');
%     [rho22 p4]=corr(data2(ind)',Q(ind,2),'type','Spearman');
%     Strain_E_corr=[Strain_E_corr;rho11,rho12,rho21,rho22];
%     All_P=[All_P;p1,p2,p3,p4];
% end
% xlswrite('Corr_strain_R_K_Q_acrossE.xlsx',[{'rhoR_QR'},{'rhoR_QK'},{'rhoK_QR'},{'rhoK_QK'},{'pR_QR'},{'pR_QK'},{'pK_QR'},{'pK_QK'}],1,'A1');
% xlswrite('Corr_strain_R_K_Q_acrossE.xlsx',Strain_E_corr,1,'A2');
% xlswrite('Corr_strain_R_K_Q_acrossE.xlsx',All_P,1,'E2');


    
