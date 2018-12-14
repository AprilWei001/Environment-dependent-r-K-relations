clear all
close all
clc

[Geno txt]=xlsread('SNPQTL36_genotype_all.xlsx');
SNP_all=Geno(1,:);
Geno=Geno(3:end,:);
[Phenotype txt]=xlsread('R_K_strain_indexed.xlsx');
Er = [];
for i = 1:9;
    ind=find(Phenotype(:,i*2-1)>=0);
    Er = [Er,mean(Phenotype(ind,i*2-1))];
end
maxEr = max(Er);
minEr = min(Er);
[X indEr]=sort(Er);


Envs =[{'HU'},{'NaCl'},{'all'},{'caf'},{'gal'},{'gly'},{'ile'},{'phl'},{'rap'}]
Colors = [1 0.2 0.2;
   0 0.5 0;
   0.4 0.4 1;
   1 0.6 0.2;
   1 1 0.2;
   1 0.4 1;
   1 0 0.5;
   0 1 0.5;
   0.6 0.2 1];
subplot('Position',[0.12 0.97 0.75 0.02])
hold on
for j=1:9;
    plot(j,2,'.','Color',Colors(j,:),'MarkerSize',70);
    text(j-0.2,0.5,char(Envs(j)),'FontSize',25)
end

box off
axis off

SNPS_inQTL36=[66, 986, 1544, 1613, 2971, 4450,  6585, 7198, 7574, 7844,  9358, 226, 3967, 4022, 4129, 5173, 9195, 10905];
for i=1:length(SNPS_inQTL36);
    snpi=SNPS_inQTL36(i);
    subplot(6,3,i)
    hold on
    for j=1:9;
        phenoR=Phenotype(:,j*2-1);
        phenoK=Phenotype(:,j*2);
        ind=find(phenoR~=-1);
        phenoR=phenoR(ind);
        phenoK=phenoK(ind);
        snp_ind=find(SNP_all==snpi);
        geno=Geno(ind,snp_ind(1));
        ind0=find(geno==0);
        ind2=find(geno==2);
        phenoR0=phenoR(ind0);
        phenoK0=phenoK(ind0);
        phenoR2=phenoR(ind2);
        phenoK2=phenoK(ind2); 
        sR = mean(phenoR2)-mean(phenoR0);
        sK = mean(phenoK2)-mean(phenoK0);
        seR = sqrt((var(phenoR2)+var(phenoR0))/(2/(1/length(ind0)+1/length(ind2))));
        seK = sqrt((var(phenoK2)+var(phenoK0))/(2/(1/length(ind0)+1/length(ind2))));
        plot(sR,sK,'.','Color',Colors(j,:),'MarkerSize',70);
        plot([sR-seR,sR+seR]',[sK,sK]','k','LineWidth',2);
        plot([sR,sR]',[sK-seK,sK+seK]','k','LineWidth',2);
    end
    set(gca,'FontSize',20);
    t=title(strcat('SNP ',num2str(SNPS_inQTL36(i))),'FontSize',20,'FontWeight','normal')
    set(t, 'horizontalAlignment', 'left')
    xlabel('Effect size on \itr','FontSize',19);
    ylabel('Effect size on \itK','FontSize',19);

end
subplot(6,3,1) %66

plot([-0.01 0.01]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-5 5]*10^5,'--b','LineWidth',1.5);
xlim([-0.008 0.008])
ylim([-5 5]*10^5)
set(gca,'XTick',[-0.008 0 0.008])
set(gca,'XTickLabel',[{'-0.008'},{'0'},{'0.008'}]);
text(-0.014,7*10^5,'a','FontSize',40)

subplot(6,3,2) %986

plot([-0.015 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[0 1.5]*10^6,'--b','LineWidth',1.5);
xlim([-0.015 0.005])
ylim([0 1.5]*10^6)
set(gca,'XTick',[-0.015 0.005])
set(gca,'XTickLabel',[{'-0.015'},{'0.005'}]);
text(-0.023,1.8*10^6,'b','FontSize',40)
subplot(6,3,3) %1544

plot([-0.015 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[0 1]*10^6,'--b','LineWidth',1.5);
xlim([-0.015 0.01])
ylim([0. .8]*10^6)
text(-0.024,.96*10^6,'c','FontSize',40)

subplot(6,3,4) % 1613
plot([-0.005 0.01]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-5 5]*10^5,'--b','LineWidth',1.5);
set(gca,'XTick',[-0.005 0 0.01])
set(gca,'XTickLabel',[{'-0.005'},{'0'},{'0.01'}]);
text(-0.0108,7*10^5,'d','FontSize',40)
xlim([-0.005 0.01])
ylim([-5 5]*10^5)

subplot(6,3,5) % 2971
plot([-0.01 0.01]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-4 4]*10^5,'--b','LineWidth',1.5);
xlim([-0.01 0.005])
ylim([-3 3]*10^5)
set(gca,'XTick',[-0.01 0 0.005])
set(gca,'XTickLabel',[{'-0.01'},{'0'},{'0.005'}]);
text(-0.016,4.2*10^5,'e','FontSize',40)

subplot(6,3,6) %4450
plot([-0.015 0.015]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-5 5]*10^5,'--b','LineWidth',1.5);
xlim([-0.005 0.015])
ylim([-2 5]*10^5)
set(gca,'XTick',[-0.005 0.015])
set(gca,'XTickLabel',[{'-0.005'},{'0.015'}]);
text(-0.0126,.64*10^6,'f','FontSize',40)

subplot(6,3,7) %6585
plot([-0.01 0.01]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-2 5]*10^5,'--b','LineWidth',1.5);
xlim([-0.01 0.01])
ylim([-1 5]*10^5)
text(-0.017,.64*10^6,'g','FontSize',40)

subplot(6,3,8)% 7198
plot([-0.02 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[0 1.5]*10^6,'--b','LineWidth',1.5);
xlim([-0.02 0.005])
ylim([0 1.5]*10^6)
set(gca,'XTick',[-0.02 0.005])
set(gca,'XTickLabel',[{'-0.02'},{'0.005'}]);
text(-0.0305,1.8*10^6,'h','FontSize',40)

subplot(6,3,9)%7574
plot([-0.01 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-3 5]*10^5,'--b','LineWidth',1.5);
xlim([-0.008 0.005])
ylim([-3 5]*10^5)
set(gca,'XTick',[-0.008 0 0.005])
set(gca,'XTickLabel',[{'-0.008'},{'0'},{'0.005'}]);
text(-0.0123,.66*10^6,'i','FontSize',40)

subplot(6,3,10)% 7844
plot([-0.007 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[0 3]*10^5,'--b','LineWidth',1.5);
xlim([-0.007 0.005])
ylim([-0 3]*10^5)
set(gca,'XTick',[-0.007 0 0.005])
set(gca,'XTickLabel',[{'-0.007'},{'0'},{'0.005'}]);
text(-0.011,.36*10^6,'j','FontSize',40)
subplot(6,3,11) %9358
plot([-0.001 0.006]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-3 3]*10^5,'--b','LineWidth',1.5);
xlim([-0.001 0.006])
ylim([-3 3]*10^5)
set(gca,'XTick',[-0.001 0.006])
set(gca,'XTickLabel',[{'-0.001'},{'0.006'}]);
text(-0.00345,.42*10^6,'k','FontSize',40)

subplot(6,3,12)
plot([-0.01 0.]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-0 5]*10^5,'--b','LineWidth',1.5);
xlim([-0.01 0.])
ylim([-0 5]*10^5)
text(-0.013, .6*10^6,'l','FontSize',40)
subplot(6,3,13)
plot([-0.01 0.002]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-0 6]*10^5,'--b','LineWidth',1.5);
xlim([-0.01 0.002])
ylim([-0 6]*10^5)
set(gca,'XTick',[-0.01 0.002])
set(gca,'XTickLabel',[{'-0.01'},{'0.002'}]);
text(-0.014,.72*10^6,'m','FontSize',40)
subplot(6,3,14)
plot([-0.003 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-0 6]*10^5,'--b','LineWidth',1.5);
xlim([-0.005 0.005])
ylim([-0 6]*10^5)

set(gca,'XTick',[-0.005 0 0.005])
set(gca,'XTickLabel',[{'-0.005'},{'0'},{'0.005'}]);
text(-0.0085, .72*10^6,'n','FontSize',40)
subplot(6,3,15)
plot([-0.015 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-0 4]*10^5,'--b','LineWidth',1.5);
xlim([-0.012 0.005])
ylim([-0 4]*10^5)
set(gca,'XTick',[-0.012 0 0.005])
set(gca,'XTickLabel',[{'-0.012'},{'0'},{'0.005'}]);
text(-0.018, 4.8*10^5,'o','FontSize',40)
subplot(6,3,16)
plot([-0.015 0.005]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-0 8]*10^5,'--b','LineWidth',1.5);
xlim([-0.015 0.005])
ylim([-0 8]*10^5)
set(gca,'XTick',[-0.015 0.005])
set(gca,'XTickLabel',[{'-0.015'},{'0.005'}]);
text(-0.022, .96*10^6,'p','FontSize',40)
subplot(6,3,17)
plot([-0.003 0.003]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-3 3]*10^5,'--b','LineWidth',1.5);
xlim([-0.003 0.003])
ylim([-3 3]*10^5)
set(gca,'XTick',[-0.003 0 0.003])
set(gca,'XTickLabel',[{'-0.003'},{'0'},{'0.003'}]);
text(-0.0053, .42*10^6,'q','FontSize',40)
subplot(6,3,18)
plot([-0.008 0.002]',[0 0]','--b','LineWidth',1.5);
plot([0 0]',[-2 4]*10^5,'--b','LineWidth',1.5);
xlim([-0.008 0.002])
ylim([-1 4]*10^5)
set(gca,'XTick',[-0.008 0.002])
set(gca,'XTickLabel',[{'-0.008'},{'0.002'}]);
text(-0.0115, .52*10^6,'r','FontSize',40)
set(gcf,'paperposition',[0 0 12 22]);
saveas(1,'Fig4','png');

