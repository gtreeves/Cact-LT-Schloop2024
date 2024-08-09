clc
clear
close all

load("Results/two_component_model_simple_normalized_GFPsd_rng2.mat","datnew")
load("Results/predictions.mat")


for i=1:length(datnew)

    XB(i,:) = datnew(i).xb;

end

%%

figure

histogram(log10(XB(:,3)),'Normalization','probability','NumBins',10,'FaceColor',[0.6 0.6 0.6])
xlabel('log_{10} k_{on} [nM^{-1} min^{-1}]')
ylabel('pdf')
xticks(-1.8:0.1:-1.2)
xlim([-1.8 -1.2])
ylim([0 0.3])
set(gca,'FontSize',20,'FontName','Arial')

figure

histogram(XB(:,4),'Normalization','probability','NumBins',10,'FaceColor',[0.6 0.6 0.6])
xlabel('K_{D} [nM]')
ylabel('pdf')
xticks(5:5:30)
xlim([5 30])
ylim([0 0.3])
set(gca,'FontSize',20,'FontName','Arial')

figure

histogram(Cactc,'Normalization','probability','NumBins',10,'FaceColor',"#4DBEEE", ...
    FaceAlpha=1)
hold on
histogram(XB(:,5),'Normalization','probability','NumBins',10,'FaceColor',[0.6 0.6 0.6], ...
    FaceAlpha=0.5)

legend('cytoplasmic','total')
xlabel('Cactus [nM]')
ylabel('pdf')
xticks(25:50:250)
xlim([25 250])
ylim([0 0.3])

set(gca,'FontSize',20,'FontName','Arial')

figure

histogram(Cactn,'Normalization','probability','NumBins',10,'FaceColor',"#EDB120", ...
    FaceAlpha=1)


xlabel('Cactus(nuclear) [nM]')
ylabel('pdf')
xticks(5:5:25)
xlim([5 30])
ylim([0 0.3])

set(gca,'FontSize',20,'FontName','Arial')

%% figure 2

figure

histogram(XB(:,5),'Normalization','probability','NumBins',10,'FaceColor',[0.6 0.6 0.6], ...
    FaceAlpha=0.5)


xlabel('Total Cactus [nM]')
ylabel('pdf')
xticks(50:50:200)
xlim([50 200])
ylim([0 0.3])

set(gca,'FontSize',20,'FontName','Arial')

figure

histogram(Cactn./Cactc,'Normalization','probability','BinWidth',0.0015,'FaceColor',[0.6 0.6 0.6], ...
    FaceAlpha=0.5)


xlabel('NCR of Cactus')
ylabel('pdf')
xticks(0.07:0.02:0.16)
xlim([0.07 0.16])
ylim([0 0.3])

set(gca,'FontSize',20,'FontName','Arial')