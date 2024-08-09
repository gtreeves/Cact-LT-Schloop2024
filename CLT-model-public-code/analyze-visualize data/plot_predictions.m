clc
clear
close all


load("two_component_model_simple_normalized_GFPsd_rng2.mat","datnew")
load("predictions3.mat")
load("phi.mat")
load("NCRI.mat")


%%
figure
y = [0 1 1 0];

c1(1) = cdfplot(phi_wt_C);
hold on
xline(mean(phi_1x1x),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
x1 = [min(phi_1x1x) min(phi_1x1x) max(phi_1x1x) max(phi_1x1x)];
fill(x1,y,[0.9290 0.6940 0.1250],'FaceAlpha',0.25)

c1(2) = cdfplot(phi_2xGFP_C);
xline(mean(phi_2xGFP),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
x1 = [min(phi_2xGFP) min(phi_2xGFP) max(phi_2xGFP) max(phi_2xGFP)];
fill(x1,y,[0.8500 0.3250 0.0980],'FaceAlpha',0.25)

c1(3) = cdfplot(phi_2xLT_C);
xline(mean(phi_2xLT),'Color',[0 0.4470 0.7410],'LineWidth',2)
x1 = [min(phi_2xLT) min(phi_2xLT) max(phi_2xLT) max(phi_2xLT)];
fill(x1,y,[0 0.4470 0.7410],'FaceAlpha',0.25)

set(c1,'LineWidth',2,{'Color'},{[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]}',{'LineStyle'},{'--','--','--'}')
lgd = legend('1xGFP1xLT','','','2xGFP','','','2xLT','','','Location','best');
title(lgd,['—- Model',' - Data'])
title('\Phi')
set(gca,'FontSize',20,'FontName','Arial')


%%
figure
y = [0 1 1 0];

c1(1) = cdfplot(NCRI_1x1x_C);
hold on
xline(mean(NCRI_1x1x),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
x1 = [min(NCRI_1x1x) min(NCRI_1x1x) max(NCRI_1x1x) max(NCRI_1x1x)];
fill(x1,y,[0.9290 0.6940 0.1250],'FaceAlpha',0.25)

c1(2) = cdfplot(NCRI_2xGFP_C);
xline(mean(NCRI_2xGFP),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
x1 = [min(NCRI_2xGFP) min(NCRI_2xGFP) max(NCRI_2xGFP) max(NCRI_2xGFP)];
fill(x1,y,[0.8500 0.3250 0.0980],'FaceAlpha',0.25)

c1(3) = cdfplot(NCRI_2xLT_C);
xline(mean(NCRI_2xLT),'Color',[0 0.4470 0.7410],'LineWidth',2)
x1 = [min(NCRI_2xLT) min(NCRI_2xLT) max(NCRI_2xLT) max(NCRI_2xLT)];
fill(x1,y,[0 0.4470 0.7410],'FaceAlpha',0.25)

set(c1,'LineWidth',2,{'Color'},{[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]}',{'LineStyle'},{'--','--','--'}')
lgd = legend('1xGFP1xLT','','','2xGFP','','','2xLT','','','Location','best');
title(lgd,['—- Model',' - Data'])
title('NCRI')
set(gca,'FontSize',20,'FontName','Arial')
