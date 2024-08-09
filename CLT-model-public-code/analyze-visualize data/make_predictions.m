clc
clear
close all

load("../Results/two_component_model_simple_normalized_GFPsd_rng2.mat","datnew")
% load("predictions2.mat")
phi_wt = 0.488107;
phi_2xLT = 0.611322;
phi_2xGFP = 0.207205;
NCRI_1x1x = 0.273980;
NCRI_2xLT = 0.331515;
NCRI_2xGFP = 0.502820;
NCRI_free_GFP = 1.244363; 
% KnucG = NCRI_free_GFP;

% Gtot_exp = 30;
sd_phi_wt = 0.156322;
sd_phi_2xLT = 0.059251;
sd_phi_2xGFP = 0.103427;
sd_NCRI_1x1x = 0.031633;
sd_NCRI_2xLT = 0.026768;
sd_NCRI_2xGFP = 0.109118;
sd_NCRI_free_GFP = 0.068004;

stage = 'interphase'; m = 51;
[~,~,~,Vn,Vc] = nuclearSize(1,'static',m,stage);
Vtot = Vn + Vc;

for i=1:length(datnew)

    xb = datnew(i).xb; XB(i,:) = xb;
    g = datnew(i).Gtot; G(i) = g;
    k = datnew(i).KnucG; K(i) = k;
    Y0 = datnew.Y0;

    %1x1x
    [phi_wt_c,NCRI_1x1x_c,Y,~,e1,~] = calc_phi_NCRI(xb,Y0,g,k);
    Cactc(i) = Y(4) + Y(5);
    Cactn(i) = (Vtot*xb(5) - Vc*Cactc(i))/Vn;

    phi_wt_C(i) = phi_wt_c; NCRI_1x1x_C(i) = NCRI_1x1x_c;

    %2xGFP
    [phi_2xGFP_c,NCRI_2xGFP_c,~,~,e2,~] = calc_phi_NCRI(xb,Y0,2*g,k);
    phi_2xGFP_C(i) = phi_2xGFP_c; NCRI_2xGFP_C(i) = NCRI_2xGFP_c;

    %2xLT
    is2xLT = 1;
    [phi_2xLT_c,NCRI_2xLT_c,~,~,e3,~] = calc_phi_NCRI(xb,Y0,g,k,is2xLT);
    phi_2xLT_C(i) = phi_2xLT_c; NCRI_2xLT_C(i) = NCRI_2xLT_c;

end
save('predictions.mat')
%%
load("predictions.mat")
figure
y = [0 1 1 0];

c1(1) = cdfplot(phi_wt_C);
hold on
xline(phi_wt,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
x1 = [phi_wt-sd_phi_wt phi_wt-sd_phi_wt phi_wt+sd_phi_wt phi_wt+sd_phi_wt];
fill(x1,y,[0.9290 0.6940 0.1250],'FaceAlpha',0.25)

c1(2) = cdfplot(phi_2xGFP_C);
xline(phi_2xGFP,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
x1 = [phi_2xGFP-sd_phi_2xGFP phi_2xGFP-sd_phi_2xGFP phi_2xGFP+sd_phi_2xGFP phi_2xGFP+sd_phi_2xGFP];
fill(x1,y,[0.8500 0.3250 0.0980],'FaceAlpha',0.25)

c1(3) = cdfplot(phi_2xLT_C);
xline(phi_2xLT,'Color',[0 0.4470 0.7410],'LineWidth',2)
x1 = [phi_2xLT-sd_phi_2xLT phi_2xLT-sd_phi_2xLT phi_2xLT+sd_phi_2xLT phi_2xLT+sd_phi_2xLT];
fill(x1,y,[0 0.4470 0.7410],'FaceAlpha',0.25)

set(c1,'LineWidth',2,{'Color'},{[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]}',{'LineStyle'},{'--','--','--'}')
lgd = legend('1xGFP1xLT','','','2xGFP','','','2xLT','','','Location','best');
title(lgd,['—- Model',' - Data'])
title('\Phi')
set(gca,'FontSize',16)

%%
figure
y = [0 1 1 0];

c1(1) = cdfplot(NCRI_1x1x_C);
hold on
xline(NCRI_1x1x,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
x1 = [NCRI_1x1x-sd_NCRI_1x1x NCRI_1x1x-sd_NCRI_1x1x NCRI_1x1x+sd_NCRI_1x1x NCRI_1x1x+sd_NCRI_1x1x];
fill(x1,y,[0.9290 0.6940 0.1250],'FaceAlpha',0.25)

c1(2) = cdfplot(NCRI_2xGFP_C);
xline(NCRI_2xGFP,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
x1 = [NCRI_2xGFP-sd_NCRI_2xGFP NCRI_2xGFP-sd_NCRI_2xGFP NCRI_2xGFP+sd_NCRI_2xGFP NCRI_2xGFP+sd_NCRI_2xGFP];
fill(x1,y,[0.8500 0.3250 0.0980],'FaceAlpha',0.25)

c1(3) = cdfplot(NCRI_2xLT_C);
xline(NCRI_2xLT,'Color',[0 0.4470 0.7410],'LineWidth',2)
x1 = [NCRI_2xLT-sd_NCRI_2xLT NCRI_2xLT-sd_NCRI_2xLT NCRI_2xLT+sd_NCRI_2xLT NCRI_2xLT+sd_NCRI_2xLT];
fill(x1,y,[0 0.4470 0.7410],'FaceAlpha',0.25)

set(c1,'LineWidth',2,{'Color'},{[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]}',{'LineStyle'},{'--','--','--'}')
lgd = legend('1xGFP1xLT','','','2xGFP','','','2xLT','','','Location','best');
title(lgd,['—- Model',' - Data'])
title('NCRI')
set(gca,'FontSize',16)
%%
figure;
histogram(XB(:,5),'Normalization','probability','BinWidth',5)
hold on; histogram(Cactc,'Normalization','probability','BinWidth',5)
hold on;histogram(Cactn,'Normalization','probability','BinWidth',5)
% xticks([0 5 20 50 100 150 200 250])
legend('Cact_{tot}','Cact_{cyt}','Cact_{nuc}')
set(gca,'FontSize',14)