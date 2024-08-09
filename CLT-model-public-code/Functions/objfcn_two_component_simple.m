%%
function [F,Phi] = objfcn_two_component_simple(k,data)
% clc
% clear
% close all

k = 10.^k;
nsets = size(k,1);
F     = zeros(nsets,1);
Phi   = zeros(nsets,1); %no penalty

%% Data

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
Y0 = data.Y0;
Gtot_exp = data.Gtot;
KnucG = data.KnucG;

%% Volumes



%parameters
% p = [kin kout kon koff Ctot]
for i = 1:nsets
% i = 1;
% k = rand(1,5);

    kin = k(i,1);
    kout = k(i,2);
    kon = k(i,3);
    KD = k(i,4);
    Ctot = k(i,5);
    %

    p = [kin kout kon KD Ctot];
    
    %1x1x
    [phi_wt_c,NCRI_1x1x_c,~,~,e1,~] = calc_phi_NCRI(p,Y0,Gtot_exp,KnucG);
    %2xGFP
    [phi_2xGFP_c,NCRI_2xGFP_c,~,~,e2,~] = calc_phi_NCRI(p,Y0,2*Gtot_exp,KnucG);
    %2xLT
    is2xLT = 1;
    [phi_2xLT_c,NCRI_2xLT_c,~,~,e3,~] = calc_phi_NCRI(p,Y0,Gtot_exp,KnucG,is2xLT);

    if (e1>0) && (e2>0) && (e3>0)
        F(i) = ((phi_wt - phi_wt_c)/sd_phi_wt)^2 + ((phi_2xGFP - phi_2xGFP_c)/sd_phi_2xGFP)^2 + ((phi_2xLT - phi_2xLT_c)/sd_phi_2xLT)^2 + ...
            ((NCRI_1x1x - NCRI_1x1x_c)/sd_NCRI_1x1x)^2 + ((NCRI_2xGFP - NCRI_2xGFP_c)/sd_NCRI_2xGFP)^2 + ((NCRI_2xLT - NCRI_2xLT_c)/sd_NCRI_2xLT)^2;
    else
        F(i) = 1e5;
    end

end

end