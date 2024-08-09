
clc
clear
close all

addpath /scratch/user/razeen/CLT_Gtot/ISRES_plus/
addpath /scratch/user/razeen/CLT_Gtot/ISRES_plus/Algorithm/
addpath /scratch/user/razeen/CLT_Gtot/ISRES_plus/Functions/
addpath /scratch/user/razeen/CLT_Gtot/ISRES_plus/Methods/
addpath /scratch/user/razeen/CLT_Gtot/ISRES_plus/Methods/Hessian/

%% Data 


%
% CONTROLS
%
modelname           = 'CLT_two_component_simple';

% evo options
evo.G               = 1000;        % maximum number of generations
evo.lambda          = 150;          % population size (number of offspring) (100 to 400) 
evo.mue             = 20;          % ~round((lambda)/7)
evo.nIslands        = 1;           % number of islands
evo.migGen          = -1;           % ~10 times migration is executed

evo.alphaa          = 0.2;         % Smoothing factor
evo.gammaa          = 0.85;        % Recombination parameter
evo.pf              = .45;         % pressure on fitness in [0 0.5] try 0.45
evo.varphi          = 1;           % expected rate of convergence (usually 1)
evo.tmax            = 72*3600;     % hrs * seconds/hr
evo.mm              = 'min';       % Minimize ('min') or Maximize('max')

% plus options
plus.yeslin         = true;
plus.yesnew         = true;
plus.nlin           = 2;        % 1-2       % 2  
plus.nnewt          = 1;        % 1-2       % 1
plus.nlinPar        = 1;        % 1-2       % 1
plus.nnewtPar       = 1;        % 1-2       % 1
plus.startlin       = 1;
plus.endlin         = evo.G;
plus.startnew       = 1;
plus.endnew         = evo.G;
plus.useFullNewtonStep      = true;     % applies only to newton step.
plus.sortPrevParamsByError  = false;     % applies to both linstep and newton step
plus.betalin = 1;

% model options
model.modelname     =  modelname;

% other options
options.reshuffGen  = 0; 
options.liveUpdates = true;
options.restart     = false;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% Create options super-structure
options.evo         = evo;
options.plus        = plus;
%options.model       = model;


%
% Run isres
%

rng('shuffle')                  
randomseed = rng;

a = (15-0.1)*rand + 0.1;
b = (15-0.1)*rand + 0.1;
c = (15-0.1)*rand + 0.1;

Y0 = [a b c];
data.Y0 = Y0;

%
Gtot = (72.5 - 23)*rand + 23;%23.3 - 72.2
KnucG = (1.3 - 1.1)*rand + 1.1;
data.Gtot = Gtot;
data.KnucG = KnucG;

fhandle  = @(k)objfcn_two_component_simple(k,data);

% p = [kin kout kon KD Ctot]
%min^-1
lb = [0.016  0.1676 1e-3  1e-2   1];
ub = [0.0680  0.3101  1e2   1e2   1e3];

lu = [lb;ub];
lu       = log10(lu);

[xb,BestMin,Gm,Stats,opts] = isres_plus(fhandle,lu,options);
xb = 10.^xb;
BestMin;
Gm;
Stats;
opts;

if ~exist("Results","dir")
    mkdir("Results")
end
filename = ['.',filesep,'Results',filesep,'Results_isres-plus_clt',datestr(now,'yyyy-mm-dd-HH-MM-SS'),'-',char(randi([97,106],1,2)),num2str(randi(10))];
save([filename,'.mat'])   