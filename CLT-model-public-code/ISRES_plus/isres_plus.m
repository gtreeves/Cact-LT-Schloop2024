function [xb,BestMin,Gm,Stats,options] = isres_plus(fcn,lu,options)
% ISRESI "Improved" Evolution Strategy using Stochastic Ranking & Island
% hopping
% usage : 
%   [xb,Gm,Mean,Min,X,F,Phi,Isort,Eta] = isresi(fcn,nIslands,lu,lambda,G,...
%                                               mu,pf,varphi,tmax,varargin)
% where
%       nIslands         : number of islands where isres is run
%
% ______________________ORIGINAL DESCRIPTION_______________________________ 
%
% ISRES "Improved" Evolution Strategy using Stochastic Ranking
% usage:
%        [xb,Stats,Gm] = isres(fcn,mm,lu,lambda,G,mu,pf,varphi);
% where
%        fcn       : name of function to be optimized (string)
%        mm        : 'max' or 'min' (for maximization or minimization)
%        lu        : parameteric constraints (lower and upper bounds)
%        lambda    : population size (number of offspring) (100 to 400)
%        G         : maximum number of generations
%        mu        : parent number (mu/lambda usually 1/7)
%        pf        : pressure on fitness in [0 0.5] try 0.45
%        varphi    : expected rate of convergence (usually 1)
%
%        xb        : best feasible individual found
%        Stats     : [min(f(x)) mean(f(x)) number_feasible(x)]
%        Gm        : the generation number when "xb" was found
%
% Copyleft (C) 2003-2004 Thomas Philip Runarsson (e-mail: tpr@hi.is)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% PB: 03/14/2020
% Code now implements a).Island Hopping and b). a Multiparent strategy - for
% recombination of the top mu individuals, it uses the top 3 individuals as
% parents.
% PB: 07/15/2020
% Current version implements lin_step as an option. 

%
% Default values of algorithm options if options is not specified
%
evo.G               = 75;           % maximum number of generations
evo.lambda          = 250;          % population size (number of offspring) (100 to 400) 
evo.alphaa          = 0.2;          % Smoothing factor
evo.gammaa          = 0.85;         % Recombination parameter
evo.mue             = 35;
evo.nIslands        = 4;            % number of islands
evo.migGen          = 10;
evo.pf              = .45;          % pressure on fitness in [0 0.5] try 0.45
evo.varphi          = 1;            % expected rate of convergence (usually 1)
evo.tmax            = 72*3600;      % hrs * seconds/hr
evo.mm              = 'min';        % Minimize ('min') or Maximize('max')
plus.yeslin         = true;
plus.yesnew         = true;
plus.nlin           = 3;
plus.nnewt          = 1;
plus.nlinPar        = 1;
plus.nnewtPar       = 1;
plus.useFullNewtonStep      = true;
plus.sortPrevParamsByError  = true;
plus.startnew       = 1;
plus.endnew         = evo.G;
plus.startlin       = 1;
plus.endlin         = evo.G;
plus.betalin        = 1;
model.modelname     = [];
model.addParams     = [];

reshuffGen          = 0;
liveUpdates         = true;
restart             = false;

% Extract algorithm and model options from options
if exist('options','var')
    struct2vars(options);
end
struct2vars(evo);
struct2vars(plus);
struct2vars(model);
clear options evo plus model


%
% Check algorithm parameters for consistency
%
if nIslands == 1, migGen = 0;   end
if ~yeslin && ~yesnew   
   nlin     = 0;
   nnewt    = 0;
   nlinPar  = -1;
   nnewtPar = -1;
   startnew = -1;
   endnew   = -1;
   startlin = -1;
   endlin   = -1;
   useFullNewtonStep = false;
   sortPrevParamsByError = false;
elseif yeslin && ~yesnew
   nnewt    = 0;
   nnewtPar = -1;
   startnew = -1;
   endnew   = -1;
   useFullNewtonStep = false;
elseif ~yeslin && yesnew 
   nlin     = 0;
   nlinPar  = -1;
   startlin = -1;
   endlin   = -1;
end
if nlinPar  > nlin,  nlinPar  = nlin; end
if nnewtPar > nnewt, nnewtPar = nnewt; end


if      strcmpi(mm,'max'), mm = -1; 
elseif  strcmpi(mm,'min'), mm = 1;
else,   error('mm should either be "min" or "max"')
end


warning off MATLAB:ode15s:IntegrationTolNotMet
rng('shuffle')                  
randomseed = rng;


% Define search parameters 
n        = size(lu,2);                    % Number of parameters in a set
chi      = (1/(2*n)+1/(2*sqrt(n)));
varphi   = sqrt((2/chi) * log((1/alphaa)*(exp(varphi^2*chi/2)-(1-alphaa))));
tau      = varphi/(sqrt(2*sqrt(n)));      % learning rate for a parameter 
tau_     = varphi/(sqrt(2*n));            % learning rate for an individual 
ub       = ones(lambda,1)*lu(2,:);        % Upper bound on parameters
lb       = ones(lambda,1)*lu(1,:);        % Lower bound on parameters                    
sI       = (1:mue)'*ones(1,ceil(lambda/(mue))); 
sI       = sI(1:lambda);


% Initialize population and other associated variables
if restart
   file = extractFileLocations('./','mat',false);
   file = file(contains(file,'restart'));
   if length(file) > 1
      error('There can only be one restart mat file') 
   else
       load(file,'x','eta')
       eta_u  = (lu(2,:)-lu(1,:))/sqrt(n);
       lambda = size(x,1)/nIslands;
       if rem(lambda,1) ~= 0, error('Restart problem: Lambda is not a whole number'), end
   end
else
    x       = ones(nIslands*lambda,1)*lu(1,:) + rand(nIslands*lambda,n).* ...
              (ones(nIslands*lambda,1)*(lu(2,:)-lu(1,:)));
    eta     = ones(nIslands*lambda,1)*(lu(2,:)-lu(1,:))/sqrt(n);
    eta_u   = eta(1,:);
end

% divide population into islands
x               = divideIntoIslands(x,nIslands);
x_add           = repmat(addParams,nIslands*lambda,1);


%
% Store variables for analysis
%
X                    = zeros(nIslands*lambda,n,G);        
Eta                  = X;
F                    = zeros(nIslands*lambda,G);         % Errors
Isort                = F;                                % Rankings
Min                  = zeros(G,nIslands);   
Mean                 = Min;    
bestID               = Min; 
NrFeas               = Min; 
bestMethod           = strings(G,nIslands);
bestIsland           = nan(G,nIslands);
BestMinI             = inf*ones(nIslands,1); 
xbI                  = nan(nIslands,n);
etabI                = xbI;
%
% evo params
Nreco                = zeros(G,nIslands); 
Nmuta                = Nreco;
Sreco                = Nreco;
Smuta                = Nreco;
nRetryX              = nan(G,nIslands);
nRetryEta            = nRetryX;

%
% linstep params
vlin                 = false(G,nIslands);   
Nlin                 = Nreco;
Slin                 = Nreco;
nlinIndUsed4Grad     = Nreco;
nRetryLin            = nan(G,nIslands);
linBlocks            = nRetryLin;
linStepLen           = Nreco;
Dlin                 = Nreco;
betaLinMax           = cell(G,nIslands);
nLinCorrectedBounds  = Nreco;
linReason4Failure    = betaLinMax;

%
% newton step params
vnewt                = vlin;
Nnewt                = Nreco; 
Snewt                = Nreco;
nNewtIndUsed4Hess    = Nreco;
nRetryHessian        = nRetryLin;
newtBlocks           = nRetryLin;
eigratio             = nRetryLin;
isPD                 = vlin;
newtStepLen          = Nreco;
betaNewtMax          = betaLinMax;
nNewtCorrectedBounds = Nreco;
newtReason4Failure   = betaLinMax;


%
% other variables of importance
vstati               = 2*ones(lambda,nIslands);
bestMeth             = 'None';
BestMin              = inf;
nretry               = 20;
bestUpdated          = false;




%
% LOOP
%
t1 = tic;       tgen = 0;            
xb = [];        nl   = 0;    
Gm = [];        gprevparam = 1000;     
g  = 1;   
while (g <= G) && ((toc(t1)+tgen)<tmax)
    
    t2 = tic;
    if liveUpdates
        disp(strcat('gen:',num2str(g),'/',num2str(G)))
    end    
    % 
    % Run code
    %
    if isempty(modelname)
        [f, phi]        = fcn([x, x_add]);
    else
        [f, phi]        = fcn([x, x_add],modelname);
    end
    f               = mm*f;
    X(:,:,g)        = x;
	F(:,g)          = f; 
    if g == 1
        Phi = nan(nIslands*lambda,size(phi,2),G);              % Penalties
    end
	Phi(:,:,g)      = phi; 
    Eta(:,:,g)      = eta;
    phi(phi<=0)     = 0;
    phi             = sum(phi.^2,2);
    
    
    %
    % First, calculate performance across all islands
    %
    posI = 1:lambda;
    for i=1:nIslands    
        %
        % Get island statistics
        %
        xI    = x(posI,:);    
        etaI  = eta(posI,:);
        fI    = f(posI);
        phiI  = phi(posI);
        % 
        % Find best fit individual 
        %
        v           = isnan(fI);
        fI(v)       = inf;
        phiI(v,:)   = inf;
        phiI        = sum(phiI,2);
        Feasible    = find(phiI<=0);
        NrFeas(g,i) = length(Feasible);
        if ~isempty(Feasible)
            Mean(g,i)       = mean(fI(Feasible));     
            [minF,minID]    = min(fI(Feasible));
            minID           = Feasible(minID); 
            Min(g,i)        = minF;
            bestID(g,i)     = minID;               
            if minF < BestMinI(i)           % best individual in the island
                xbI(i,:)    = xI(minID,:);
                etabI(i,:)  = etaI(minID,:);
                BestMinI(i) = minF;
            end
            if minF < BestMin               % best individual in the population               
                xb          = xI(minID,:);
                BestMin     = minF;
                Gm          = g;              
                if vstati(minID,i) == 1,          bestMeth = 'Recombination';
                elseif vstati(minID,i) == 2,      bestMeth = 'Mutation';
                elseif vstati(minID,i) == 3,      bestMeth = 'linstep';
                elseif vstati(minID,i) == 4,      bestMeth = 'newton step';
                end
                bestI           = i;
                bestMethod(g,i) = bestMeth;
                bestIsland(g,i) = bestI;
                bestUpdated     = true;
            end
        else
            Mean(g,i) = NaN;  Min(g,i) = NaN;   bestID(g,i) = NaN;
        end
        %
        % Rank and score population on island
        %
        I              = stoch(f(posI),phi(posI),pf);
        Isort(posI,g)  = I;
        score          = scorepop(I(sI),vstati(:,i),lambda);
        Sreco(g,i)     = score(1); 
        Smuta(g,i)     = score(2);
        Slin(g,i)      = score(3);
        Snewt(g,i)     = score(4);                          
        if liveUpdates
           fprintf('Score: island %i - Reco: %.2f, Muta: %.2f, Lin: %.2f, Newt: %.2f\n',i,score(1),score(2),score(3),score(4)) 
        end       
        %
        % Prepare population for the next generation
        %
        x(posI,:)    = xI(I(sI),:);
        eta(posI,:)  = etaI(I(sI),:); 
        phi(posI,:)  = phiI(I(sI));
        f(posI)      = fI(I(sI));
        posI         = posI + lambda;
    end
    
    if liveUpdates && bestUpdated
        fprintf('Best fit individual - Island: %i, Method: %s, error: %.2f\n',bestI,bestMeth,mm*BestMin)
        bestUpdated  = false;
    end 
    
    if mod(g,reshuffGen)==0
        [x,eta,nmuta,retryX,retryEta] = reshufflepop(X,Phi,Eta,F,lambda,mue,nIslands,lu,eta_u,tau,tau_);
        nRetryEta(g,:)  = retryEta;
        nRetryX(g,:)    = retryX;
        Nmuta(g,:)      = nmuta/nIslands;
        vstati          = 2*ones(lambda,nIslands);
        Eta(:,:,g)      = eta; 
        g               = g + 1;  
        tgen            = toc(t2);
        fprintf('Population reshuffled\n\n')
        continue;
    end


    %
    % Build new population
    %
    posI = 1:lambda;
    for i=1:nIslands 
        %
        % Get island statistics
        %
        xI    = x(posI,:);         
        etaI  = eta(posI,:);
        %
        % Island hopping: Best individual in the island's population 
        % migrates to other islands
        %
        if mod(g,migGen) == 0      
           posR = 1:nIslands;
           posR = circshift(posR,-(i-1));
           posR = posR(2:end);
           posL = (mue-nIslands+2):mue;
           if ~any(any(isnan(xbI(posR,:)),2))
                xI(posL,:)   = xbI(posR,:);
                etaI(posL,:) = etabI(posR,:);
                xI           = xI(sI,:);
                etaI         = etaI(sI,:);  
           end
        end   
        x_                = xI;                % Create a copy of the present generation
        eta_              = etaI;
        etaI(mue:end,:)   = eta_(mue:end,:).*exp(tau_*randn(lambda-mue+1,1)*ones(1,n) + tau*randn(lambda-mue+1,n));
        %
        % Check upper bound on eta by retry      
        %
        v0      = any((etaI>eta_u)==true,2);
        retry   = 1;
        while any(v0>=1) && retry<=nretry
            ntimes       = length(find(v0 == true));
            etaI(v0,:)   = eta_(v0,:).*exp(tau_*randn(ntimes,1)*ones(1,n) + tau*randn(ntimes,n));
            v0           = any((etaI>eta_u)==true,2);
            retry        = retry + 1;
        end  
        nRetryEta(g,i)   = retry;
        for j=1:n 
            I         = find(etaI(:,j)>eta_u(j)); 
            etaI(I,j) = eta_u(j)*ones(size(I));
        end  
        %
        % Method 1: Recombination
        %
        xI(1:mue-1,:)     = xI(1:mue-1,:)  + gammaa*(ones(mue-1,1)*xI(1,:) - xI(2:mue,:));
        Nreco(g,i)        = mue-1;
        vstati(1:mue-1,i) = 1;
        %
        % Method 2: Mutation  
        %
        nmuta             = lambda - mue + 1;          
        xI(mue:end,:)     = xI(mue:end,:) + etaI(mue:end,:).*randn(nmuta,n);
        Nmuta(g,i)        = nmuta;
        vstati(mue:end,i) = 2;
        %
        % Update population by newton's step and linstep if asked for
        %
        if (yesnew || yeslin) && ~isempty(xbI(i,:)) && ((g>=startlin && g<=endlin) || (g>=startnew && g<=endnew))  
            if g<=gprevparam, fromgen = 1; else, fromgen = g-gprevparam; end
            [xbub,fbub,etabub]  = previous_params(X(:,:,fromgen:g),F(:,fromgen:g),Phi(:,:,fromgen:g),Eta(:,:,fromgen:g),xbI(i,:),sortPrevParamsByError);
            nViable             = size(xbub,1);
            %
            % Method: linstep 
            %
            if  nViable>=(n+1)  && yeslin && (g>=startlin && g<=endlin)
                if nViable < nlinPar
                    xPar   = xbub(1:nViable,:);
                    etaPar = etabub(1:nViable,:);
                else
                    xPar   = xbub(1:nlinPar,:);
                    etaPar = etabub(1:nlinPar,:);
                end
                [xmeth,etameth,statsmeth]   = step_lin(xbub,fbub,xPar,etaPar,nlin,lu,betalin);                             
                nl                          = statsmeth.nlin;
                vlin(g,i)                   = nl>0;
                Nlin(g,i)                   = nl;                 
                Nmuta(g,i)                  = Nmuta(g,i) - nl; 
                nlinIndUsed4Grad(g,i)       = statsmeth.nIndUsed;
                nRetryLin(g,i)              = statsmeth.retry;
                linBlocks(g,i)              = statsmeth.nBlocks;
                linStepLen(g,i)             = statsmeth.gradStepLen;
                Dlin(g,i)                   = statsmeth.Dpar;
                betaLinMax{g,i}             = statsmeth.betamax;
                nLinCorrectedBounds(g,i)    = statsmeth.nCorrectedBounds;
                linReason4Failure{g,i}      = statsmeth.reason4failure;
                if nl > 0
                    i0              = lambda-nl+1:lambda;
                    xI(i0,:)        = xmeth;
                    etaI(i0,:)      = etameth;
                    vstati(i0,i)    = 3;
                end
            end
            %
            % Method: Newton's step
            %
            if  nViable>=((n^2+n)/2 +1) && yesnew && (g>=startnew && g<=endnew) 
                if nViable < nnewtPar
                    xPar   = xbub(1:nViable,:);
                    etaPar = etabub(1:nViable,:);
                else
                    xPar   = xbub(1:nnewtPar,:);
                    etaPar = etabub(1:nnewtPar,:);
                end
                [xmeth,etameth,statsmeth]    = step_newton(xbub,fbub,xPar,etaPar,nnewt,lu,fcn,useFullNewtonStep);           
                nn                           = statsmeth.nnewt;
                vnewt(g,i)                   = nn>0; 
                Nnewt(g,i)                   = nn;  
                Nmuta(g,i)                   = Nmuta(g,i) - nn;
                nNewtIndUsed4Hess(g,i)       = statsmeth.nIndUsed;
                nRetryHessian(g,i)           = statsmeth.nRetryHessian;
                newtBlocks(g,i)              = statsmeth.nBlocks;
                eigratio(g,i)                = statsmeth.eigratio;
                isPD(g,i)                    = statsmeth.pd;
                newtStepLen(g,i)             = statsmeth.gradStepLen;
                betaNewtMax{g,i}             = statsmeth.betamax;
                nNewtCorrectedBounds(g,i)    = statsmeth.nCorrectedBounds;
                newtReason4Failure{g,i}      = statsmeth.reason4failure;
                if nn > 0
                    i0               = lambda-nl-nn+1:lambda-nl;
                    xI(i0,:)         = xmeth;
                    etaI(i0,:)       = etameth;
                    vstati(i0,i)     = 4;
                end
            end
        end
        
 
        % Correct individuals that are out of bounds parameter-wise.
        I     = find((xI>ub) | (xI<lb));
        retry = 1 ;
        while ~isempty(I)
            xI(I)   = x_(I) + etaI(I).*randn(length(I),1);
            I       = find((xI>ub) | (xI<lb));
            if (retry>nretry) 
                break; 
            end
            retry   = retry + 1;
        end
        if ~isempty(I)
            xI(I) = x_(I);              % ignore failures
        end
        nRetryX(g,i) = retry;

        
        % Reduce eta for next generations
        %nevo              = (mue-1) + Nmuta(g,i);
	    %etaI(mue:nevo,:)  = eta_(mue:nevo,:) + alphaa*(etaI(mue:nevo,:) - eta_(mue:nevo,:));
        etaI  = eta_ + alphaa*(etaI - eta_);
        
        % Update population and eta for next gen
        nl            = 0;
        x(posI,:)     = xI;         
        eta(posI,:)   = etaI;
        posI          = posI + lambda;
    end
    
    Eta(:,:,g)   = eta; 
    g            = g + 1;  
    tgen         = toc(t2);
    if liveUpdates
        fprintf('\n')
    end
end



%
% Stats struct
%
Stats.Mean                  = Mean;
Stats.Min                   = Min;
Stats.bestID                = bestID;
Stats.NrFeas                = NrFeas;
Stats.X                     = X;
Stats.F                     = F;
Stats.Phi                   = Phi;
Stats.Isort                 = Isort;
Stats.Eta                   = Eta;
Stats.BestMinI              = BestMinI;
Stats.xbI                   = xbI;
Stats.bestMethod            = bestMethod;

% Evolution stats 
Stats.evo.Nreco             = Nreco;
Stats.evo.Nmuta             = Nmuta;
Stats.evo.Sreco             = Sreco;
Stats.evo.Smuta             = Smuta;
Stats.evo.nRetryX           = nRetryX;
Stats.evo.nRetryEta         = nRetryEta;

% linstep stats
Stats.lin.vlin              = vlin;
Stats.lin.Nlin              = Nlin;
Stats.lin.Slin              = Slin;
Stats.lin.nBlocks           = linBlocks;
Stats.lin.nIndUsed4Grad     = nlinIndUsed4Grad;
Stats.lin.nRetry            = nRetryLin;
Stats.lin.steplen           = linStepLen;
Stats.lin.Dlin              = Dlin;
Stats.lin.betamax           = betaLinMax;
Stats.lin.nCorrectedBounds  = nLinCorrectedBounds;
Stats.lin.linReason4Failure = linReason4Failure;


% newton step stats
Stats.newt.vnewt            = vnewt;
Stats.newt.Nnewt            = Nnewt;
Stats.newt.Snewt            = Snewt;
Stats.newt.nBlocks          = newtBlocks;
Stats.newt.nIndUsed4Hess    = nNewtIndUsed4Hess;
Stats.newt.nRetryHessian    = nRetryHessian;
Stats.newt.eigratio         = eigratio;
Stats.newt.isPD             = isPD;
Stats.newt.steplen          = newtStepLen;
Stats.newt.betamax          = betaNewtMax;
Stats.newt.nCorrectedBounds = nNewtCorrectedBounds;
Stats.newt.newtReason4Failure = newtReason4Failure;


% metadata
Stats.meta.randomseed       = randomseed;



% 
% recreate options struct for output
%
if      mm == 1,  mm = 'min'; 
elseif  mm == -1, mm = 'max';
end
evo.G                       = G;              
evo.lambda                  = lambda;         
evo.alphaa                  = alphaa ;        
evo.gammaa                  = gammaa; 
evo.mue                     = mue;
evo.nIslands                = nIslands;       
evo.migGen                  = migGen;
evo.pf                      = pf;         
evo.varphi                  = varphi;           
evo.tmax                    = tmax;     
evo.mm                      = mm;       
plus.nlin                   = nlin;
plus.nnewt                  = nnewt;
plus.nlinPar                = nlinPar;
plus.nnewtPar               = nnewtPar;
plus.useFullNewtonStep      = useFullNewtonStep;
plus.sortPrevParamsByError  = sortPrevParamsByError;
plus.yeslin                 = yeslin;
plus.yesnew                 = yesnew;
plus.startnew               = startnew;
plus.endnew                 = endnew;
plus.startlin               = startlin;
plus.endlin                 = endlin;
plus.betalin                = betalin;
model.modelname             = modelname;
model.addParams             = addParams;


options.evo         = evo;
options.plus        = plus;
options.model       = model;
options.reshuffGen  = reshuffGen;
options.liveUpdates = liveUpdates;




