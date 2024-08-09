function [xlin,etalin,stats] = step_lin(xbub,fbub,xPar,etalin,nlin,lu,betalin)
% This function calculates the best-fit hyper-plane surface (dimensionality
% n, where n = number of parameters) using a group of "np" parameter
% sets and their corresponding objective function values, "f".

%load('linstep.mat')

warning('on','verbose')
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:singularMatrix

if ~exist('betalin','var')
    betalin = 1;
end

[np,n]   = size(xbub);
mu       = size(xPar,1);
% xPar     = 10.^xPar;
% xbub     = 10.^xbub;
% lb       = 10.^lu(1,:);
% ub       = 10.^lu(2,:);

lb       = lu(1,:);
ub       = lu(2,:);

% Calculate linstep direction
flag        = false;
nInd        = n+1;
nretry      = 4;
retry       = 1;
xincParamBy = linspace(0,2,nretry);
nIndUsed    = nInd;
while ~flag && retry<nretry   
    A       = [xbub(1:nIndUsed,:) ones(nIndUsed,1)];
    lastwarn('')
    lam     = A \ fbub(1:nIndUsed);
    [~, warnID] = lastwarn;
    if contains(warnID,'Matrix') || any(isnan(lam))
       nIndUsed  = round(nInd*exp(xincParamBy(retry)));
       retry     = retry + 1;
       if nIndUsed > np || retry >= nretry
           xlin                     = [];
           stats.nlin               = 0;
           stats.nIndUsed           = nIndUsed;
           stats.retry              = retry;
           stats.nBlocks            = 0;
           stats.gradStepLen        = nan;
           stats.Dpar               = nan;
           stats.betamax            = nan;
           stats.nCorrectedBounds   = nan;
           stats.reason4failure     = warnID;
           disp(['Linstep failed. nIndUsed:', num2str(nIndUsed)])
           return;
       end
    else
        flag = true;
    end  
end
fdot     = lam(1:n);

disp(['# of individuals used in successful  linstep is: ', num2str(nIndUsed)])

% Calculate how much distance to translate linstep parents
D        = sum((xbub(1:nIndUsed,:)-xPar(1,:)).^2,2).^0.5;
D(D==0)  = [];
Dpar     = mean(D);


% Now the linear explanatory model suggests that the gradient of f is
% most steeply descending in the "-fdot" direction.
Q       = sum(fdot.^2).^0.5;
grad    = Dpar*betalin*(-fdot'/norm(fdot));         % direction w magnitude!


% Ready inidividuals for Linstep
nBlocks    = ceil(nlin/mu);
xPar       = repmat(xPar,nBlocks,1);
xPar       = xPar(1:nlin,:);
etalin     = repmat(etalin,nBlocks,1);
etalin     = etalin(1:nlin,:);
grad       = repmat(grad,nlin,1);


% Calculate: At what distance in the linstep direction do the parents
% hit the boundary following the gradient?
llimit   = (repmat(lb,nlin,1)- xPar)./(grad);
ulimit   = (repmat(ub,nlin,1)- xPar)./(grad);
betamax  = zeros(nlin,1);
for i=1:nlin
    lim         = llimit(i,:);
    lbeta       = min(lim(lim > 0));
    lim         = ulimit(i,:);
    ubeta       = min(lim(lim > 0));
    if isempty(lbeta)
       lbeta = inf; 
    end
    if isempty(ubeta)
        ubeta = inf; 
    end
    betamax(i) = min(lbeta,ubeta);
end     
% The individuals that have betamax<1 will be unsuccessful in taking a full
% step
v = betamax < 1; 



%
% LINSTEP
%
% Note that the direction is "-fdot". Therefore, we use xnew = xold +
% grad*d instead of xnew = xold - grad*d.
betarand    = rand(nlin,1);
xlin        = xPar + grad;
if any(v)
    xlin(v,:)  = xPar(v,:) + betarand(v).*betamax(v).*grad(v,:);  
end
%xlin        = log10(xlin);
nlin        = size(xlin,1);


% Linstep stats
stats.nlin              = nlin;
stats.nIndUsed          = nIndUsed;
stats.retry             = retry - 1;
stats.nBlocks           = nBlocks;
stats.gradStepLen       = Q;
stats.Dpar              = Dpar;
stats.betamax           = betamax;
stats.nCorrectedBounds  = length(find(v));
stats.reason4failure    = '';








