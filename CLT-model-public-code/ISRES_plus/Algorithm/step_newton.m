function [xnewt,etaPar,stats] = step_newton(xbub,fbub,xPar,etaPar,nnewt,lu,fcn,useFullNewtonStep)


[np,n]  = size(xbub);
lb      = lu(1,:);
ub      = lu(2,:);
mu      = size(xPar,1);


%
% Calculate approx. Hessian
%
% If flag is false, then the Hessian calculation is inaccurate. 
% Re-calculate using more individuals by increasing the number of 
% individuals used in its calculation. If the Hessian calculation is
% incorrect for whatever reason, then don't execute newton step. 
% If the Hessian is non-singular but newton direction has nan values 
% (due to some zeros is the gradient) then don't execute newton's step. 
%
flag    = false;
nInd    = (n^2+n)/2 + 1;
if mu > nInd
   mu = nInd; 
end
nretry      = 4;
retry       = 1;
xincParamBy = linspace(0,2,nretry);
nIndUsed    = nInd;
while ~flag && retry<nretry   
    lastwarn('')
    [H,Grad,isSingular,er,pd]  = calc_hessian_taylor(xbub(1:nIndUsed,:),fbub(1:nIndUsed),fcn);
    [~, warnID] = lastwarn;
    if isSingular || contains(warnID,'Matrix')
        nIndUsed  = round(nInd*exp(xincParamBy(retry)));
        retry     = retry + 1;
    else
        newtdir = (H\Grad(:,1:mu))';
        Q       = norm(newtdir(1,:));
        if ~any(isnan(Q))
           flag = true;  
        else 
            retry = nretry;
        end     
    end
    if nIndUsed > np || retry >= nretry
       xnewt                   = [];
       stats.nnewt             = 0;
       stats.nBlocks           = 0;
       stats.nIndUsed          = nIndUsed;
       stats.nRetryHessian     = retry;           
       stats.eigratio          = nan;
       stats.pd                = false;  
       stats.gradStepLen       = nan;
       stats.betamax           = nan;
       stats.nCorrectedBounds  = nan;
       stats.reason4failure    = warnID;
       disp(['Newton step failed. nIndUsed:', num2str(nIndUsed)])
       return;
    end
end

disp(['# of individuals used in successful Newton step is: ', num2str(nIndUsed)])

% Get Newton step direction
if ~useFullNewtonStep
    newtdir  = newtdir./vecnorm(newtdir')';     % get direction only
end


% Ready individuals for newton's step
nBlocks     = ceil(nnewt/mu);
xPar        = repmat(xPar,nBlocks,1);
xPar        = xPar(1:nnewt,:);
etaPar      = repmat(etaPar,nBlocks,1);
etaPar      = etaPar(1:nnewt,:);
newtdir     = repmat(newtdir,nBlocks,1);
newtdir     = newtdir(1:nnewt,:);


% Calculate: At what distance in the newton step direction do the parents
% hit the boundary
llimit   = (xPar - repmat(lb,nnewt,1))./(newtdir);
ulimit   = (xPar - repmat(ub,nnewt,1))./(newtdir);
betamax  = zeros(nnewt,1);
for i=1:nnewt
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
% NEWTON'S STEP
%
betarand   = rand(nnewt,1);
if nnewt>1
    xnewt = xPar - betarand.*newtdir;
else
    xnewt = xPar - newtdir;
end
if any(v)
    xnewt(v,:) = xPar(v,:) - betarand(v).*betamax(v).*newtdir(v,:);
end
nnewt      = size(xnewt,1);


% Newton's step stats
stats.nnewt             = nnewt;
stats.nBlocks           = nBlocks;
stats.nIndUsed          = nIndUsed;
stats.nRetryHessian     = retry - 1;
stats.eigratio          = er;
stats.pd                = pd;
stats.gradStepLen       = Q;
stats.betamax           = betamax;
stats.nCorrectedBounds  = length(find(v));
stats.reason4failure    = '';
end








