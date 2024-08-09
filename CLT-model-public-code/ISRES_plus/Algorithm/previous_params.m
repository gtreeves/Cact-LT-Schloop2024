function [X,F,Eta,D] = previous_params(X,F,Phi,Eta,xPoint,sortByError)
% This function collects parameters from the all previous generations that
% are within a close distance from the centroid.

if ~exist('sortByError','var')
    sortByError = false;
end

nParams     = size(X,2);
nPenalty    = size(Phi,2);
X           = permute(X,[1,3,2]);
X           = reshape(X,[],nParams); 
Eta         = permute(Eta,[1,3,2]);
Eta         = reshape(Eta,[],nParams); 
F           = reshape(F,[],1);
Phi         = permute(Phi,[1,3,2]);
Phi         = reshape(Phi,[],nPenalty); 
Phi(Phi<0)  = 0;
Phi         = sum(Phi.^2,2);  
v           = find(Phi <= 0);
X           = X(v,:);
Eta         = Eta(v,:);
F           = F(v);
v           = ~isnan(F);
X           = X(v,:);
Eta         = Eta(v,:);
F           = F(v);



% Calculate distance
D   = sum((X-xPoint).^2,2).^0.5;


% v      = D == 0;
% D(v)   = [];
% X(v,:) = [];
% F(v)   = [];


% Sort
if sortByError
    [F,i0] = sort(F);
    D      = D(i0);
else
    [D,i0] = sort(D);
    F      = F(i0);
end
X      = X(i0,:);
Eta    = Eta(i0,:);


