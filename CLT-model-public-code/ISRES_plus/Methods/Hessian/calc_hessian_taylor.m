function [H,Fdot,isSingular,eigratio,pd,flag] = calc_hessian_taylor(xs,fs,fcn)

warning('on','verbose')
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:singularMatrix

[nParams, np] = size(xs);


%
% If a single parameter-set is in input, then create a cluster around xs
% and run the error function on all sets
%
if nParams == 1
    nNeed = (np*np-np)/2 + np + 1;
    xpert = xs + (-1).^randi(2,nNeed,1).* randrange(0.1,0.2,nNeed).*xs;
    xs    = [xs;xpert];
    if ~exist('fs','var')
        fs = fcn(xs); 
    else
        fpert = fcn(xs(2:end,:));
        fs    = [fs;fpert];
    end
end

if ~exist('fs','var')
   fs = fcn(xs); 
end

%xs = 10.^xs;



%
% Start constructing the Hessian matrix
%

% 1. Get xixj and xi parts of the vector
Mat             = zeros(nParams,np*(np+1)/2 + np);
O               = ones(np);
O(1:np+1:np*np) = 0.5;
for i = 1:nParams
    x1          = xs(i,:);
    x2          = (x1'*x1).*O;
    x2          = triu(x2)';
    x2          = x2(:);
    x2          = x2(x2~=0);
    Mat(i,:)    = [x2;x1'];
end                                         



% 2. add to that the ones and evaluate it!
Mat  = [Mat,ones(nParams,1)];
vec  = Mat\fs;



% 3. get hessian and eigenvalues of the matrix
H           = zeros(np,np);
isSingular  = false; 
ibegin      = 1;
for i = 1:np
   iend       = ibegin + np - i; 
   H(i,i:end) = vec(ibegin:iend);   
   ibegin     = iend + 1;
end
H       = (H+H') - eye(np).*diag(H);
a       = vec(end-np:end-1);
if isnan(norm(H))
   Fdot        = zeros(np,nParams);
   isSingular  = true;
   eigratio    = inf;
   pd          = false;
   flag        = false;
   return
end
e       = eig(H);
if any(e==0)
   isSingular = true; 
end
if all(e > 0)
   pd   = true;
else
   pd   = false;
end
eigratio = max(abs(e))/min(abs(e));



%
% Now get the value of the Gradient from the vector "a"
%
Fdot        = zeros(np,nParams);
for i = 1:nParams
   x1           = xs(i,:);   
   Fdot(:,i)    = a + sum(H.*x1,2);   
end



%
% check to see if calculation is correct by comparing the values of the
% constant. Flag must be true for a successful calculation. 
%
for i = 1:nParams
    x1      = xs(i,:);   
    t2      = sum(H.*x1,2);
    t2      = sum(t2.*x1');
    T2(i)   = t2/2;
    fdot    = Fdot(:,i);
    fdot    = sum(fdot.*x1');
    T3(i)   = fdot;
end
ccal    = fs + T2' - T3';
ccal    = mean(ccal);
flag    = false;
if abs(ccal - vec(end)) < 1e-4
    flag = true;
end

end

