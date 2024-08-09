function dlambdaTdt = dlambdaT_by_dt(t,lambda,M,dudY,pp,soln,fhandle,...
            params,Vn,Vc,An,Am,stage,sumD2,sumDC,beta,nlen,tspan,Yobs,P)

% get y and yobs
y     = deval(soln,t);

[~,i0] = min(abs(tspan - t));
if i0 == length(tspan)
    i1 = [i0-1,i0];
elseif i0 == 1
    i1 = [i0, i0+1];
else
    i1 = [i0-1, i0, i0+1];
end
yobs = interp1(tspan(i1),Yobs(:,i1)',t,'spline');   

    
    
%yobs  = interp1(tspan,Yobs',t,'spline');
%yobs  = ppval(pp,t)';
% dyobs = ppval(dpp,t)';


% get total Dl
dltot = y(1:M) + y(2*M+1:3*M);
dltot = dltot';
 

% calculate djdy
T1 = sumD2*(yobs*dudY);
T2 = sumDC*(2*dltot*dudY);
dbetadY = (T1 - T2)./sumD2^2;

djdy  = -2 * (yobs-beta*dltot) *(beta*dudY + dltot'*dbetadY);


% Get lambda*jacobian by element-wise multiplication
lamjac   = fhandle(y,params,Vn,Vc,An,Am,stage,M,lambda);

% jac     = jac_dsat2(t,y,params,[],P,Vn,Vc,An,Am,stage);
% lamjac = lambda'*jac;
% 
% a = norm(lamjac - lamjac)


dlambdaTdt = (djdy - lamjac)';

end










